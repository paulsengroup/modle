#include <absl/hash/hash.h>
#include <absl/time/clock.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

#include "modle/config.hpp"
#include "modle/cu/simulation.hpp"
#include "modle/dna.hpp"
#include "modle/extrusion_barriers.hpp"
#include "modle/libmodle_io.hpp"

namespace modle::cu {

Simulation::Simulation(const modle::Config& config, size_t grid_size, size_t block_size,
                       size_t max_grid_size, size_t device_heap_size)
    : modle::Config(config),
      _config(&config),
      _genome(path_to_chrom_sizes, path_to_extr_barriers, path_to_chrom_subranges,
              write_contacts_for_ko_chroms),
      _grid_size(grid_size),
      _block_size(block_size) {
  assert(grid_size != 0);   // NOLINT
  assert(block_size != 0);  // NOLINT
  if (max_grid_size == std::numeric_limits<size_t>::max()) {
    max_grid_size = grid_size;
  }
  cuda::device::current::get().set_resource_limit(cudaLimitMallocHeapSize, device_heap_size);
  this->_global_state_host = std::make_unique<GlobalStateHost>(
      this->_grid_size, this->_block_size, max_grid_size, this->_config->to_cuda_config(),
      this->_genome.longest_chromosome().nlefs(this->_config->number_of_lefs_per_mbp),
      this->_genome.chromosome_with_max_nbarriers().num_valid_barriers());
}

/*
void Simulation::init_global_state(size_t grid_size, size_t block_size, size_t max_nlefs,
                                   size_t max_nbarriers, cuda::device_t dev) {
  this->_grid_size = grid_size;
  this->_block_size = block_size;
  this->_global_state_host = std::make_unique<GlobalStateHost>(
      grid_size, block_size, this->_config->to_cuda_config(), max_nlefs, max_nbarriers, dev);
}

void Simulation::update_tasks(std::vector<Task>&& new_tasks) noexcept {
  this->_tasks = std::move(new_tasks);
  this->_grid_size = this->_tasks.size();

  this->_global_state_host->write_tasks_to_device(this->_tasks);
}

void Simulation::update_barriers(std::vector<uint32_t>&& barrier_pos,
                                 std::vector<dna::Direction>&& barrier_dir) noexcept {
  assert(barrier_pos.size() == barrier_dir.size());  // NOLINT
  this->_barrier_positions = std::move(barrier_pos);
  this->_barrier_directions = std::move(barrier_dir);

  this->_global_state_host->write_barriers_to_device(this->_barrier_positions,
                                                     this->_barrier_directions);
}

void Simulation::run_batch() {
  assert(this->_grid_size == this->_tasks.size() && this->_grid_size != 0);  // NOLINT
  // TODO handle case where nlefs and nbarriers can be different across states

  this->_global_state_host->init(this->_grid_size, this->_block_size);
  this->_global_state_host->setup_burnin_phase(this->_tasks);

  // const auto shared_memory_size = std::max(nlefs, nbarriers) * sizeof(uint32_t);

  bool end_of_simulation{false};

  for (this->_current_epoch = 0UL; !end_of_simulation; ++this->_current_epoch) {
    end_of_simulation = this->simulate_one_epoch();
    if (this->_current_epoch == 10) {
      return;
    }
  }
}
 */

void Simulation::run() {
  std::atomic<bool> end_of_simulation = false;
  std::mutex progress_queue_mutex;  // Protect rw access to progress_queue
  std::deque<std::pair<Chromosome*, size_t>> progress_queue;

  // These variables are here just for compatibility with the function used to write contacts to
  // disk
  std::mutex barrier_mutex;
  absl::flat_hash_map<Chromosome*, std::unique_ptr<std::vector<ExtrusionBarrier>>> barriers;
  /*
    const auto genome =
        this->instantiate_genome(this->path_to_chrom_sizes, this->path_to_extr_barriers,
                                 this->path_to_chrom_subranges, this->write_contacts_for_ko_chroms);
                                 */

  const auto& longest_chrom = this->_genome.chromosome_with_longest_name().name().size();
  const auto max_nbarriers =
      static_cast<uint32_t>(this->_genome.chromosome_with_max_nbarriers().num_valid_barriers());
  const auto max_nlefs = static_cast<uint32_t>(
      std::round(this->number_of_lefs_per_mbp *
                 (static_cast<double>(this->_genome.longest_chromosome().size()) / 1.0e6)));

  const auto nthreads_ = 208'896;  // 48 * 32 * 68UL;
  const auto block_size = 192UL;   // 128U;
  const auto grid_size = nthreads_ / block_size;
  const auto max_grid_size = grid_size;
  const auto task_batch_size = grid_size;

  auto device = cuda::device::current::get();
  auto global_state_host =
      modle::cu::GlobalStateHost(grid_size, block_size, max_grid_size,
                                 this->_config->to_cuda_config(), max_nlefs, max_nbarriers, device);

  std::vector<Chromosome*> chroms(static_cast<size_t>(this->_genome.size()));
  std::transform(this->_genome.begin(), this->_genome.end(), chroms.begin(),
                 [](auto& chrom) { return &chrom; });

  std::vector<uint32_t> barrier_pos_buff(chroms.front()->nbarriers());
  std::vector<modle::cu::dna::Direction> barrier_direction_buff(chroms.front()->nbarriers());

  for (const auto& chrom : chroms) {
    if (!chrom->ok()) {
      fmt::print(stderr, FMT_STRING("SKIPPING '{}'...\n"), chrom->name());
      if (this->write_contacts_for_ko_chroms) {
        chrom->allocate_contacts(this->bin_size,
                                 this->diagonal_width);  // TODO: Can we remove this?
        std::scoped_lock l(barrier_mutex, progress_queue_mutex);
        progress_queue.emplace_back(chrom, ncells);
        barriers.emplace(chrom, nullptr);
      }
      continue;
    }

    const auto t0 = absl::Now();
    // Allocate the contact matrix. Once it's not needed anymore, the contact matrix will be
    // de-allocated by the thread that is writing contacts to disk
    chrom->allocate_contacts(this->bin_size, this->diagonal_width);

    modle::cu::Simulation::Task base_task{};
    base_task.chrom_id = chrom->id();
    base_task.chrom_start = chrom->start_pos();
    base_task.chrom_end = chrom->end_pos();

    // Compute # of LEFs to be simulated based on chrom. sizes
    base_task.nlefs = chrom->nlefs(this->number_of_lefs_per_mbp);

    base_task.n_target_contacts = 0;
    if (this->target_contact_density != 0) {  // Compute the number of simulation rounds required
      // to reach the target contact density
      base_task.n_target_contacts = static_cast<uint32_t>(
          std::max(1.0, std::round((this->target_contact_density *
                                    static_cast<double>(chrom->contacts().npixels())) /
                                   static_cast<double>(this->ncells))));
    }
    base_task.n_target_epochs = base_task.n_target_contacts == 0
                                    ? static_cast<uint32_t>(this->simulation_iterations)
                                    : std::numeric_limits<uint32_t>::max();

    base_task.nbarriers = chrom->num_valid_barriers();

    this->_tasks.resize(task_batch_size);
    size_t cellid = 0;
    const auto nbatches = (this->ncells + task_batch_size - 1) / task_batch_size;
    fmt::print(stderr, FMT_STRING("Processing '{}' in {} batches of up to {} blocks...\n"),
               chrom->name(), nbatches, grid_size);
    for (auto batchid = 0UL; batchid < nbatches; ++batchid) {
      // Generate a batch of tasks for all the simulations involving the current chrom
      this->_tasks.resize(std::min(this->ncells - cellid, task_batch_size));
      if (this->_tasks.size() == 0) {
        break;
      }
      std::generate_n(this->_tasks.begin(), this->_tasks.size(), [&]() {
        auto t = base_task;
        t.cell_id = cellid++;
        t.seed = absl::Hash<uint64_t>{}(this->seed) +
                 absl::Hash<std::string_view>{}(chrom->name()) +
                 absl::Hash<size_t>{}(chrom->size()) + absl::Hash<size_t>{}(t.cell_id);
        return t;
      });

      modle::cu::Simulation::run_batch(this->_tasks, this->_barrier_positions,
                                       this->_barrier_directions);
    }
    fmt::print(stderr, FMT_STRING("Done simulating '{}'. Simulation took {}.\n"), chrom->name(),
               absl::FormatDuration(absl::Now() - t0));
  }
}
}  // namespace modle::cu
