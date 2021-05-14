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
      _block_size(block_size),
      _global_state_host(
          this->_grid_size, this->_block_size,
          max_grid_size == std::numeric_limits<size_t>::max() ? this->_grid_size : max_grid_size,
          device_heap_size, this->_config->to_cuda_config(),
          this->_genome.longest_chromosome().nlefs(this->_config->number_of_lefs_per_mbp),
          this->_genome.chromosome_with_max_nbarriers().num_valid_barriers()) {
  assert(grid_size != 0);   // NOLINT
  assert(block_size != 0);  // NOLINT
}

void Simulation::run() {
  std::mutex progress_queue_mutex;  // Protect rw access to progress_queue
  std::deque<std::pair<Chromosome*, size_t>> progress_queue;

  const auto task_batch_size = this->_grid_size;

  for (auto& chrom : this->_genome) {
    if (!chrom.ok()) {
      fmt::print(stderr, FMT_STRING("SKIPPING '{}'...\n"), chrom.name());
      if (this->write_contacts_for_ko_chroms) {
        chrom.allocate_contacts(this->bin_size,
                                this->diagonal_width);  // TODO: Can we remove this?
        std::scoped_lock l(progress_queue_mutex);
        progress_queue.emplace_back(&chrom, ncells);
      }
      continue;
    }

    const auto t0 = absl::Now();
    // Allocate the contact matrix. Once it's not needed anymore, the contact matrix will be
    // de-allocated by the thread that is writing contacts to disk
    chrom.allocate_contacts(this->bin_size, this->diagonal_width);

    modle::cu::Task base_task{};
    base_task.chrom_id = chrom.id();
    base_task.chrom_start = chrom.start_pos();
    base_task.chrom_end = chrom.end_pos();

    // Compute # of LEFs to be simulated based on chrom. sizes
    base_task.nlefs = chrom.nlefs(this->number_of_lefs_per_mbp);

    base_task.n_target_contacts = 0;
    if (this->target_contact_density != 0) {  // Compute the number of simulation rounds required
      // to reach the target contact density
      base_task.n_target_contacts = static_cast<uint32_t>(
          std::max(1.0, std::round((this->target_contact_density *
                                    static_cast<double>(chrom.contacts().npixels())) /
                                   static_cast<double>(this->ncells))));
    }
    base_task.n_target_epochs = base_task.n_target_contacts == 0
                                    ? static_cast<uint32_t>(this->simulation_iterations)
                                    : std::numeric_limits<uint32_t>::max();

    const auto barriers = this->_genome.generate_vect_of_barriers(
        chrom.name(), ctcf_occupied_self_prob, ctcf_not_occupied_self_prob);
    this->_barrier_positions.resize(barriers.size());
    this->_barrier_directions.resize(barriers.size());

    base_task.nbarriers = barriers.size();
    std::transform(barriers.begin(), barriers.end(), this->_barrier_positions.begin(),
                   [](const auto& b) { return b.pos(); });
    std::transform(barriers.begin(), barriers.end(), this->_barrier_directions.begin(),
                   [](const auto& b) {
                     return b.blocking_direction_major() == modle::dna::Direction::fwd
                                ? cu::dna::Direction::rev
                                : cu::dna::Direction::fwd;
                   });

    auto cellid = 0UL;
    const auto nbatches = (this->ncells + task_batch_size - 1) / task_batch_size;
    fmt::print(stderr, FMT_STRING("Processing '{}' in {} batches of up to {} blocks...\n"),
               chrom.name(), nbatches, this->_grid_size);
    for (auto batchid = 0UL; batchid < nbatches; ++batchid) {
      // Generate a batch of tasks for all the simulations involving the current
      // chrom
      this->_tasks.resize(std::min(this->ncells - cellid, task_batch_size));
      if (this->_tasks.size() == 0) {
        break;
      }
      std::generate_n(this->_tasks.begin(), this->_tasks.size(), [&]() {
        auto t = base_task;
        t.cell_id = cellid++;
        t.seed = absl::Hash<uint64_t>{}(this->seed) + absl::Hash<std::string_view>{}(chrom.name()) +
                 absl::Hash<size_t>{}(chrom.size()) + absl::Hash<size_t>{}(t.cell_id);
        return t;
      });

      modle::cu::Simulation::run_batch(this->_tasks, this->_barrier_positions,
                                       this->_barrier_directions);
    }
    fmt::print(stderr, FMT_STRING("Done simulating '{}'. Simulation took {}.\n"), chrom.name(),
               absl::FormatDuration(absl::Now() - t0));
  }
}
}  // namespace modle::cu
