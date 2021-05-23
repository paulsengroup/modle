#include <absl/time/clock.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cuda/runtime_api.hpp>
#include <mutex>
#include <numeric>
#include <string>
#include <vector>

#include "modle/config.hpp"
#include "modle/cooler.hpp"
#include "modle/cu/simulation.hpp"
#include "modle/dna.hpp"
#include "modle/extrusion_barriers.hpp"
#include "modle/setup.hpp"

namespace modle::cu {

Simulation::Simulation(const modle::Config& config_, size_t grid_size, size_t block_size,
                       size_t max_grid_size, size_t device_heap_size)
    : modle::Config(config_),
      _config(&config_),
      _genome(path_to_chrom_sizes, path_to_extr_barriers, path_to_chrom_subranges,
              write_contacts_for_ko_chroms),
      _grid_size(grid_size),
      _block_size(block_size),
      _global_state_host(
          _grid_size, _block_size,
          max_grid_size == std::numeric_limits<size_t>::max() ? _grid_size : max_grid_size,
          device_heap_size, _config->to_cuda_config(),
          _genome.longest_chromosome().nlefs(_config->number_of_lefs_per_mbp),
          _genome.chromosome_with_max_nbarriers().num_valid_barriers(),
          _genome.max_target_contacts(bin_size, diagonal_width, target_contact_density,
                                      simulation_iterations, lef_fraction_contact_sampling,
                                      number_of_lefs_per_mbp, ncells)) {
  assert(grid_size != 0);   // NOLINT
  assert(block_size != 0);  // NOLINT
}

void Simulation::run() {
  std::atomic<bool> end_of_simulation{false};
  std::mutex progress_queue_mutex;  // Protect rw access to progress_queue
  std::deque<std::pair<Chromosome*, size_t>> progress_queue;

  std::thread writer_thread([&]() {  // This thread is in charge of writing contacts to disk
    this->write_contacts_to_disk(progress_queue, progress_queue_mutex, end_of_simulation);
  });

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

    base_task.nbarriers = barriers.size();

    this->_barrier_positions.resize(barriers.size());
    this->_barrier_directions.resize(barriers.size());
    this->_barrier_probs_occ_to_occ.resize(barriers.size());
    this->_barrier_probs_nocc_to_nocc.resize(barriers.size());

    for (auto i = 0UL; i < barriers.size(); ++i) {
      const auto& b = barriers[i];
      this->_barrier_positions[i] = b.pos();
      this->_barrier_directions[i] = b.blocking_direction_major() == modle::dna::Direction::fwd
                                         ? modle::cu::dna::Direction::fwd
                                         : modle::cu::dna::Direction::rev;
      this->_barrier_probs_occ_to_occ[i] = static_cast<float>(b.prob_occupied_to_occupied());
      this->_barrier_probs_nocc_to_nocc[i] =
          static_cast<float>(b.prob_not_occupied_to_not_occupied());
    }

    auto cellid = 0UL;
    const auto nbatches = (this->ncells + task_batch_size - 1) / task_batch_size;
    fmt::print(stderr, FMT_STRING("Processing '{}' in {} batches of up to {} blocks...\n"),
               chrom.name(), nbatches, this->_grid_size);
    for (auto batchid = 0UL; batchid < nbatches; ++batchid) {
      // Generate a batch of tasks for all the simulations involving the current
      // chrom
      this->_tasks.resize(std::min(this->ncells - cellid, task_batch_size));
      if (this->_tasks.empty()) {
        break;
      }
      std::generate_n(this->_tasks.begin(), this->_tasks.size(), [&]() {
        auto t = base_task;
        t.cell_id = cellid++;
        t.seed = chrom.hash(this->seed, t.cell_id);
        return t;
      });

      modle::cu::Simulation::run_batch(this->_tasks, this->_barrier_positions,
                                       this->_barrier_directions, this->_barrier_probs_occ_to_occ,
                                       this->_barrier_probs_nocc_to_nocc);

      modle::cu::Simulation::update_contacts_for_chrom(chrom);
    }
    std::scoped_lock l(progress_queue_mutex);
    progress_queue.emplace_back(&chrom, this->ncells);
    fmt::print(stderr, FMT_STRING("Done simulating '{}'. Simulation took {}.\n"), chrom.name(),
               absl::FormatDuration(absl::Now() - t0));
  }
  {
    std::scoped_lock l(progress_queue_mutex);
    progress_queue.emplace_back(nullptr, 0);
  }
  writer_thread.join();
}

void Simulation::write_contacts_to_disk(
    std::deque<std::pair<Chromosome*, size_t>>& progress_queue, std::mutex& progress_queue_mutex,
    std::atomic<bool>& end_of_simulation) {  // This thread is in charge of writing contacts to disk
  Chromosome* chrom_to_be_written = nullptr;
  const auto max_str_length =
      std::max_element(  // Find chrom with the longest name
          this->_genome.begin(), this->_genome.end(),
          [](const auto& c1, const auto& c2) { return c1.name().size() < c2.name().size(); })
          ->name()
          .size();

  auto c =
      this->path_to_output_file.empty()
          ? nullptr
          : std::make_unique<cooler::Cooler>(this->path_to_output_file, cooler::Cooler::WRITE_ONLY,
                                             this->bin_size, max_str_length);

  auto sleep_us = 100;
  while (true) {  // Structuring the loop in this way allows us to sleep without holding the mutex
    sleep_us = std::min(500000, sleep_us * 2);  // NOLINT
    std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
    {
      std::scoped_lock l(progress_queue_mutex);
      if (progress_queue.empty()) {
        // There are no contacts to write to disk at the moment. Go back to sleep
        continue;
      }

      // chrom == nullptr is the end-of-queue signal
      if (auto& [chrom, count] = progress_queue.front(); chrom == nullptr) MODLE_UNLIKELY {
          end_of_simulation = true;
          return;
        }
      // count == ncells signals that we are done simulating the current chromosome
      else if (count == ncells) {
        chrom_to_be_written = chrom;
        progress_queue.pop_front();
      } else
        MODLE_LIKELY {
          assert(count < ncells);  // NOLINT
          continue;
        }
    }
    sleep_us = 100;
    try {
      if (c) {  // c == nullptr only when --skip-output is used
        fmt::print(stderr, "Writing contacts for '{}' to file {}...\n", chrom_to_be_written->name(),
                   c->get_path());
        c->write_or_append_cmatrix_to_file(
            chrom_to_be_written->contacts(), chrom_to_be_written->name(),
            chrom_to_be_written->start_pos(), chrom_to_be_written->end_pos(),
            chrom_to_be_written->size(), true);
        fmt::print(stderr, "Written {} contacts for '{}' in {:.2f}M pixels to file {}.\n",
                   chrom_to_be_written->contacts().get_tot_contacts(), chrom_to_be_written->name(),
                   static_cast<double>(chrom_to_be_written->contacts().npixels()) / 1.0e6,
                   c->get_path());
      }
      // Deallocate the contact matrix to free up unused memory
      chrom_to_be_written->deallocate_contacts();
    } catch (const std::runtime_error& err) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("The following error occurred while writing contacts for '{}' to file {}: {}"),
          chrom_to_be_written->name(), c->get_path(), err.what()));
    }
  }
}

void Simulation::update_contacts_for_chrom(Chromosome& chrom) {
  BlockState block_state;
  std::vector<uint2> contacts;
  const auto global_state = this->_global_state_host.get_copy_of_device_instance();
  for (auto i = 0UL; i < this->_grid_size; ++i) {
    cuda::memory::copy_single(&block_state, global_state.block_states + i);
    contacts.resize(block_state.contact_local_buff_size);
    cuda::memory::copy(contacts.data(), block_state.contact_local_buff,
                       block_state.contact_local_buff_size * sizeof(uint2));

    for (const auto& [row, col] : contacts) {
      chrom.increment_contacts(row, col);
    }
  }
}

}  // namespace modle::cu
