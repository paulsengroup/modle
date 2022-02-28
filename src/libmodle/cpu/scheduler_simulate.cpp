// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "modle/simulation.hpp"
// clang-format on

#include <absl/container/btree_map.h>    // for btree_iterator
#include <absl/container/fixed_array.h>  // for FixedArray
#include <absl/types/span.h>             // for MakeConstSpan, Span
#include <fmt/format.h>                  // for make_format_args, vformat_to, FMT_STRING
#include <fmt/ostream.h>
#include <moodycamel/blockingconcurrentqueue.h>  // for BlockingConcurrentQueue
#include <moodycamel/concurrentqueue.h>          // for ConsumerToken, ProducerToken
#include <spdlog/spdlog.h>                       // for info

#include <algorithm>                        // for max, copy, min, find_if, generate
#include <atomic>                           // for atomic
#include <boost/filesystem/operations.hpp>  // for exists, remove
#include <boost/filesystem/path.hpp>        // for path
#include <cassert>                          // for assert
#include <chrono>                           // for microseconds, milliseconds
#include <cmath>                            // for round
#include <deque>                            // for deque, operator-, operator!=, _Deque_ite...
#include <exception>                        // for exception_ptr, exception, current_exception
#include <iterator>                         // for move_iterator, make_move_iterator
#include <limits>                           // for numeric_limits
#include <mutex>                            // for mutex, scoped_lock
#include <stdexcept>                        // for runtime_error
#include <thread>                           // IWYU pragma: keep for sleep_for
#include <thread_pool/thread_pool.hpp>      // for thread_pool
#include <utility>                          // for pair
#include <vector>                           // for vector

#include "modle/common/common.hpp"                      // for u64
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...
#include "modle/genome.hpp"                             // for Chromosome, Genome
#include "modle/interval_tree.hpp"                      // for IITree, IITree::data

namespace modle {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Simulation::run_simulate() {
  if (!this->skip_output) {  // Write simulation params to file
    assert(boost::filesystem::exists(this->path_to_output_prefix.parent_path()));
    if (this->force) {
      boost::filesystem::remove(this->path_to_output_file_cool);
    }
  }
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SHORTEN_64_TO_32
  this->_tpool.reset(this->nthreads + 1);
  DISABLE_WARNING_POP

  // These are the threads spawned by run_simulate:
  // - 1 thread to write contacts to disk. This thread pops Chromosome* from a std::deque once the
  //   simulation on the Chromosome* has been ultimated. This thread is also responsible of
  //   freeing memory as soon as it is not needed by any other thread.
  // - The main thread (i.e. this thread), which loops over chromosomes, instantiates data
  //   structures that are shared across simulation threads (e.g. the vector of extr. barriers),
  //   then it submits simulation tasks to a concurrent queue
  // - N simulation threads. These threads pop tasks from the concurrent queue fed by the main
  //   thread and carry out the actual simulation
  //
  // Switching to this thread-layout allowed us to significantly reduce memory allocations, as by
  // running few simulation threads for the whole simulation allows us to recycle buffers more
  // efficiently. Most of the time, the only data moved are tasks, which are implemented as
  // light-weight structs

  std::mutex progress_queue_mutex;  // Protect rw access to progress_queue
  std::deque<std::pair<Chromosome*, usize>> progress_queue;

  const auto task_batch_size_deq = [this]() -> usize {
    if (this->num_cells <= this->nthreads) {
      return 1;
    }
    return std::min(usize(16), this->num_cells / this->nthreads);
  }();

  const usize task_batch_size_enq = 2 * task_batch_size_deq;
  const usize queue_capacity = 2 * this->nthreads * task_batch_size_deq;

  // Queue used to submit simulation tasks to the thread pool
  moodycamel::BlockingConcurrentQueue<Simulation::Task> task_queue(queue_capacity, 1, 0);
  moodycamel::ProducerToken ptok(task_queue);

  std::mutex model_state_logger_mtx;  // Protect rw access to the log file located at
                                      // this->path_to_model_state_log_file
  if (this->log_model_internal_state && !this->skip_output) {
    compressed_io::Writer(this->path_to_model_state_log_file).write(model_state_log_header);
  }

  try {
    this->_tpool.push_task([&]() {  // This thread is in charge of writing contacts to disk
      this->write_contacts_to_disk(progress_queue, progress_queue_mutex);
    });

    for (u64 tid = 0; tid < this->nthreads; ++tid) {  // Start simulation threads
      this->_tpool.push_task([&, tid]() {
        this->simulate_worker(tid, task_queue, progress_queue, progress_queue_mutex,
                              model_state_logger_mtx, task_batch_size_deq);
      });
    }

    // The remaining code submits simulation tasks to the queue. Then it waits until all the tasks
    // have been completed and contacts have been written to disk

    absl::FixedArray<Task> tasks(task_batch_size_enq);
    usize taskid = 0;

    // Loop over chromosomes
    for (auto& chrom : this->_genome) {
      if (!this->ok()) {
        this->handle_exceptions();
      }
      // Don't simulate KO chroms (but write them to disk if the user desires so)
      if (!chrom.ok() || chrom.num_barriers() == 0) {
        spdlog::info(FMT_STRING("SKIPPING \"{}\"..."), chrom.name());
        if (this->write_contacts_for_ko_chroms) {
          std::scoped_lock lck(progress_queue_mutex);
          progress_queue.emplace_back(&chrom, num_cells);
        }
        continue;
      }

      {
        std::scoped_lock lck(progress_queue_mutex);

        // Signal that we have started processing the current chrom
        progress_queue.emplace_back(&chrom, usize(0));
      }

      // Compute # of LEFs to be simulated based on chrom.sizes
      const auto nlefs = this->compute_num_lefs(chrom.simulated_size());

      const auto npixels = chrom.npixels(this->diagonal_width, this->bin_size);
      const auto target_contacts =
          (this->compute_tot_target_contacts(npixels) + this->num_cells - 1) / this->num_cells;
      const auto target_epochs = this->compute_tot_target_epochs();

      usize cellid = 0;
      const auto nbatches = (this->num_cells + task_batch_size_enq - 1) / task_batch_size_enq;
      for (usize batchid = 0; batchid < nbatches; ++batchid) {
        // Generate a batch of tasks for the current chrom
        std::generate(tasks.begin(), tasks.end(), [&]() {
          return Task{{taskid++, &chrom, cellid++, target_epochs, target_contacts, nlefs,
                       chrom.barriers().data()}};
        });
        const auto ntasks =
            cellid > this->num_cells ? tasks.size() - (cellid - this->num_cells) : tasks.size();
        auto sleep_us = 100;
        while (  // Enqueue tasks
            !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), ntasks)) {
          if (!this->ok()) {
            this->handle_exceptions();
          }
          sleep_us = std::min(100000, sleep_us * 2);
          std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
        }
      }
    }

    {  // Signal end of simulation to the thread that is writing contacts to disk
      std::scoped_lock lck(progress_queue_mutex);
      progress_queue.emplace_back(nullptr, usize(0));
    }
    this->_tpool.wait_for_tasks();
    assert(this->_end_of_simulation);
    assert(!this->_exception_thrown);
  } catch (...) {
    this->_exception_thrown = true;
    this->_tpool.paused = true;
    this->_tpool.wait_for_tasks();
    throw;
  }
}

void Simulation::simulate_worker(const u64 tid,
                                 moodycamel::BlockingConcurrentQueue<Simulation::Task>& task_queue,
                                 std::deque<std::pair<Chromosome*, usize>>& progress_queue,
                                 std::mutex& progress_queue_mtx, std::mutex& model_state_logger_mtx,
                                 const usize task_batch_size) {
  spdlog::info(FMT_STRING("Spawning simulation thread {}..."), tid);

  moodycamel::ConsumerToken ctok(task_queue);
  absl::FixedArray<Task> task_buff(task_batch_size);  // Tasks are dequeued in batch.
  // This is to reduce contention
  // when accessing the queue

  // This state object owns all the buffers and PRNG + seed required in order to simulate loop
  // extrusion for a single cell. Buffers are allocated once and resized, cleared and reused
  // throughout the simulation
  Simulation::State local_state;

  const auto tmp_model_internal_state_log_path = [&]() {
    auto p = this->path_to_model_state_log_file;
    if (this->log_model_internal_state && !this->skip_output) {
      p.replace_extension(fmt::format(FMT_STRING("{}{}"), tid, p.extension().string()));
      local_state.model_state_logger = std::make_unique<compressed_io::Writer>(p);
    }
    return p;
  }();

  try {
    while (this->ok()) {
      const auto avail_tasks = this->consume_tasks_blocking(task_queue, ctok, task_buff);
      if (avail_tasks == 0) {
        assert(this->_end_of_simulation);
        // Reached end of simulation (i.e. all tasks have been processed)
        if (this->log_model_internal_state && !this->skip_output) {
          assert(!tmp_model_internal_state_log_path.empty());
          // Ensure all log records have been written to disk
          local_state.model_state_logger = nullptr;
          std::scoped_lock<std::mutex> lck(model_state_logger_mtx);
          utils::concatenate_files<true>(this->path_to_model_state_log_file,
                                         tmp_model_internal_state_log_path);
        }
        return;
      }

      // Loop over new tasks
      for (const auto& task : absl::MakeConstSpan(task_buff.data(), avail_tasks)) {
        if (!this->ok()) {
          return;
        }
        // Resize and reset buffers
        task.chrom->allocate_contacts(this->bin_size, this->diagonal_width);
        local_state = task;            // Set simulation state based on task data
        local_state.resize_buffers();  // This resizes buffers based on the nlefs to be simulated
        local_state.reset_buffers();   // Clear all buffers

        if (task.cell_id == 0) {
          // Print a status update when we are processing cell #0 for a given chromosome
          const auto tot_target_epochs = [&]() {
            if (this->target_contact_density == 0.0) {
              return this->compute_tot_target_epochs();
            }

            const auto tot_target_contacts =
                this->compute_tot_target_contacts(task.chrom->npixels());
            const auto new_contacts_per_epoch =
                static_cast<double>(task.num_lefs) * this->lef_fraction_contact_sampling *
                static_cast<double>(this->number_of_tad_contacts_per_sampling_event +
                                    this->number_of_loop_contacts_per_sampling_event);
            return std::max(usize(1), static_cast<usize>(static_cast<double>(tot_target_contacts) /
                                                         new_contacts_per_epoch));
          }();
          spdlog::info(
              FMT_STRING(
                  "Begin processing \"{}\": simulating ~{} epochs across {} cells using {} LEFs "
                  "(~{} epochs per cell)..."),
              task.chrom->name(), tot_target_epochs, this->num_cells, task.num_lefs,
              tot_target_epochs / this->num_cells);
        }

        // Start the simulation kernel
        Simulation::simulate_one_cell(local_state);

        // Update progress for the current chrom
        std::scoped_lock lck(progress_queue_mtx);
        auto progress = std::find_if(progress_queue.begin(), progress_queue.end(),
                                     [&](const auto& p) { return task.chrom == p.first; });
        assert(progress != progress_queue.end());

        if (++progress->second == num_cells && !this->_exception_thrown) {
          // We are done simulating loop-extrusion on task.chrom: print a status update
          spdlog::info(FMT_STRING("Simulation of \"{}\" successfully completed."),
                       task.chrom->name());
        }
      }
    }
  } catch (const std::exception& e) {
    std::scoped_lock<std::mutex> l(this->_exceptions_mutex);

    const auto excp = std::runtime_error(fmt::format(
        FMT_STRING("Detected an error in worker thread {}:\n{}\n{}"), tid, local_state, e.what()));
    this->_exceptions.emplace_back(std::make_exception_ptr(excp));
    this->_exception_thrown = true;

  } catch (...) {
    std::scoped_lock<std::mutex> l(this->_exceptions_mutex);
    this->_exceptions.emplace_back(std::current_exception());
    this->_exception_thrown = true;
  }
}

}  // namespace modle
