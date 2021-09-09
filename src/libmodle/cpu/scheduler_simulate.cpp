// clang-format off
#include "modle/simulation.hpp"
// clang-format on

#include <absl/container/fixed_array.h>          // for FixedArray
#include <absl/types/span.h>                     // for MakeConstSpan, Span
#include <fmt/format.h>                          // for make_format_args, vformat_to, format
#include <moodycamel/blockingconcurrentqueue.h>  // for BlockingConcurrentQueue
#include <moodycamel/concurrentqueue.h>          // for ConsumerToken, ProducerToken
#include <spdlog/spdlog.h>                       // for info

#include <algorithm>                        // for copy, max, min, find_if, generate
#include <atomic>                           // for atomic
#include <boost/filesystem/operations.hpp>  // for exists, remove
#include <boost/filesystem/path.hpp>        // for path
#include <cassert>                          // for assert
#include <chrono>                           // for microseconds, milliseconds
#include <cmath>                            // for round, floor, log
#include <cstddef>                          // for size_t
#include <cstdint>                          // for uint64_t
#include <deque>                            // for deque, operator-, operator!=, _Deque_ite...
#include <exception>                        // for exception_ptr, exception, current_exception
#include <iterator>                         // for move_iterator, make_move_iterator
#include <limits>                           // for numeric_limits
#include <mutex>                            // for mutex, scoped_lock
#include <stdexcept>                        // for runtime_error
#include <thread>                           // for sleep_for
#include <thread_pool/thread_pool.hpp>      // for thread_pool
#include <utility>                          // for pair
#include <vector>                           // for vector

#include "modle/genome.hpp"         // for Chromosome, Genome
#include "modle/interval_tree.hpp"  // for IITree, IITree::data

namespace modle {

void Simulation::run_simulate() {
  if (!this->skip_output) {  // Write simulation params to file
    assert(boost::filesystem::exists(this->path_to_output_prefix.parent_path()));  // NOLINT
    if (this->force) {
      boost::filesystem::remove(this->path_to_output_file_cool);
    }
  }

  this->_tpool.reset(this->nthreads + 1);

  // These are the threads spawned by simulate_extrusion:
  // - 1 thread to write contacts to disk. This thread pops Chromosome* from a std::deque once the
  //   simulation on the Chromosome* has been ultimated. This thread is also responsible of
  // freeing memory as soon as it is not needed by any other thread.
  // - The main thread (i.e. this thread), which loops over chromosomes, instantiates data
  //   structures that are shared across simulation threads (e.g. the vector of extr. barriers),
  // then it submits simulation tasks to a concurrent queue
  // - N simulation threads. These threads pop tasks from the concurrent queue fed by the main
  //   thread and carry out the actual simulation
  //
  // Switching to this thread-layout allowed us to significantly reduce memory allocations, as by
  // running few simulation threads for the whole simulation allows us to recycle buffers more
  // efficiently. Most of the time, the only data moved are tasks, which are structs consisting of
  // few numbers, a pointer and a Span (i.e. ptr + size)

  std::mutex progress_queue_mutex;  // Protect rw access to progress_queue
  std::deque<std::pair<Chromosome*, size_t>> progress_queue;

  const auto task_batch_size_deq = [this]() -> size_t {
    if (this->num_cells <= this->nthreads) {
      return 1;
    }
    return std::min(16UL, this->num_cells / this->nthreads);  // NOLINT
  }();
  const size_t task_batch_size_enq = 2 * task_batch_size_deq;  // NOLINTNEXTLINE
  const size_t queue_capacity = 2 * this->nthreads * task_batch_size_deq;
  // Queue used to submit simulation tasks to the thread pool
  moodycamel::BlockingConcurrentQueue<Simulation::Task> task_queue(queue_capacity, 1, 0);
  moodycamel::ProducerToken ptok(task_queue);

  try {
    this->_tpool.push_task([&]() {  // This thread is in charge of writing contacts to disk
      this->write_contacts_to_disk(progress_queue, progress_queue_mutex);
    });

    for (uint64_t tid = 0; tid < this->nthreads; ++tid) {  // Start simulation threads
      this->_tpool.push_task([&, tid]() {
        this->simulate_worker(tid, task_queue, progress_queue, progress_queue_mutex,
                              task_batch_size_deq);
      });
    }

    // The remaining code submits simulation tasks to the queue. Then it waits until all the tasks
    // have been completed, and contacts have been written to disk

    absl::FixedArray<Task> tasks(task_batch_size_enq);
    auto taskid = 0UL;

    // Loop over chromosomes
    for (auto& chrom : this->_genome) {
      if (!this->ok()) {
        this->handle_exceptions();
      }
      // Don't simulate KO chroms (but write them to disk if the user desires so)
      if (!chrom.ok() || chrom.num_barriers() == 0) {
        spdlog::info(FMT_STRING("SKIPPING '{}'..."), chrom.name());
        if (this->write_contacts_for_ko_chroms) {
          std::scoped_lock l(progress_queue_mutex);
          progress_queue.emplace_back(&chrom, num_cells);
        }
        continue;
      }

      {
        std::scoped_lock l(progress_queue_mutex);

        // Signal that we have started processing the current chrom
        progress_queue.emplace_back(&chrom, 0UL);
      }

      // Compute # of LEFs to be simulated based on chrom. sizes
      const auto nlefs = static_cast<size_t>(
          std::max(1.0, std::round(this->number_of_lefs_per_mbp *
                                   (static_cast<double>(chrom.simulated_size()) / Mbp))));

      auto target_contacts = 0UL;
      if (this->target_contact_density != 0) {  // Compute the number of simulation rounds required
        // to reach the target contact density
        target_contacts = static_cast<size_t>(std::max(
            1.0,
            std::round((this->target_contact_density *
                        static_cast<double>(chrom.npixels(this->diagonal_width, this->bin_size))) /
                       static_cast<double>(this->num_cells))));
      }
      auto target_epochs = target_contacts == 0UL ? this->simulation_iterations
                                                  : (std::numeric_limits<size_t>::max)();

      size_t cellid = 0;
      const auto nbatches = (this->num_cells + task_batch_size_enq - 1) / task_batch_size_enq;
      for (auto batchid = 0UL; batchid < nbatches; ++batchid) {
        // Generate a batch of tasks for all the simulations involving the current chrom
        std::generate(tasks.begin(), tasks.end(), [&]() {
          return Task{{taskid++, &chrom, cellid++, target_epochs, target_contacts, nlefs,
                       chrom.barriers().data()}};
        });
        const auto ntasks =
            cellid > this->num_cells ? tasks.size() - (cellid - this->num_cells) : tasks.size();
        auto sleep_us = 100;  // NOLINT(readability-magic-numbers)
        while (               // Enqueue tasks
            !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), ntasks)) {
          if (!this->ok()) {
            this->handle_exceptions();
          }
          sleep_us = std::min(100000, sleep_us * 2);  // NOLINT(readability-magic-numbers)
          std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
        }
      }
    }

    {  // Signal end of simulation to the thread that is writing contacts to disk
      std::scoped_lock l(progress_queue_mutex);
      progress_queue.emplace_back(nullptr, 0UL);
    }
    this->_tpool.wait_for_tasks();
    assert(this->_end_of_simulation);  // NOLINT
    assert(!this->_exception_thrown);  // NOLINT
  } catch (...) {
    this->_exception_thrown = true;
    this->_tpool.paused = true;
    this->_tpool.wait_for_tasks();
    throw;
  }
}

void Simulation::simulate_worker(const uint64_t tid,
                                 moodycamel::BlockingConcurrentQueue<Simulation::Task>& task_queue,
                                 std::deque<std::pair<Chromosome*, size_t>>& progress_queue,
                                 std::mutex& progress_queue_mutex, const size_t task_batch_size) {
  spdlog::info(FMT_STRING("Spawning simulation thread {}..."), tid);

  moodycamel::ConsumerToken ctok(task_queue);
  absl::FixedArray<Task> task_buff(task_batch_size);  // Tasks are dequeued in batch.
  // This is to reduce contention
  // when accessing the queue

  // This state object owns all the buffers and PRNG + seed required in order to simulate loop
  // extrusion for a single cell. Buffers are allocated once and resized, cleared and reused
  // throughout the simulation
  Simulation::State local_state;

  try {
    while (this->ok()) {  // Try to dequeue a batch of tasks
      const auto avail_tasks = task_queue.wait_dequeue_bulk_timed(
          ctok, task_buff.begin(), task_buff.size(), std::chrono::milliseconds(10));
      // Check whether dequeue operation timed-out before any task became available
      if (avail_tasks == 0) {
        if (this->_end_of_simulation) {
          // Reached end of simulation (i.e. all tasks have been processed)
          return;
        }
        // Keep waiting until one or more tasks become available
        continue;
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
          const auto target_epochs = static_cast<size_t>(std::max(
              1.0, std::round(
                       (this->target_contact_density * static_cast<double>(task.chrom->npixels(
                                                           this->diagonal_width, this->bin_size))) /
                       (this->lef_fraction_contact_sampling *
                        static_cast<double>(this->num_cells * task.num_lefs)))));

          spdlog::info(FMT_STRING("Simulating ~{} epochs for '{}' across {} cells..."),
                       target_epochs, task.chrom->name(), num_cells);
        }

        // Start the simulation kernel
        Simulation::simulate_one_cell(local_state);

        // Update progress for the current chrom
        std::scoped_lock l1(progress_queue_mutex);
        auto progress = std::find_if(progress_queue.begin(), progress_queue.end(),
                                     [&](const auto& p) { return task.chrom == p.first; });
        assert(progress != progress_queue.end());  // NOLINT

        if (++progress->second == num_cells && !this->_exception_thrown) {
          // We are done simulating loop-extrusion on task.chrom: print a status update
          spdlog::info(FMT_STRING("Simulation for '{}' successfully completed."),
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
