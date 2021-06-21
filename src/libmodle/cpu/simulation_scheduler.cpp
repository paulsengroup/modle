#include <absl/container/fixed_array.h>          // for FixedArray
#include <absl/container/flat_hash_map.h>        // for flat_hash_map, BitMask, raw_hash_set<...
#include <absl/strings/str_join.h>               // for StrJoin
#include <absl/types/span.h>                     // for Span, MakeConstSpan, MakeSpan
#include <fmt/format.h>                          // for print, FMT_STRING
#include <fmt/ostream.h>                         // for formatbuf<>::int_type, print
#include <moodycamel/blockingconcurrentqueue.h>  // for BlockingConcurrentQueue
#include <moodycamel/concurrentqueue.h>          // for ConsumerToken, ProducerToken

#include <algorithm>  // for max, copy, find_if, generate, sort
#include <atomic>     // for atomic
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>               // for thread_pool
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/exception/get_error_info.hpp>       // for get_error_info
#include <boost/stacktrace/stacktrace.hpp>          // for operator<<
#include <cassert>                                  // for assert
#include <chrono>                                   // for milliseconds
#include <cmath>                                    // for round
#include <cstdio>                                   // for stderr
#include <deque>                                    // for deque, operator-, operator!=, _Deque_...
#include <filesystem>                               // for create_directories, operator<<, remov...
#include <iostream>                                 // for size_t, streamsize, ofstream, operator<<
#include <iterator>                                 // for move_iterator, make_move_iterator
#include <limits>                                   // for numeric_limits
#include <memory>                                   // for unique_ptr, make_unique
#include <mutex>                                    // for mutex
#include <sstream>                                  // for basic_stringbuf<>::int_type, basic_st...
#include <string>                                   // for basic_string
#include <thread>                                   // for operator<<, get_id, sleep_for
#include <utility>                                  // for pair, addressof
#include <vector>                                   // for vector

#include "modle/common/config.hpp"       // for Config
#include "modle/common/utils.hpp"        // for traced
#include "modle/contacts.hpp"            // for ContactMatrix
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier
#include "modle/genome.hpp"              // for Chromosome
#include "modle/simulation.hpp"

namespace modle {

void Simulation::run_base() {
  if (!this->skip_output) {  // Write simulation params to file
    if (this->force) {
      std::filesystem::remove_all(this->path_to_output_file_cool);
    }
    std::filesystem::create_directories(this->path_to_output_file_cool.parent_path());
    std::ofstream log_file(this->path_to_log_file);
    if (log_file) {
      fmt::print(log_file, FMT_STRING("{}\n{}\n"), Config::to_string(),
                 absl::StrJoin(this->argv, this->argv + this->argc, " "));
    } else {
      fmt::print(
          stderr,
          FMT_STRING("WARNING: Unable to open log file {} for writing. Continuing anyway..."),
          this->path_to_log_file);
    }
  }

  auto tpool = Simulation::instantiate_thread_pool(this->nthreads + 1, false);

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

  std::atomic<bool> end_of_simulation = false;
  std::mutex progress_queue_mutex;  // Protect rw access to progress_queue
  std::deque<std::pair<Chromosome*, size_t>> progress_queue;

  // Barriers are shared across threads that are simulating loop-extrusion on a given chromosome.
  // Extrusion barriers are placed on this std::deque, and are passed to simulation threads
  // through a Span of const ExtrusionBarriers
  std::mutex barrier_mutex;
  absl::flat_hash_map<Chromosome*, std::unique_ptr<std::vector<ExtrusionBarrier>>> barriers;

  constexpr size_t task_batch_size = 128;
  // Queue used to submit simulation tasks to the thread pool
  moodycamel::BlockingConcurrentQueue<Simulation::Task> task_queue(
      std::min(static_cast<size_t>(static_cast<double>(this->ncells) * 1.1),
               this->nthreads * task_batch_size),
      1, 0);
  moodycamel::ProducerToken ptok(task_queue);

  boost::asio::post(tpool, [&]() {  // This thread is in charge of writing contacts to disk
    this->write_contacts_to_disk(progress_queue, progress_queue_mutex, end_of_simulation);
  });

  for (auto i = 0UL; i < this->nthreads; ++i) {  // Start simulation threads
    boost::asio::post(tpool, [&]() {
      this->worker(task_queue, progress_queue, progress_queue_mutex, barrier_mutex, barriers,
                   end_of_simulation);
    });
  }

  // The remaining code submits simulation tasks to the queue. Then it waits until all the tasks
  // have been completed, and contacts have been written to disk

  absl::Span<const ExtrusionBarrier> extr_barriers_buff{};
  absl::FixedArray<Task> tasks(task_batch_size);
  auto taskid = 0UL;

  // Loop over chromosomes
  for (auto& chrom : this->_genome) {
    // Don't simulate KO chroms (but write them to disk if the user desires so)
    if (!chrom.ok()) {
      fmt::print(stderr, "SKIPPING '{}'...\n", chrom.name());
      if (this->write_contacts_for_ko_chroms) {
        std::scoped_lock l(barrier_mutex, progress_queue_mutex);
        progress_queue.emplace_back(&chrom, ncells);
        barriers.emplace(&chrom, nullptr);
      }
      continue;
    }

    // Allocate the contact matrix. Once it's not needed anymore, the contact matrix will be
    // de-allocated by the thread that is writing contacts to disk
    chrom.allocate_contacts(this->bin_size, this->diagonal_width);
    {
      // For consistency, it is important that both locks are held while new items are added to the
      // two queues
      std::scoped_lock l(barrier_mutex, progress_queue_mutex);

      // Signal that we have started processing the current chrom
      progress_queue.emplace_back(&chrom, 0UL);

      // Allocate extr. barriers mapping on the chromosome that is being simulated
      auto node = barriers.emplace(
          &chrom,
          std::make_unique<std::vector<ExtrusionBarrier>>(this->_genome.generate_vect_of_barriers(
              chrom.name(), this->_config->ctcf_occupied_self_prob,
              this->_config->ctcf_not_occupied_self_prob)));
      extr_barriers_buff = absl::MakeConstSpan(*node.first->second);
    }

    // Compute # of LEFs to be simulated based on chrom. sizes
    const auto nlefs = static_cast<size_t>(std::round(
        this->number_of_lefs_per_mbp * (static_cast<double>(chrom.simulated_size()) / 1.0e6)));

    auto target_contacts = 0UL;
    if (this->target_contact_density != 0) {  // Compute the number of simulation rounds required
      // to reach the target contact density
      target_contacts = static_cast<size_t>(
          std::max(1.0, std::round((this->target_contact_density *
                                    static_cast<double>(chrom.contacts().npixels())) /
                                   static_cast<double>(this->ncells))));
    }
    auto target_epochs =
        target_contacts == 0UL ? this->simulation_iterations : std::numeric_limits<size_t>::max();

    size_t cellid = 0;
    const auto nbatches = (this->ncells + task_batch_size - 1) / task_batch_size;
    for (auto batchid = 0UL; batchid < nbatches; ++batchid) {
      // Generate a batch of tasks for all the simulations involving the current chrom
      std::generate(tasks.begin(), tasks.end(), [&]() {
        return Task{taskid++,        &chrom, cellid++,          target_epochs,
                    target_contacts, nlefs,  extr_barriers_buff};
      });
      const auto ntasks =
          cellid > this->ncells ? tasks.size() - (cellid - this->ncells) : tasks.size();
      auto sleep_us = 100;  // NOLINT
      while (               // Enqueue tasks
          !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), ntasks)) {
        sleep_us = std::min(100000, sleep_us * 2);  // NOLINT
        std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
      }
    }
  }

  {  // Signal end of simulation to the thread that is writing contacts to disk
    std::scoped_lock l(progress_queue_mutex);
    progress_queue.emplace_back(nullptr, 0UL);
  }
  tpool.join();
  assert(end_of_simulation);  // NOLINT
}

void Simulation::worker(
    moodycamel::BlockingConcurrentQueue<Simulation::Task>& task_queue,
    std::deque<std::pair<Chromosome*, size_t>>& progress_queue, std::mutex& progress_queue_mutex,
    std::mutex& barrier_mutex,
    absl::flat_hash_map<Chromosome*, std::unique_ptr<std::vector<ExtrusionBarrier>>>& barriers,
    std::atomic<bool>& end_of_simulation, size_t task_batch_size) const {
  fmt::print(stderr, FMT_STRING("Spawning simulation thread {}...\n"), std::this_thread::get_id());

  moodycamel::ConsumerToken ctok(task_queue);
  absl::FixedArray<Task> task_buff(task_batch_size);  // Tasks are dequeue in batch.
  // This is to reduce contention when
  // accessing the queue

  // This state object owns all the buffers and PRNG + seed required in order to simulate loop
  // extrusion for a single cell. Buffers are allocated once and resized, cleared and reused
  // throughout the simulation
  Simulation::State state_buff;

  try {
    while (true) {  // Try to dequeue a batch of tasks
      const auto avail_tasks = task_queue.wait_dequeue_bulk_timed(
          ctok, task_buff.begin(), task_buff.size(), std::chrono::milliseconds(10));
      // Check whether dequeue operation timed-out before any task became available
      if (avail_tasks == 0) {
        if (end_of_simulation) {
          // Reached end of simulation (i.e. all tasks have been processed)
          return;
        }
        // Keep waiting until one or more tasks become available
        continue;
      }

      // Loop over new tasks
      auto tasks = absl::MakeSpan(task_buff.data(), avail_tasks);
      for (auto& task : tasks) {
        // Resize and reset buffers
        state_buff = task;    // Set simulation state based on task data
        state_buff.resize();  // This resizes buffers based on the nlefs to be simulated
        state_buff.reset();   // Clear all buffers

        if (task.cell_id == 0) {
          // Print a status update when we are processing cell #0 for a given chromosome
          const auto target_epochs = static_cast<size_t>(
              std::max(1.0, std::round((this->target_contact_density *
                                        static_cast<double>(task.chrom->contacts().npixels())) /
                                       (this->lef_fraction_contact_sampling *
                                        static_cast<double>(this->ncells * task.nlefs)))));

          fmt::print(stderr, FMT_STRING("Simulating ~{} epochs for '{}' across {} cells...\n"),
                     target_epochs, task.chrom->name(), ncells);
        }

        // Start the simulation kernel
        Simulation::simulate_extrusion_kernel(state_buff);

        // Update progress for the current chrom
        std::scoped_lock l1(progress_queue_mutex);
        auto progress = std::find_if(progress_queue.begin(), progress_queue.end(),
                                     [&](const auto& p) { return task.chrom == p.first; });
        assert(progress != progress_queue.end());  // NOLINT

        if (++progress->second == ncells) {
          // We are done simulating loop-extrusion on task.chrom
          // Print a status update and deallocate extr. barriers
          fmt::print(stderr, "Simulation for '{}' successfully completed.\n", task.chrom->name());
          std::scoped_lock l2(barrier_mutex);
          barriers.erase(task.chrom);
        }
      }
    }
  } catch (const std::exception& err) {
    // This is needed, as exceptions don't seem to always propagate to the main()
    // TODO: Find a better way to propagate exception up to the main thread. Maybe use a
    // concurrent queue?
    fmt::print(stderr, FMT_STRING("Detected an error in thread {}:\n{}\n{}\n"),
               std::this_thread::get_id(), state_buff.to_string(), err.what());
#ifndef BOOST_STACKTRACE_USE_NOOP
    const auto* st = boost::get_error_info<modle::utils::traced>(err);
    if (st) {
      std::cerr << *st << '\n';
    } else {
      fmt::print(stderr, "Stack trace not available!\n");
    }
#endif
    throw;
  }
}

void Simulation::run_pairwise() {
  assert(std::filesystem::exists(this->path_to_output_file_cool));  // NOLINT
}

}  // namespace modle
