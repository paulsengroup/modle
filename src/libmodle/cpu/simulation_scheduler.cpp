#include <absl/container/fixed_array.h>          // for FixedArray
#include <absl/container/flat_hash_map.h>        // for flat_hash_map, BitMask, raw_hash_set<...
#include <absl/strings/str_join.h>               // for StrJoin
#include <absl/types/span.h>                     // for Span, MakeConstSpan, MakeSpan
#include <fmt/compile.h>                         // for FMT_COMPILE
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
#include <boost/filesystem/path.hpp>                // for create_directories, operator<<, remov...
#include <boost/stacktrace/stacktrace.hpp>          // for operator<<
#include <cassert>                                  // for assert
#include <chrono>                                   // for milliseconds
#include <cmath>                                    // for round
#include <cstdio>                                   // for stderr
#include <deque>                                    // for deque, operator-, operator!=, _Deque_...
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

#include "modle/common/config.hpp"  // for Config
#include "modle/common/utils.hpp"   // for traced
#include "modle/compressed_io.hpp"  // for Writer, write
#include "modle/contacts.hpp"       // for ContactMatrix
#include "modle/cooler.hpp"
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier
#include "modle/genome.hpp"              // for Chromosome
#include "modle/simulation.hpp"

namespace modle {

void Simulation::run_simulation() {
  if (!this->skip_output) {  // Write simulation params to file
    assert(boost::filesystem::exists(this->path_to_output_prefix.parent_path()));
    if (this->force) {
      boost::filesystem::remove(this->path_to_output_file_cool);
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

  boost::asio::post(tpool, [&]() {  // This thread is in charge of writing contacts to disk
    this->write_contacts_to_disk(progress_queue, progress_queue_mutex, end_of_simulation);
  });

  for (auto i = 0UL; i < this->nthreads; ++i) {  // Start simulation threads
    boost::asio::post(tpool, [&]() {
      this->worker(task_queue, progress_queue, progress_queue_mutex, end_of_simulation,
                   task_batch_size_deq);
    });
  }

  // The remaining code submits simulation tasks to the queue. Then it waits until all the tasks
  // have been completed, and contacts have been written to disk

  absl::FixedArray<Task> tasks(task_batch_size_enq);
  auto taskid = 0UL;

  // Loop over chromosomes
  for (auto& chrom : this->_genome) {
    // Don't simulate KO chroms (but write them to disk if the user desires so)
    if (!chrom.ok() || chrom.num_barriers() == 0) {
      fmt::print(stderr, "SKIPPING '{}'...\n", chrom.name());
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
    const auto nlefs = static_cast<size_t>(std::round(
        this->number_of_lefs_per_mbp * (static_cast<double>(chrom.simulated_size()) / Mbp)));

    auto target_contacts = 0UL;
    if (this->target_contact_density != 0) {  // Compute the number of simulation rounds required
      // to reach the target contact density
      target_contacts = static_cast<size_t>(std::max(
          1.0,
          std::round((this->target_contact_density *
                      static_cast<double>(chrom.npixels(this->diagonal_width, this->bin_size))) /
                     static_cast<double>(this->num_cells))));
    }
    auto target_epochs =
        target_contacts == 0UL ? this->simulation_iterations : std::numeric_limits<size_t>::max();

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

void Simulation::worker(moodycamel::BlockingConcurrentQueue<Simulation::Task>& task_queue,
                        std::deque<std::pair<Chromosome*, size_t>>& progress_queue,
                        std::mutex& progress_queue_mutex, std::atomic<bool>& end_of_simulation,
                        size_t task_batch_size) const {
  fmt::print(stderr, FMT_STRING("Spawning simulation thread {}...\n"), std::this_thread::get_id());

  moodycamel::ConsumerToken ctok(task_queue);
  absl::FixedArray<Task> task_buff(task_batch_size);  // Tasks are dequeued in batch.
                                                      // This is to reduce contention
                                                      // when accessing the queue

  // This state object owns all the buffers and PRNG + seed required in order to simulate loop
  // extrusion for a single cell. Buffers are allocated once and resized, cleared and reused
  // throughout the simulation
  Simulation::State local_state;

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

          fmt::print(stderr, FMT_STRING("Simulating ~{} epochs for '{}' across {} cells...\n"),
                     target_epochs, task.chrom->name(), num_cells);
        }

        // Allocate the contact matrix. Once it's not needed anymore, the contact matrix will be
        // de-allocated by the thread that is writing contacts to disk
        local_state.chrom->allocate_contacts(this->bin_size, this->diagonal_width);
        // Start the simulation kernel
        Simulation::simulate_one_cell(local_state);

        // Update progress for the current chrom
        std::scoped_lock l1(progress_queue_mutex);
        auto progress = std::find_if(progress_queue.begin(), progress_queue.end(),
                                     [&](const auto& p) { return task.chrom == p.first; });
        assert(progress != progress_queue.end());  // NOLINT

        if (++progress->second == num_cells) {
          // We are done simulating loop-extrusion on task.chrom: print a status update
          fmt::print(stderr, "Simulation for '{}' successfully completed.\n", task.chrom->name());
        }
      }
    }
  } catch (const std::exception& err) {
    // This is needed, as exceptions don't seem to always propagate to the main()
    // TODO: Find a better way to propagate exception up to the main thread. Maybe use a
    // concurrent queue?
    fmt::print(stderr, FMT_STRING("Detected an error in thread {}:\n{}\n{}\n"),
               std::this_thread::get_id(), local_state.to_string(), err.what());
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

bed::BED_tree<> Simulation::import_deletions() const {
  fmt::print(stderr, FMT_STRING("Importing deletions from file {}...\n"),
             this->path_to_deletion_bed);
  const auto t0 = absl::Now();
  auto deletions = bed::BED_tree<>{this->path_to_deletion_bed, bed::BED::BED3};
  fmt::print(stderr, FMT_STRING("Imported {} deletions in {}\n"), deletions.size(),
             absl::FormatDuration(absl::Now() - t0));
  return deletions;
}

bed::BED_tree<> Simulation::generate_deletions() const {
  fmt::print(stderr, FMT_STRING("Generating deletions for {} chromosomes...\n"),
             this->_genome.size());
  const auto t0 = absl::Now();
  bed::BED_tree<> deletions;
  for (const auto& chrom : this->_genome) {
    for (const auto& barrier : chrom.barriers().data()) {
      // Compute the span of the deletion
      const auto deletion_size_ = std::min(barrier.pos() + 1, this->deletion_size);
      const auto deletion_begin = barrier.pos() + 1 - deletion_size_;
      const auto deletion_end = deletion_begin + deletion_size_;
      deletions.insert(std::string{chrom.name()}, deletion_begin, deletion_end);
    }
  }
  deletions.index();
  fmt::print(stderr, FMT_STRING("Generated {} deletions in {}\n"), deletions.size(),
             absl::FormatDuration(absl::Now() - t0));
  return deletions;
}

void Simulation::run_perturbate() {
  if (!this->skip_output) {  // Write simulation params to file
    assert(boost::filesystem::exists(this->path_to_output_prefix.parent_path()));
    if (this->force) {
      boost::filesystem::remove(this->path_to_output_file_bedpe);
    }
  }
  // TODO Do proper error handling
  // For now we support only the case where exactly two files are specified
  assert(this->path_to_feature_bed_files.size() == 2);  // NOLINT
  const auto regions_for_contact_output =
      bed::Parser(this->path_to_regions_for_contact_output_bed, bed::BED::BED3)
          .parse_all_in_interval_tree();

  constexpr size_t task_batch_size_enq = 32;  // NOLINTNEXTLINE
  moodycamel::BlockingConcurrentQueue<TaskPW> task_queue(this->nthreads * 2, 1, 0);
  moodycamel::ProducerToken ptok(task_queue);
  std::array<TaskPW, task_batch_size_enq> tasks;

  std::mutex out_stream_mutex;
  std::mutex cooler_mutex;
  auto out_stream = compressed_io::Writer(this->path_to_output_file_bedpe);

  std::atomic<bool> end_of_simulation = false;

  auto tpool = Simulation::instantiate_thread_pool();
  for (auto i = 0UL; i < this->nthreads; ++i) {  // Start simulation threads
    boost::asio::post(tpool, [&]() {
      this->worker(task_queue, out_stream, out_stream_mutex, cooler_mutex, end_of_simulation);
    });
  }

  const auto all_deletions =
      this->path_to_deletion_bed.empty() ? this->generate_deletions() : this->import_deletions();

  size_t task_id = 0;
  size_t num_tasks = 0;
  for (auto& chrom : this->_genome) {
    const auto chrom_name = std::string{chrom.name()};
    const auto& features = chrom.get_features();
    const auto& barriers = chrom.barriers();
    if (features.size() != 2 || barriers.empty()) {
      continue;  // Skip chromosomes that have less than two "kinds" of features or no barriers
    }

    if (!all_deletions.contains(chrom_name)) {
      continue;
    }

    // The idea here is that given a diagonal width D and a feature F1, in order to track all
    // possible interactions between F1 and nearby features we need to simulate at least window
    // 2*D wide centered around F1. If we now consider another feature F2 located 100 bp
    // downstream of F1, in order to track all possible interactions we need to simulate a 2*D
    // window centered around F2. Given how close F1 and F2 we basically end up simulating the
    // same region twice. This can be avoided by extending the window we are simulating by 1*D
    // left and right, increasing the window width to 4*D. We then define an active window, which
    // is a 2*D wide region centered inside the 4*D window. Now we simulate loop extrusion over a
    // 4*D region and record all possible, non-zero contacts between features falling in the 2*D
    // region. Next, we advance both windows by 1*D, and repeat the same procedure. With this
    // logic the overlap between subsequent windows is always 3*D. Deletions are performed on the
    // entire 4*D region, this is because deletions outside of the active window can still affect
    // interactions between features in the active window. The last problem we have to solve, is
    // to figure out how to ensure that we don't output the number of contacts for a pair of
    // features and a given barrier configuration twice. This is likely to happen, as the active
    // window is 2*D and windows are advanced by only 1*D at a time. The solution to this problem
    // was to use the center of the two windows as partition point. For features upstream of the
    // partition point we output all possible pairwise contacts, while for features located
    // downstream of the partition point we output contacts only if one of the feature
    // participating in the pair is upstream of the partition point. Consider the following
    // example: Given four features F1, F2, F3 and F4, where F1 and F2 are located upstream of the
    // partition point P, while F3 and F4 are located downstream of said point. All features fall
    // in a D wide region in the active window. For each distinct extrusion barrier configuration
    // we output the number of contacts for the following pairs (assumin nnz contacts):
    //
    // - F1:F2 - both are upstream of P
    // - F1:F3 - F1 is upstream of P
    // - F1:F4 - F1 is upstream of P
    // - F2:F3 - F2 is upstream of P
    // - F2:F4 - F2 is upstream of P
    //
    // - Contacts between F3 and F4 will be produced when simulating the next window, as in that
    // case they will both be located upstream of the new partition point
    //
    // It should be noted that the first and last window on a chromosome need special handling, as
    // the outer window cannot extend past the chromosomal boundaries. In these cases the active
    // window is extended such that it goes from 5'-end->3*D for the first window and from P +
    // 1*D->3'-end, where P is the partition point of the previous window.

    size_t cell_id = 0;

    // Setup the
    TaskPW base_task{};

    base_task.chrom = &chrom;
    // Initialize task with the initial window
    base_task.window_start = 0;
    base_task.window_end = 4 * this->diagonal_width;
    base_task.active_window_start = 0;
    base_task.active_window_end = 3 * this->diagonal_width;

    do {
      // Find all barriers falling within the outer window. Skip over windows with 0 barriers
      if (!Simulation::map_barriers_to_window(base_task, chrom)) {
        continue;
      }

      // Find all features of type 1 falling within the outer window. Skip over windows with 0
      // features
      if (!Simulation::map_features_to_window(base_task, chrom)) {
        continue;
      }

      base_task.write_contacts_to_disk = regions_for_contact_output.contains_overlap(
          chrom_name, base_task.window_start, base_task.window_end);

      assert(!base_task.barriers.empty());  // NOLINT
      assert(!base_task.feats1.empty());    // NOLINT
      assert(!base_task.feats2.empty());    // NOLINT

      const auto deletions = [&]() {  // Find deletions mapping to the outer window
        const auto [first_deletion, last_deletion] =
            all_deletions.at(chrom_name).equal_range(base_task.window_start, base_task.window_end);
        return absl::MakeConstSpan(&(*first_deletion), &(*last_deletion));
      }();

      for (const auto& deletion : deletions) {
        // Add task to the current batch
        auto& t = (tasks[num_tasks++] = base_task);  // NOLINT

        // Complete task setup
        t.id = task_id++;
        t.cell_id = cell_id++;
        t.deletion_begin = deletion.chrom_start;
        t.deletion_size = deletion.chrom_end - deletion.chrom_start;
        t.window_end = std::min(t.window_end, chrom.end_pos());
        t.active_window_end = std::min(t.active_window_end, chrom.end_pos());

        // Compute the number of simulation rounds required to reach the target contact density
        if (this->target_contact_density != 0) {
          // Compute the number of pixels mapping to the outer window
          const auto npix1 = (t.window_end - t.window_start + this->bin_size - 1) / this->bin_size;
          const auto npix2 = (this->diagonal_width + this->bin_size - 1) / this->bin_size;

          t.num_target_contacts = static_cast<size_t>(std::max(
              1.0, std::round(this->target_contact_density * static_cast<double>(npix1 * npix2))));
        }

        // Compute the number of LEFs based on the window size
        t.num_lefs = static_cast<size_t>(
            std::round((static_cast<double>(t.window_end - t.window_start) / Mbp) *
                       this->number_of_lefs_per_mbp));

        // Compute the target number of epochs based on the target number of contacts
        t.num_target_epochs = t.num_target_contacts == 0UL ? this->simulation_iterations
                                                           : std::numeric_limits<size_t>::max();

        if (num_tasks == tasks.size()) {  // Enqueue a batch of tasks
          auto sleep_us = 100;            // NOLINT
          while (!task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()),
                                              num_tasks)) {
            sleep_us = std::min(100000, sleep_us * 2);  // NOLINT
            std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
          }
          num_tasks = 0;
        }
      }
    } while (Simulation::advance_window(base_task, chrom));
  }

  // Submit any remaining task
  if (num_tasks != 0) {
    while (!task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), num_tasks)) {
      std::this_thread::sleep_for(std::chrono::microseconds(100));  // NOLINT
    }
  }

  end_of_simulation = true;
  tpool.join();  // Wait on worker threads
}

void Simulation::worker(moodycamel::BlockingConcurrentQueue<Simulation::TaskPW>& task_queue,
                        compressed_io::Writer& out_stream, std::mutex& out_stream_mutex,
                        std::mutex& cooler_mutex, std::atomic<bool>& end_of_simulation,
                        size_t task_batch_size) const {
  fmt::print(stderr, FMT_STRING("Spawning simulation thread {}...\n"), std::this_thread::get_id());
  moodycamel::ConsumerToken ctok(task_queue);

  absl::FixedArray<TaskPW> task_buff(task_batch_size);  // Tasks are dequeue in batch.

  Simulation::StatePW state{};

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
      for (const auto& task : absl::MakeSpan(task_buff.data(), avail_tasks)) {
        assert(!task.barriers.empty());  // NOLINT
        assert(!task.feats1.empty());    // NOLINT
        assert(!task.feats2.empty());    // NOLINT

        state = task;  // Set simulation state based on task data
        state.contacts.resize(state.window_end - state.window_start, this->diagonal_width,
                              this->bin_size);

        Simulation::simulate_window(state, out_stream, out_stream_mutex, cooler_mutex);
      }
    }
  } catch (const std::exception& err) {
    // This is needed, as exceptions don't seem to always propagate to main()
    // TODO: Find a better way to propagate exception up to the main thread. Maybe use a
    // concurrent queue?
    fmt::print(stderr, FMT_STRING("Detected an error in thread {}:\n{}\n{}\n"),
               std::this_thread::get_id(), state.to_string(), err.what());
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

void Simulation::simulate_window(Simulation::StatePW& state, compressed_io::Writer& out_stream,
                                 std::mutex& out_stream_mutex, std::mutex& cooler_mutex) const {
  fmt::print(stderr, FMT_STRING("Processing {}[{}-{}]; outer_window=[{}-{}]; deletion=[{}-{}];\n"),
             state.chrom->name(), state.active_window_start, state.active_window_end,
             state.window_start, state.window_end, state.deletion_begin,
             state.deletion_begin + state.deletion_size);
  // Simulating a window consists in generating all valid combinations of barriers after deleting
  // a portion of the DNA from the window. The idea here is that given that biologically relevant
  // boundaries are often marked by multiple extrusion barriers, deleting one barrier at a time is
  // not very interesting, while by deleting regions spanning several kbs we are able to delete
  // clusters of barriers, which should produce a more dramatic change in the number of contacts.
  // Deletions are performed according to the following algorithm:
  //
  // Going from the first barrier at the 5-end of the window and moving towards the last barrier
  // at the 3-'end of the window do:
  // - Delete a region going from the BP + 1 - DS to BP, where BP is the position of the extrusion
  //   barrier and DS is the deletion size. This will produce a novel extrusion barrier
  //   configuration
  // - Simulate loop extrusion on the current window given the above barrier configuration
  // - Move to the next barrier and start over from step #1
  //
  // Thew idea here is that if there's a cluster of barriers going from BP1 until BP10, and
  // BP10 - BP1 < DS, when BP points to BP10, the deletion will include every barrier in the
  // cluster.

  const auto write_to_stdout = out_stream.path().empty();

  std::string out_buffer;
  std::string barrier_str_buff;
  // Figure out whether we are processing the first or last window and compute the partition
  // point
  // const auto first_window = state.active_window_start == state.chrom->start_pos();
  const auto last_window = state.active_window_end == state.chrom->end_pos();
  const auto partition_point = state.active_window_start + this->diagonal_width;

  // Generate barrier configuration
  size_t num_active_barriers = 0;
  state.barrier_tmp_buff.resize(state.barriers.size());
  std::transform(state.barriers.begin(), state.barriers.end(), state.barrier_tmp_buff.begin(),
                 [&](const auto& barrier) {
                   if (barrier.pos() >= state.deletion_begin &&
                       barrier.pos() < state.deletion_begin + state.deletion_size) {
                     // Do not change ( -> {. We don't want to construct the object using the
                     // initializer list
                     return ExtrusionBarrier(barrier.pos(), 0.0, 1.0,
                                             barrier.blocking_direction_major());
                   }
                   ++num_active_barriers;
                   return barrier;
                 });

  state.barriers = absl::MakeConstSpan(state.barrier_tmp_buff);

  // Resize and reset state buffers
  state.resize_buffers();
  state.reset_buffers();
  state.contacts.resize(state.window_end - state.window_start, this->diagonal_width,
                        this->bin_size);
  state.contacts.reset();

  Simulation::simulate_one_cell(state);

  // Output contacts for valid pairs of features
  for (auto i = 0UL; i < state.feats1.size(); ++i) {
    const auto& feat1 = state.feats1[i];
    const auto feat1_abs_center_pos = (feat1.chrom_start + feat1.chrom_end + 1) / 2;
    for (auto j = i; j < state.feats2.size(); ++j) {
      const auto& feat2 = state.feats2[j];
      const auto feat2_abs_center_pos = (feat2.chrom_start + feat2.chrom_end + 1) / 2;

      // Don't output contacts when both feat1 and feat2 are located downstream of the partition
      // point. This does not apply when we are processing the last window
      if (!last_window && feat1.chrom_start >= partition_point &&
          feat2.chrom_start >= partition_point) {
        continue;
      }

      // Convert absolute positions to relative positions, as the contact matrix does not refer
      // to an entire chromosome, but only to a 4*D wide region
      const auto feat1_rel_center_pos = feat1_abs_center_pos - state.window_start;
      const auto feat2_rel_center_pos = feat2_abs_center_pos - state.window_start;

      const auto feat1_rel_bin = feat1_rel_center_pos / this->bin_size;
      const auto feat2_rel_bin = feat2_rel_center_pos / this->bin_size;

      const auto contacts = state.contacts.get(feat1_rel_bin, feat2_rel_bin, this->block_size);
      if (contacts == 0) {  // Don't output entries with 0 contacts
        continue;
      }

      // Compute the absolute bin coordinates. This are used to compute the positions shown in
      // the BEDPE file
      const auto feat1_abs_bin = (feat1_abs_center_pos + this->bin_size - 1) / this->bin_size;
      const auto feat2_abs_bin = (feat2_abs_center_pos + this->bin_size - 1) / this->bin_size;

      // Generate the name field. The field will be "none;none" in case both features don't have
      // a name
      const auto name = absl::StrCat(feat1.name.empty() ? "none" : feat1.name, ";",
                                     feat2.name.empty() ? "none" : feat2.name);
      absl::StrAppend(  // Append a BEDPE record to the local buffer
          &out_buffer,
          fmt::format(FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tdeletion={}-{}"
                                  "\tnum_barr={}/{}\tactive_window={}-{}\touter_window={}-{}\n"),
                      feat1.chrom, feat1_abs_bin * bin_size, (feat1_abs_bin + 1) * bin_size,
                      feat2.chrom, feat2_abs_bin * bin_size, (feat2_abs_bin + 1) * bin_size, name,
                      contacts, feat1.strand, feat2.strand, state.deletion_begin,
                      state.deletion_begin + state.deletion_size, num_active_barriers,
                      state.barriers.size(), state.active_window_start, state.active_window_end,
                      state.window_start, state.window_end));
    }
  }

  // Write the buffer to the appropriate stream
  if (!out_buffer.empty()) {
    if (std::scoped_lock l(out_stream_mutex); write_to_stdout) {
      fmt::print(stdout, out_buffer);
    } else {
      out_stream.write(out_buffer);
    }
    out_buffer.clear();
  }

  if (state.write_contacts_to_disk) {
    const auto file_name = fmt::format(
        FMT_STRING("{}_{:06d}_{}_window_{}-{}_deletion_{}-{}.cool"),
        this->path_to_output_prefix.string(), state.id, state.chrom->name(), state.window_start,
        state.window_end, state.deletion_begin, state.deletion_begin + state.deletion_size);

    std::scoped_lock l(cooler_mutex);
    const auto t0 = absl::Now();
    fmt::print(stderr, FMT_STRING("Writing contacts for {} to file \"{}\"..."), state.chrom->name(),
               file_name);
    {
      auto c = cooler::Cooler(file_name, cooler::Cooler::WRITE_ONLY, this->bin_size,
                              this->_genome.chromosome_with_longest_name().name().size());
      for (const auto& chrom : this->_genome) {
        if (&chrom == state.chrom) {
          c.write_or_append_cmatrix_to_file(state.contacts, state.chrom->name(), state.window_start,
                                            state.window_end, state.chrom->size(), true);
          continue;
        }
        assert(chrom.contacts_ptr() == nullptr);  // NOLINT
        c.write_or_append_empty_cmatrix_to_file(chrom.name(), chrom.start_pos(), chrom.end_pos(),
                                                chrom.size(), true);
      }
    }
    fmt::print(stderr, FMT_STRING(" DONE in {}\n"), absl::FormatDuration(absl::Now() - t0));
  }
}

bool Simulation::advance_window(TaskPW& base_task, const Chromosome& chrom) const {
  base_task.active_window_start += this->diagonal_width;
  base_task.active_window_end =
      std::min(base_task.active_window_start + (2 * this->diagonal_width), chrom.end_pos());

  base_task.window_start = base_task.active_window_start - this->diagonal_width;
  base_task.window_end =
      std::min(base_task.active_window_end + this->diagonal_width, chrom.end_pos());

  return base_task.active_window_start < chrom.end_pos();
}

bool Simulation::map_features_to_window(TaskPW& base_task, const Chromosome& chrom) {
  const auto chrom_name = std::string{chrom.name()};
  const auto& features = chrom.get_features();
  auto print_status_update = [](const auto& t) {
    fmt::print(stderr, "Skipping {}[{}-{}]...\n", t.chrom->name(), t.active_window_start,
               t.active_window_end);
  };

  // Find all features of type 1 falling within the outer window. Skip over windows with 0
  // features
  auto [first_feat1, last_feat1] =
      features[0].equal_range(base_task.active_window_start, base_task.active_window_end);
  if (first_feat1 == features[0].data_end()) {
    print_status_update(base_task);
    return false;
  }

  // Find all features of type 2 falling within the outer window. Skip over windows with 0
  // features
  auto [first_feat2, last_feat2] =
      features[1].equal_range(base_task.active_window_start, base_task.active_window_end);
  if (first_feat2 == features[1].data_end()) {
    print_status_update(base_task);
    return false;
  }

  base_task.feats1 = absl::MakeConstSpan(&(*first_feat1), &(*last_feat1));
  base_task.feats2 = absl::MakeConstSpan(&(*first_feat2), &(*last_feat2));
  return true;
}

bool Simulation::map_barriers_to_window(TaskPW& base_task, const Chromosome& chrom) {
  const auto& barriers = chrom.barriers();

  // Find all barriers falling within the outer window. Skip over windows with 0 barriers
  auto [first_barrier, last_barrier] =
      barriers.equal_range(base_task.window_start, base_task.window_end);
  if (first_barrier == barriers.data_end()) {
    fmt::print(stderr, "Skipping {}[{}-{}]...\n", base_task.chrom->name(),
               base_task.active_window_start, base_task.active_window_end);
    return false;
  }

  base_task.barriers = absl::MakeConstSpan(&(*first_barrier), &(*last_barrier));
  return true;
}

}  // namespace modle
