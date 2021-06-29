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

#include "modle/common/config.hpp"  // for Config
#include "modle/common/utils.hpp"   // for traced
#include "modle/compressed_io.hpp"  // for Writer, write
#include "modle/contacts.hpp"       // for ContactMatrix
#include "modle/cooler.hpp"
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

  constexpr size_t task_batch_size_enq = 128;  // NOLINTNEXTLINE
  const auto queue_capacity =
      std::min(static_cast<size_t>(static_cast<double>(this->num_cells) * 1.1),
               this->nthreads * task_batch_size_enq);
  // Queue used to submit simulation tasks to the thread pool
  moodycamel::BlockingConcurrentQueue<Simulation::Task> task_queue(queue_capacity, 1, 0);
  moodycamel::ProducerToken ptok(task_queue);

  boost::asio::post(tpool, [&]() {  // This thread is in charge of writing contacts to disk
    this->write_contacts_to_disk(progress_queue, progress_queue_mutex, end_of_simulation);
  });

  const auto task_batch_size_deq = std::min(32UL, this->num_cells / this->nthreads);
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
    if (!chrom.ok()) {
      fmt::print(stderr, "SKIPPING '{}'...\n", chrom.name());
      if (this->write_contacts_for_ko_chroms) {
        std::scoped_lock l(progress_queue_mutex);
        progress_queue.emplace_back(&chrom, num_cells);
      }
      continue;
    }

    // Allocate the contact matrix. Once it's not needed anymore, the contact matrix will be
    // de-allocated by the thread that is writing contacts to disk
    chrom.allocate_contacts(this->bin_size, this->diagonal_width);
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
      target_contacts = static_cast<size_t>(
          std::max(1.0, std::round((this->target_contact_density *
                                    static_cast<double>(chrom.contacts().npixels())) /
                                   static_cast<double>(this->num_cells))));
    }
    auto target_epochs =
        target_contacts == 0UL ? this->simulation_iterations : std::numeric_limits<size_t>::max();

    size_t cellid = 0;
    const auto nbatches = (this->num_cells + task_batch_size_enq - 1) / task_batch_size_enq;
    for (auto batchid = 0UL; batchid < nbatches; ++batchid) {
      // Generate a batch of tasks for all the simulations involving the current chrom
      std::generate(tasks.begin(), tasks.end(), [&]() {
        return Task{{taskid++, &chrom, cellid++, target_epochs, target_contacts, nlefs}};
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
                                        static_cast<double>(this->num_cells * task.num_lefs)))));

          fmt::print(stderr, FMT_STRING("Simulating ~{} epochs for '{}' across {} cells...\n"),
                     target_epochs, task.chrom->name(), num_cells);
        }

        // Start the simulation kernel
        Simulation::simulate_extrusion_kernel(state_buff);

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
  // TODO Do proper error handling
  // For now we support only the case where exactly two files are specified
  assert(this->path_to_feature_bed_files.size() == 2);  // NOLINT

  constexpr size_t task_batch_size_enq = 32;  // NOLINTNEXTLINE
  moodycamel::BlockingConcurrentQueue<TaskPW> task_queue(this->nthreads * 2, 1, 0);
  moodycamel::ProducerToken ptok(task_queue);
  std::array<StatePW, task_batch_size_enq> tasks;

  std::mutex out_stream_mutex;
  auto out_stream = compressed_io::Writer(this->path_to_output_file_bedpe);

  std::atomic<bool> end_of_simulation = false;

  auto tpool = Simulation::instantiate_thread_pool();
  for (auto i = 0UL; i < this->nthreads; ++i) {  // Start simulation threads
    boost::asio::post(tpool, [&]() {
      this->worker(task_queue, out_stream, out_stream_mutex, end_of_simulation);
    });
  }

  auto print_status_update = [](const auto& t) {
    fmt::print(stderr, "Skipping {}[{}-{}]...\n", t.chrom->name(), t.active_window_start,
               t.active_window_end);
  };

  size_t task_id = 0;
  size_t num_tasks = 0;
  for (auto& chrom : this->_genome) {
    const auto& features = chrom.get_features();
    const auto& barriers = chrom.barriers();
    if (features.size() != 2 || barriers.empty()) {
      continue;
    }

    size_t cell_id = 0;
    TaskPW t{};
    t.chrom = &chrom;
    t.window_start = 0;
    t.window_end = 4 * this->diagonal_width;
    t.active_window_start = 0;
    t.active_window_end = 3 * this->diagonal_width;

    auto advance_window = [&]() -> bool {
      t.window_start += this->diagonal_width;
      t.active_window_start = t.window_start + this->diagonal_width;

      t.active_window_end = std::min(t.active_window_start + this->diagonal_width, chrom.end_pos());
      t.window_end = std::min(t.active_window_end + this->diagonal_width, chrom.end_pos());

      return t.active_window_start >= chrom.end_pos();
    };

    do {
      if (this->target_contact_density != 0) {  // Compute the number of simulation rounds required
                                                // to reach the target contact density

        const auto num_pixels = [&]() {
          const auto npix1 =
              (t.window_end - this->deletion_size - t.window_start + this->bin_size - 1) /
              this->bin_size;
          const auto npix2 = (this->diagonal_width + this->bin_size - 1) / this->bin_size;

          return npix1 * npix2;
        }();
        t.num_target_contacts = static_cast<size_t>(std::max(
            1.0, std::round(this->target_contact_density * static_cast<double>(num_pixels))));
      }
      t.num_target_epochs = t.num_target_contacts == 0UL ? this->simulation_iterations
                                                         : std::numeric_limits<size_t>::max();

      t.num_lefs = static_cast<size_t>(static_cast<double>(t.window_end - t.window_start) *
                                       this->number_of_lefs_per_mbp);

      auto [first_barrier, last_barrier] = barriers.equal_range(t.window_start, t.window_end);
      if (first_barrier == barriers.data_end()) {
        print_status_update(t);
        continue;
      }

      auto [first_feat1, last_feat1] =
          features[0].equal_range(t.active_window_start, t.active_window_end);
      if (first_feat1 == features[0].data_end()) {
        print_status_update(t);
        continue;
      }

      auto [first_feat2, last_feat2] =
          features[1].equal_range(t.active_window_start, t.active_window_end);
      if (first_feat2 == features[1].data_end()) {
        print_status_update(t);
        continue;
      }

      t.id = task_id++;
      t.cell_id = cell_id++;

      t.barriers =
          absl::MakeConstSpan(&(*first_barrier), static_cast<size_t>(last_barrier - first_barrier));
      t.feats1 =
          absl::MakeConstSpan(&(*first_feat1), static_cast<size_t>(last_feat1 - first_feat1));
      t.feats2 =
          absl::MakeConstSpan(&(*first_feat2), static_cast<size_t>(last_feat2 - first_feat2));

      if (num_tasks == tasks.size()) {
        auto sleep_us = 100;  // NOLINT
        while (               // Enqueue tasks
            !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), num_tasks)) {
          sleep_us = std::min(100000, sleep_us * 2);  // NOLINT
          std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
        }
        num_tasks = 0;
      }
      tasks[num_tasks++] = t;  // NOLINT
    } while (!advance_window());
  }
  if (num_tasks != 0) {
    while (  // Enqueue tasks
        !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), num_tasks)) {
      std::this_thread::sleep_for(std::chrono::microseconds(100));  // NOLINT
    }
  }
  end_of_simulation = true;
  tpool.join();
}

void Simulation::worker(moodycamel::BlockingConcurrentQueue<Simulation::TaskPW>& task_queue,
                        compressed_io::Writer& out_stream, std::mutex& out_file_mutex,
                        std::atomic<bool>& end_of_simulation, size_t task_batch_size) const {
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
      auto tasks = absl::MakeSpan(task_buff.data(), avail_tasks);
      for (const auto& task : tasks) {
        assert(!task.barriers.empty());  // NOLINT
        assert(!task.feats1.empty());    // NOLINT

        if (task.chrom->simulated_size() <= this->deletion_size) {
          continue;
        }
        state = task;  // Set simulation state based on task data
        state.contacts.resize(task.window_end - task.window_start, this->diagonal_width,
                              this->bin_size);

        Simulation::simulate_window(state, out_stream, out_file_mutex);
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
                                 std::mutex& out_stream_mutex) const {
  const auto write_to_stdout = out_stream.path().empty();
  const auto all_barriers = state.barriers;

  size_t last_barrier_deleted_idx{0};

  std::string out_buffer;
  std::string barrier_str_buff;
  do {
    const auto first_chunk = state.active_window_start == state.chrom->start_pos();
    const auto last_chunk = state.active_window_end == state.chrom->end_pos();
    const auto partition_point = (state.active_window_start + state.active_window_end + 1) / 2;

    state.deletion_size =
        std::min(all_barriers[last_barrier_deleted_idx].pos() + 1, this->deletion_size);
    state.deletion_begin = all_barriers[last_barrier_deleted_idx].pos() + 1 - deletion_size;
    const auto deletion_end = state.deletion_begin + state.deletion_size;

    // Generate barrier configuration
    state.barrier_tmp_buff.clear();
    barrier_str_buff.clear();
    // TODO Removeme
    for (const auto& barrier : all_barriers) {
      if (barrier.pos() < state.deletion_begin ||
          barrier.pos() >= state.deletion_begin + state.deletion_size) {
        absl::StrAppend(&barrier_str_buff, barrier.pos(), ",");
      }
    }
    if (!barrier_str_buff.empty()) {
      barrier_str_buff.pop_back();
    }

    std::copy_if(all_barriers.begin(), all_barriers.end(),
                 std::back_inserter(state.barrier_tmp_buff), [&](const auto& barrier) {
                   return barrier.pos() < state.deletion_begin ||
                          barrier.pos() >= state.deletion_begin + state.deletion_size;
                 });

    if (state.barrier_tmp_buff.empty()) {
      continue;
    }

    state.barriers = state.barrier_tmp_buff;
    assert(state.window_end > deletion_size);  // NOLINT
    state.num_lefs =                           // Compute # of LEFs to be simulated
        static_cast<size_t>(std::round(
            this->number_of_lefs_per_mbp *
            (static_cast<double>(state.window_end - deletion_size - state.window_start)) / Mbp));
    // Resize and reset buffers
    state.resize();
    state.reset();
    state.contacts.reset();

    Simulation::simulate_extrusion_kernel(state);

    for (auto i = 0UL; i < state.feats1.size(); ++i) {
      const auto& feat1 = state.feats1[i];
      const auto feat1_abs_center_pos = (feat1.chrom_start + feat1.chrom_end + 1) / 2;
      if (feat1_abs_center_pos >= state.deletion_begin && feat1_abs_center_pos < deletion_end) {
        continue;  // Skip deleted feat1
      }
      for (auto j = i; j < state.feats2.size(); ++j) {
        const auto& feat2 = state.feats2[j];
        const auto feat2_abs_center_pos = (feat2.chrom_start + feat2.chrom_end + 1) / 2;
        if (feat2_abs_center_pos >= state.deletion_begin && feat2_abs_center_pos < deletion_end) {
          continue;  // Skip deleted feat2
        }

        if ((first_chunk || last_chunk) &&
            (feat1_abs_center_pos >= partition_point && feat2_abs_center_pos >= partition_point)) {
          continue;
        }

        const auto feat1_rel_center_pos = feat1_abs_center_pos - state.window_start;
        const auto feat2_rel_center_pos = feat2_abs_center_pos - state.window_start;

        const auto feat1_rel_bin = feat1_rel_center_pos / this->bin_size;
        const auto feat2_rel_bin = feat2_rel_center_pos / this->bin_size;

        const auto contacts = state.contacts.get(feat1_rel_bin, feat2_rel_bin);
        if (contacts == 0) {  // Don't output entries with 0 contacts
          continue;
        }

        const auto feat1_abs_bin = (feat1_abs_center_pos + this->bin_size - 1) / this->bin_size;
        const auto feat2_abs_bin = (feat2_abs_center_pos + this->bin_size - 1) / this->bin_size;

        const auto name = feat1.name.empty() && feat2.name.empty()
                              ? ""
                              : absl::StrCat(feat1.name, ";", feat2.name);
        absl::StrAppend(
            &out_buffer,
            fmt::format(FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tdeletion={}-{}"
                                    "\tnum_barr={}/{}\tbarrs={}\n"),
                        feat1.chrom, feat1_abs_bin * bin_size, (feat1_abs_bin + 1) * bin_size,
                        feat2.chrom, feat2_abs_bin * bin_size, (feat2_abs_bin + 1) * bin_size,
                        name.size() > 1 ? name : "none", contacts, feat1.strand, feat2.strand,
                        state.deletion_begin, deletion_end, state.barriers.size(),
                        all_barriers.size(), barrier_str_buff));
      }
    }
    if (!out_buffer.empty()) {
      std::scoped_lock l(out_stream_mutex);
      if (write_to_stdout) {
        fmt::print(stdout, out_buffer);
      } else {
        out_stream.write(out_buffer);
      }
      out_buffer.clear();
    }
    /*
    if (foo == k++) {
      auto c = cooler::Cooler(fmt::format("/tmp/test_{}_{}_{:03}.cool", state.id, state.cell_id,
    foo), cooler::Cooler::WRITE_ONLY, this->bin_size, state.chrom.name().size());
      c.write_or_append_cmatrix_to_file(state.contacts, state.chrom.name(),
                                        state.chrom.start_pos(), state.chrom.end_pos(),
                                        state.chrom.size());
    }
     */
  } while (++last_barrier_deleted_idx < all_barriers.size());

  fmt::print(stderr, FMT_STRING("Done processing {}[{}-{}]!\n"), state.chrom->name(),
             state.active_window_start, state.active_window_end);
}

}  // namespace modle
