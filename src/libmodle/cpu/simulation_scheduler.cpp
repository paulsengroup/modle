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

#include "modle/common/config.hpp"  // for Config
#include "modle/common/utils.hpp"   // for traced
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

  auto query_interval_tree = [&](const auto& tree, const auto start, const auto end) {
    const auto [overlap_begin, overlap_end] = tree.find_overlaps(start, end);
    if (overlap_begin == overlap_end) {
      return absl::MakeConstSpan(&(*overlap_begin), 0UL);
    }
    return absl::MakeConstSpan(&(*overlap_begin), static_cast<size_t>(overlap_end - overlap_begin));
  };

  using Task_t = std::pair<const bed::BED*, Simulation::TaskPW>;

  constexpr size_t task_batch_size_enq = 32;  // NOLINTNEXTLINE
  moodycamel::BlockingConcurrentQueue<Task_t> task_queue(this->nthreads * task_batch_size_enq, 1,
                                                         0);
  moodycamel::ProducerToken ptok(task_queue);
  std::array<Task_t, task_batch_size_enq> tasks;

  std::ifstream out_file;  // TODO
  std::mutex out_file_mutex;
  std::atomic<bool> end_of_simulation = false;

  auto tpool = Simulation::instantiate_thread_pool();
  for (auto i = 0UL; i < this->nthreads; ++i) {  // Start simulation threads
    boost::asio::post(
        tpool, [&]() { this->worker(task_queue, out_file, out_file_mutex, end_of_simulation); });
  }

  size_t task_id = 0;
  size_t num_tasks = 0;
  for (auto& chrom : this->_genome) {
    size_t cell_id = 0;
    const auto& reference_features = chrom.get_features();
    const auto& barriers = chrom.barriers();
    if (reference_features.size() != 2 || barriers.empty()) {
      continue;
    }

    // We assume the first set of feature is the promoter/TSS, and we look around using window of
    // size ~2x diagonal_width. We want to process windows with at least 1 enhancer and 1 extr.
    // barrier
    for (const auto& feature : reference_features[0].data()) {
      const auto feature_center_pos = (feature.chrom_start + feature.chrom_end + 1) / 2;

      const auto range_start = feature_center_pos > this->diagonal_width
                                   ? feature_center_pos - this->diagonal_width
                                   : 0UL;
      const auto range_end = std::min(feature_center_pos + this->diagonal_width, chrom.end_pos());

      const auto overlapping_barriers = query_interval_tree(barriers, range_start, range_end);
      if (overlapping_barriers.empty()) {
        break;
      }

      const auto overlapping_features =
          query_interval_tree(reference_features[1], range_start, range_end);
      if (overlapping_features.empty()) {
        break;
      }
      TaskPW t{{task_id++, &chrom, cell_id++}};

      if (this->target_contact_density != 0) {  // Compute the number of simulation rounds required
        // to reach the target contact density
        t.num_target_contacts = static_cast<size_t>(std::max(
            1.0,
            std::round(this->target_contact_density *
                       static_cast<double>((range_end - range_start) * this->diagonal_width))));
      }
      t.num_target_epochs = t.num_target_contacts == 0UL ? this->simulation_iterations
                                                         : std::numeric_limits<size_t>::max();

      t.num_lefs = static_cast<size_t>(static_cast<double>(range_end - range_start) *
                                       this->number_of_lefs_per_mbp);
      t.range_start = range_start;
      t.range_end = range_end;

      t.barriers = overlapping_barriers;
      t.features = overlapping_features;

      if (num_tasks == tasks.size()) {
        auto sleep_us = 100;  // NOLINT
        while (               // Enqueue tasks
            !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), num_tasks)) {
          sleep_us = std::min(100000, sleep_us * 2);  // NOLINT
          std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
        }
        num_tasks = 0;
      }
      tasks[num_tasks++] = std::make_pair(&feature, t);

      // Put overlaps on queue
      //
    }
  }
  if (num_tasks != 0) {
    while (  // Enqueue tasks
        !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), num_tasks)) {
      std::this_thread::sleep_for(std::chrono::microseconds(100));  // NOLINT
    }
  }
}

void Simulation::worker(
    moodycamel::BlockingConcurrentQueue<std::pair<const bed::BED*, Simulation::TaskPW>>& task_queue,
    std::ifstream& out_file, std::mutex& out_file_mutex, std::atomic<bool>& end_of_simulation,
    size_t task_batch_size) const {
  fmt::print(stderr, FMT_STRING("Spawning simulation thread {}...\n"), std::this_thread::get_id());
  moodycamel::ConsumerToken ctok(task_queue);

  using Task_t = std::pair<const bed::BED*, Simulation::TaskPW>;
  absl::FixedArray<Task_t> task_buff(task_batch_size);  // Tasks are dequeue in batch.

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
      for (auto& [feature_ptr, task] : tasks) {
        assert(!task.barriers.empty());  // NOLINT
        assert(!task.features.empty());  // NOLINT

        if (task.chrom->simulated_size() <= this->deletion_size) {
          continue;
        }
        state = task;  // Set simulation state based on task data
        state.contacts.resize(task.range_end - task.range_start, this->diagonal_width,
                              this->bin_size);

        Simulation::simulate_window(*feature_ptr, state, out_file, out_file_mutex);
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

void Simulation::simulate_window(const bed::BED& reference_feature, Simulation::StatePW& state,
                                 std::ifstream& out_file, std::mutex& out_file_mutex) const {
  (void)out_file;
  const auto all_barriers = state.barriers;
  const auto full_range_start = state.range_start;
  const auto full_range_end = state.range_end;

  size_t last_barrier_deleted_idx{0};

  do {
    state.deletion_size =
        std::min(all_barriers[last_barrier_deleted_idx].pos() + 1, this->deletion_size);
    state.deletion_begin = all_barriers[last_barrier_deleted_idx].pos() + 1 - deletion_size;
    const auto deletion_end = state.deletion_begin + state.deletion_size;

    if (const auto pos = (reference_feature.chrom_start + reference_feature.chrom_end + 1) / 2;
        pos >= state.deletion_begin && pos < deletion_end) {
      continue;  // "Reference" reference_feature falls withing the deletion
    }

    // TODO: all these comparisons should be made using the bin to which the middle of the
    // reference_feature maps
    if (state.features.front().chrom_start >= state.deletion_begin &&
        state.features.back().chrom_start < deletion_end) {
      continue;  // Is it safe to break here?
    }

    [&]() {  // Generate barrier configuration
      state.barrier_tmp_buff.clear();
      for (const auto& barrier : all_barriers) {
        if (barrier.pos() < state.deletion_begin) {
          state.barrier_tmp_buff.push_back(barrier);
        } else if (barrier.pos() >= deletion_end) {
          state.barrier_tmp_buff.push_back(barrier);
          state.barrier_tmp_buff.back().pos() -= state.deletion_size;
        }
      }
    }();

    if (state.barrier_tmp_buff.empty()) {
      continue;
    }

    state.num_target_contacts = static_cast<size_t>(this->target_contact_density *
                                                    static_cast<double>(state.contacts.npixels()));
    state.barriers = state.barrier_tmp_buff;
    assert(state.range_end > deletion_size);  // NOLINT
    state.range_end -= deletion_size;
    state.num_lefs =  // Compute # of LEFs to be simulated
        static_cast<size_t>(std::round(this->number_of_lefs_per_mbp *
                                       (static_cast<double>(state.range_end - state.range_start)) /
                                       Mbp));
    // Resize and reset buffers
    state.resize();
    state.reset();
    state.contacts.reset();

    // fmt::print(stderr, "{}[{}-{}]\n", state.chrom->name(), state.range_start, state.range_end);
    Simulation::simulate_extrusion_kernel(state);

    const auto reference_feature_rel_bin =
        (reference_feature.chrom_end - reference_feature.chrom_start + this->bin_size - 1) /
        this->bin_size;

    for (const auto& target_feature : state.features) {
      // Skip deleted features
      if (target_feature.chrom_start >= state.deletion_begin &&
          target_feature.chrom_end < deletion_end) {
        continue;
      }

      const auto target_feature_rel_bin =
          (target_feature.chrom_end - target_feature.chrom_start + this->bin_size - 1) /
          this->bin_size;

      const auto contacts = state.contacts.get(reference_feature_rel_bin, target_feature_rel_bin);
      if (contacts == 0) {  // Don't output entries with 0 contacts
        continue;
      }

      const auto reference_feature_abs_bin =
          (reference_feature.chrom_start + reference_feature.chrom_start + this->bin_size - 1) /
          this->bin_size;

      const auto target_feature_abs_bin =
          (target_feature.chrom_start + target_feature.chrom_start + this->bin_size - 1) /
          this->bin_size;

      const auto name = absl::StrCat(reference_feature.name, ";", target_feature.name);
      fmt::print(stdout,
                 FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcenter={}\tdeletion={}-{}"
                            "\tnum_barr={}/{}\tnum_feats={}\n"),
                 reference_feature.chrom, reference_feature_abs_bin * bin_size,
                 (reference_feature_abs_bin + 1) * bin_size, target_feature.chrom,
                 target_feature_abs_bin * bin_size, (target_feature_abs_bin + 1) * bin_size,
                 name.size() > 1 ? name : "none", contacts, reference_feature.strand,
                 target_feature.strand, ((2 * reference_feature_abs_bin * bin_size) + bin_size) / 2,
                 state.deletion_begin, deletion_end, state.barriers.size(), all_barriers.size(),
                 state.features.size());
    }

    // auto c = cooler::Cooler(fmt::format("/tmp/test_{}_{}.cool", state.id, state.cell_id),
    //                         cooler::Cooler::WRITE_ONLY, this->bin_size);
    // c.write_or_append_cmatrix_to_file(state.contacts, state.chrom->name(),
    // state.chrom->start_pos(),
    //                                   state.chrom->end_pos(), state.chrom->size());
  } while (++last_barrier_deleted_idx < all_barriers.size());

  std::scoped_lock l(out_file_mutex);
  fmt::print(stderr, FMT_STRING("Done processing {}[{}-{}]!\n"), state.chrom->name(),
             full_range_start, full_range_end);
}

}  // namespace modle
