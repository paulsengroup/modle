#pragma once

// IWYU pragma: private, include "modle/simulation.hpp"

#include <absl/container/btree_set.h>            // for btree_multiset
#include <absl/container/flat_hash_map.h>        // for flat_hash_map
#include <absl/meta/type_traits.h>               // for remove_reference_t
#include <absl/strings/str_join.h>               // for StrJoin
#include <absl/types/span.h>                     // for Span, MakeSpan
#include <bits/exception.h>                      // for exception
#include <cpp-sort/sorter_facade.h>              // for sorter_facade
#include <cpp-sort/sorters/counting_sorter.h>    // for counting_sort
#include <cpp-sort/sorters/drop_merge_sorter.h>  // for drop_merge_sort
#include <cpp-sort/sorters/insertion_sorter.h>   // for insertion_sort
#include <cpp-sort/sorters/pdq_sorter.h>         // for pdq_sort
#include <cpp-sort/sorters/ska_sorter.h>         // for ska_sort, ska_so...
#include <fmt/format.h>                          // for print, format
#include <fmt/ostream.h>                         // for formatbuf<>::int...
#include <moodycamel/blockingconcurrentqueue.h>  // for BlockingConcurre...
#include <moodycamel/concurrentqueue.h>          // for ConsumerToken

#include <algorithm>                                // for max, fill, is_so...
#include <array>                                    // for array
#include <atomic>                                   // for atomic
#include <boost/asio/post.hpp>                      // IWYU pragma: keep for post
#include <boost/asio/thread_pool.hpp>               // for thread_pool
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/range/adaptor/reversed.hpp>         // for adaptors::reverse
#include <cassert>                                  // for assert
#include <chrono>                                   // for milliseconds
#include <cmath>                                    // for round
#include <cstddef>                                  // IWYU pragma: keep for size_t
#include <cstdint>                                  // for uint64_t, uint32_t
#include <cstdio>                                   // for stderr
#include <deque>                                    // for deque, _Deque_it...
#include <filesystem>                               // for operator<<, path
#include <functional>                               // for hash
#include <iterator>                                 // for move_iterator
#include <limits>                                   // for numeric_limits
#include <memory>                                   // for unique_ptr, make...
#include <mutex>                                    // for mutex
#include <new>                                      // for operator new
#include <numeric>                                  // for accumulate, iota
#include <random>                                   // for uniform_int_dist...
#include <sstream>                                  // for basic_stringbuf<...
#include <stdexcept>                                // for runtime_error
#include <string>                                   // for string, basic_st...
#include <string_view>                              // for hash, string_view
#include <thread>                                   // for operator<<, get_id
#include <type_traits>                              // for declval, decay_t
#include <utility>                                  // for pair, addressof
#include <vector>                                   // for vector

#include "modle/bed.hpp"                         // for BED, Parser
#include "modle/chrom_sizes.hpp"                 // for ChromSize, Parser
#include "modle/common.hpp"                      // for bp_t, fwd, rev, PRNG
#include "modle/config.hpp"                      // for Config
#include "modle/contacts.hpp"                    // for ContactMatrix
#include "modle/cooler.hpp"                      // for Cooler, Cooler::...
#include "modle/dna.hpp"                         // for Chromosome
#include "modle/extrusion_barriers.hpp"          // for ExtrusionBarrier
#include "modle/extrusion_factors.hpp"           // for Lef, ExtrusionUnit
#include "modle/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP
#include "modle/utils.hpp"                       // for ndebug_defined, traced

#ifndef BOOST_STACKTRACE_USE_NOOP
#include <boost/exception/get_error_info.hpp>  // for get_error_info
#include <boost/stacktrace/stacktrace.hpp>     // for operator<<
#include <ios>                                 // IWYU pragma: keep for streamsize
#endif

namespace modle {

Simulation::Simulation(const Config& c, bool import_chroms)
    : Config(c),
      _chromosomes(import_chroms
                       ? import_chromosomes(path_to_chrom_sizes, path_to_extr_barriers,
                                            path_to_chrom_subranges, write_contacts_for_ko_chroms)
                       : Genome{}) {}

std::size_t Simulation::size() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
      [](auto accumulator, const auto& chrom) { return accumulator + chrom.size(); });
}

std::size_t Simulation::simulated_size() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
                         [](auto accumulator, const auto& chrom) {
                           return accumulator + (chrom.end_pos() - chrom.start_pos());
                         });
}

// TODO Add flag to import skip chrom without barriers
Simulation::Genome Simulation::import_chromosomes(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_extr_barriers,
    const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms) {
  assert(!path_to_chrom_sizes.empty());    // NOLINT
  assert(!path_to_extr_barriers.empty());  // NOLINT
  Genome chroms;

  // Parse chrom subranges from BED. We parse everything at once to deal with duplicate entries
  absl::flat_hash_map<std::string, std::pair<uint64_t, uint64_t>> chrom_ranges;
  if (!path_to_chrom_subranges.empty()) {
    for (auto&& record : modle::bed::Parser(path_to_chrom_subranges).parse_all()) {
      chrom_ranges.emplace(std::move(record.chrom),
                           std::make_pair(record.chrom_start, record.chrom_end));
    }
  }

  // Parse chrom. sizes and build the set of chromosome to be simulated.
  // When the BED file with the chrom. subranges is not provided, all the chromosome in the
  // chrom.sizes file will be selected and returned. When a BED file with the chrom. subranges is
  // available, then only chromosomes that are present in both files will be selected. Furthermore
  // we are also checking that the subrange lies within the coordinates specified in the chrom.
  // sizes file
  for (auto&& chrom : modle::chrom_sizes::Parser(path_to_chrom_sizes).parse_all()) {
    if (auto match = chrom_ranges.find(chrom.name); match != chrom_ranges.end()) {
      const auto& range_start = match->second.first;
      const auto& range_end = match->second.second;
      if (range_start < chrom.start || range_end > chrom.end) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("According to the chrom.sizes file {}, chromosome '{}' should have a size "
                       "of '{}', but the subrange specified in BED file {} extends past this "
                       "region: range {}:{}-{} does not fit in range {}:{}-{}"),
            path_to_chrom_sizes, chrom.name, chrom.end, path_to_chrom_subranges, chrom.name,
            range_start, range_end, chrom.name, chrom.start, chrom.end));
      }
      chrom.start = range_start;
      chrom.end = range_end;
      chroms.emplace(std::move(chrom));
    } else if (chrom_ranges.empty() || keep_all_chroms) {
      chroms.emplace(std::move(chrom));
    }
  }

  // Parse all the records from the BED file. parse_all() will throw in case of duplicates.
  // This for loop selects extrusion barriers that fall within the chromosomes to be simulated
  for (auto&& record : modle::bed::Parser(path_to_extr_barriers).parse_all()) {
    if (record.score < 0 || record.score > 1) {
      throw std::runtime_error(
          fmt::format("Invalid score field detected for record {}[{}-{}]: expected a score "
                      "between 0 and 1, got {:.4g}.",
                      record.chrom, record.chrom_start, record.chrom_end, record.score));
    }
    if (!chrom_ranges.empty() && chrom_ranges.find(record.chrom) == chrom_ranges.end()) {
      continue;
    }

    if (auto match = chroms.find(record.chrom); match != chroms.end()) {
      match->add_extrusion_barrier(record);
    }
  }
  return chroms;
}

std::vector<ExtrusionBarrier> Simulation::allocate_barriers(const Chromosome* const chrom) {
  std::vector<ExtrusionBarrier> barriers;
  std::size_t barriers_skipped = 0;
  for (const auto& b : chrom->get_barriers()) {
    if (b.strand == '+' || b.strand == '-') MODLE_LIKELY {
        const auto pos = (b.chrom_start + b.chrom_end + 1) / 2;
        if (pos < chrom->start_pos() || pos >= chrom->end_pos()) {
          continue;
        }
        if (b.score != 0) {
          const auto pblock = b.score;
          const auto pno = this->ctcf_not_occupied_self_prob;
          const auto poo =
              ExtrusionBarrier::compute_blocking_to_blocking_transition_probabilities_from_pblock(
                  pblock, pno);
          barriers.emplace_back(pos, poo, pno, b.strand);
        } else {
          barriers.emplace_back(pos, this->ctcf_occupied_self_prob,
                                this->ctcf_not_occupied_self_prob, b.strand);
        }
      }
    else {
      ++barriers_skipped;
    }
  }

  fmt::print(stderr,
             FMT_STRING("Instantiated {} extr. barriers for '{}' ({} barriers were skipped).\n"),
             barriers.size(), chrom->name(), barriers_skipped);

  return barriers;
}

void Simulation::run() {
  if (!this->skip_output) {  // Write simulation params to file
    if (this->force) {
      std::filesystem::remove_all(this->path_to_output_file);
    }
    std::filesystem::create_directories(this->path_to_output_file.parent_path());
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
  std::deque<std::pair<Chromosome*, std::size_t>> progress_queue;

  // Barriers are shared across threads that are simulating loop-extrusion on a given chromosome.
  // Extrusion barriers are placed on this std::deque, and are passed to simulation threads
  // through a Span of const ExtrusionBarriers
  std::mutex barrier_mutex;
  absl::flat_hash_map<Chromosome*, std::unique_ptr<std::vector<ExtrusionBarrier>>> barriers;

  constexpr std::size_t task_batch_size = 128;
  // Queue used to submit simulation tasks to the thread pool
  moodycamel::BlockingConcurrentQueue<Simulation::Task> task_queue(
      std::max(static_cast<std::size_t>(static_cast<double>(ncells) * 1.25),
               this->nthreads * task_batch_size),
      1, 0);
  moodycamel::ProducerToken ptok(task_queue);

  boost::asio::post(tpool, [&]() {  // This thread is in charge of writing contacts to disk
    Chromosome* chrom_to_be_written = nullptr;
    const auto max_str_length =
        std::max_element(  // Find chrom with the longest name
            this->_chromosomes.begin(), this->_chromosomes.end(),
            [](const auto& c1, const auto& c2) { return c1.name().size() < c2.name().size(); })
            ->name()
            .size();

    auto c = this->path_to_output_file.empty()
                 ? nullptr
                 : std::make_unique<cooler::Cooler>(this->path_to_output_file,
                                                    cooler::Cooler::WRITE_ONLY, this->bin_size,
                                                    max_str_length);
    while (true) {
      std::this_thread::sleep_for(std::chrono::milliseconds(25));
      {
        std::scoped_lock l(progress_queue_mutex);
        if (progress_queue.empty()) {
          continue;
        }
        if (auto& [chrom, count] = progress_queue.front(); chrom == nullptr) MODLE_UNLIKELY {
            end_of_simulation = true;
            return;
          }
        else if (count == ncells) {
          chrom_to_be_written = chrom;
          progress_queue.pop_front();
        } else
          MODLE_LIKELY {
            assert(count < ncells);  // NOLINT
            continue;
          }
      }
      try {
        if (c) {  // c == nullptr only when --skip-output is used
          fmt::print(stderr, "Writing contacts for '{}' to file {}...\n",
                     chrom_to_be_written->name(), c->get_path());
          c->write_or_append_cmatrix_to_file(
              chrom_to_be_written->contacts(), chrom_to_be_written->name(),
              chrom_to_be_written->start_pos(), chrom_to_be_written->end_pos(),
              chrom_to_be_written->size(), true);
          fmt::print(stderr, "Written {} contacts for '{}' in {:.2f}M pixels to file {}.\n",
                     chrom_to_be_written->contacts().get_tot_contacts(),
                     chrom_to_be_written->name(),
                     static_cast<double>(chrom_to_be_written->contacts().npixels()) / 1.0e6,
                     c->get_path());
        }
        chrom_to_be_written->deallocate_contacts();
      } catch (const std::runtime_error& err) {
        throw std::runtime_error(fmt::format(
            FMT_STRING(
                "The following error occurred while writing contacts for '{}' to file {}: {}"),
            chrom_to_be_written->name(), c->get_path(), err.what()));
      }
    }
  });

  for (auto i = 0UL; i < this->nthreads; ++i) {  // Start simulation threads
    moodycamel::ConsumerToken ctok(task_queue);
    boost::asio::post(tpool, [&, ctok = std::move(ctok)]() mutable {
      fmt::print(stderr, FMT_STRING("Spawning simulation thread {}...\n"),
                 std::this_thread::get_id());
      try {
        std::array<Task, task_batch_size> task_buff;  // Tasks are dequeue in batch.
                                                      // This is to reduce contention when
                                                      // accessing the queue

        // This state object owns all the buffers and PRNG + seed required in order to simulate loop
        // extrusion for a single cell. Buffers are allocated once and resized, cleared and reused
        // throughout the simulation
        Simulation::State state_buff;

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
            state_buff.barrier_mask.resize(task.barriers.size());
            state_buff.reset();  // Clear all buffers

            if (task.cell_id == 0) {
              // Print a status update when we are processing cell #0 for a given chromosome
              const auto target_epochs = static_cast<std::size_t>(
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
              fmt::print(stderr, "Simulation for '{}' successfully completed.\n",
                         task.chrom->name());
              std::scoped_lock l2(barrier_mutex);
              barriers.erase(task.chrom);
            }
          }
        }
      } catch (const std::exception& err) {
        // This is needed, as exceptions don't seem to always propagate to the main()
        // TODO: Find a better way to propagate exception up to the main thread. Maybe use a
        // concurrent queue?
        fmt::print(stderr, FMT_STRING("Detected an error in thread {}:\n{}\n"),
                   std::this_thread::get_id(), err.what());
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
    });
  }

  // The remaining code submits simulation tasks to the queue. Then it waits until all the tasks
  // have been completed, and contacts have been written to disk

  // Create a vector of pointers to chroms, and sort the pointers by chrom size
  // TODO Update this to process chroms in the same order they are found in the chrom.sizes file
  std::vector<Chromosome*> chroms(static_cast<std::size_t>(this->_chromosomes.size()));
  std::transform(this->_chromosomes.begin(), this->_chromosomes.end(), chroms.begin(),
                 [](auto& chrom) { return &chrom; });
  std::sort(chroms.begin(), chroms.end(),
            [](const auto* const c1, const auto* const c2) { return c1->size() > c2->size(); });

  absl::Span<const ExtrusionBarrier> extr_barriers_buff{};
  std::vector<Task> tasks(ncells);
  auto taskid = 0UL;

  // Loop over chromosomes
  for (auto* chrom : chroms) {
    // Don't simulate KO chroms (but write them to disk if the user desires so)
    if (!chrom->ok()) {
      fmt::print(stderr, "SKIPPING '{}'...\n", chrom->name());
      if (this->write_contacts_for_ko_chroms) {
        chrom->allocate_contacts(this->bin_size, this->diagonal_width);
        std::scoped_lock l(barrier_mutex, progress_queue_mutex);
        progress_queue.emplace_back(chrom, ncells);
        barriers.emplace(chrom, nullptr);
      }
      continue;
    }

    // Allocate the contact matrix. Once it's not needed anymore, the contact matrix will be
    // de-allocated by the thread that is writing contacts to disk
    chrom->allocate_contacts(this->bin_size, this->diagonal_width);
    {
      // For consistency, it is important that both locks are held while new items are added to the
      // two queues
      std::scoped_lock l(barrier_mutex, progress_queue_mutex);

      // Signal that we have started processing the current chrom
      progress_queue.emplace_back(chrom, 0UL);

      // Allocate extr. barriers mapping on the chromosome that is being simulated
      auto node = barriers.emplace(chrom, std::make_unique<std::vector<ExtrusionBarrier>>(
                                              Simulation::allocate_barriers(chrom)));
      extr_barriers_buff = absl::MakeConstSpan(*node.first->second);
    }

    // Compute # of LEFs to be simulated based on chrom. sizes
    const auto nlefs = static_cast<std::size_t>(std::round(
        this->number_of_lefs_per_mbp * (static_cast<double>(chrom->simulated_size()) / 1.0e6)));

    // This resize is not needed at the moment. It's there as a safe-guard just in case we decide to
    tasks.resize(ncells);  // allow the number of simulated cells to change across different chroms

    auto target_contacts = 0UL;
    if (this->target_contact_density != 0) {  // Compute the number of simulation rounds required
                                              // to reach the target contact density
      target_contacts = static_cast<std::size_t>(
          std::max(1.0, std::round((this->target_contact_density *
                                    static_cast<double>(chrom->contacts().npixels())) /
                                   static_cast<double>(this->ncells))));
    }
    auto target_epochs = target_contacts == 0UL ? this->simulation_iterations
                                                : std::numeric_limits<std::size_t>::max();

    // Generate a batch of tasks for all the simulations involving the current chrom
    std::generate(tasks.begin(), tasks.end(), [&, cellid = 0UL]() mutable {
      return Task{taskid++,        chrom, cellid++,          target_epochs,
                  target_contacts, nlefs, extr_barriers_buff};
    });

    while (  // Enqueue tasks
        !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), tasks.size())) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  }

  {  // Signal end of simulation to the thread that is writing contacts to disk
    std::scoped_lock l(progress_queue_mutex);
    progress_queue.emplace_back(nullptr, 0UL);
  }
  tpool.join();
  assert(end_of_simulation);  // NOLINT
}

void Simulation::simulate_extrusion_kernel(Simulation::State& s) {
  assert(s.nlefs == s.lef_buff.size());     // NOLINT
  assert(s.nlefs == s.rank_buff1.size());   // NOLINT
  assert(s.nlefs == s.rank_buff2.size());   // NOLINT
  assert(s.nlefs == s.moves_buff1.size());  // NOLINT
  assert(s.nlefs == s.moves_buff2.size());  // NOLINT
  assert(s.nlefs == s.idx_buff1.size());    // NOLINT
  assert(s.nlefs == s.idx_buff2.size());    // NOLINT
  assert(s.nlefs == s.epoch_buff.size());   // NOLINT

  // Seed is computed based on chrom. name, size and cellid
  s.seed = this->seed + std::hash<std::string_view>{}(s.chrom->name()) +
           std::hash<std::size_t>{}(s.chrom->size()) + std::hash<std::size_t>{}(s.cell_id);
  s.rand_eng = modle::PRNG(s.seed);

  // Generate the epoch at which each LEF is supposed to be initially loaded
  auto& lef_initial_loading_epoch = s.epoch_buff;
  lef_initial_loading_epoch.resize(this->skip_burnin ? 0 : s.nlefs);

  if (!this->skip_burnin) {
    // TODO Consider using a Poisson process instead of sampling from an uniform distribution
    std::uniform_int_distribution<std::size_t> round_gen{
        0, (4 * this->average_lef_lifetime) / this->bin_size};
    std::generate(lef_initial_loading_epoch.begin(), lef_initial_loading_epoch.end(),
                  [&]() { return round_gen(s.rand_eng); });

    // Sort epochs in descending order
    if (round_gen.max() > 2048) {
      // Counting sort uses n + r space in memory, where r is the number of unique values in the
      // range to be sorted. For this reason it is not a good idea to use it when the sampling
      // interval is relatively large. Whether 2048 is a reasonable threshold has yet to be tested
      cppsort::ska_sort(lef_initial_loading_epoch.rbegin(), lef_initial_loading_epoch.rend());
    } else {
      cppsort::counting_sort(lef_initial_loading_epoch.rbegin(), lef_initial_loading_epoch.rend());
    }
  }

  // Shift epochs so that the first epoch == 0
  if (const auto offset = lef_initial_loading_epoch.back(); offset != 0) {
    std::for_each(lef_initial_loading_epoch.begin(), lef_initial_loading_epoch.end(),
                  [&](auto& n) { n -= offset; });
  }

  // The highest loading epoch equals to the number of burnin epochs
  const auto n_burnin_epochs = this->skip_burnin ? 0 : lef_initial_loading_epoch.front();
  s.n_target_epochs += s.n_target_contacts != 0 ? 0 : n_burnin_epochs;

  // Compute the avg. # of LEFs to use to sample contact every iterations
  const auto avg_nlefs_to_sample =
      static_cast<double>(s.nlefs) * this->lef_fraction_contact_sampling;

  // Compute the avg. # of LEFs to be released every iterations
  const auto avg_nlefs_to_release =
      static_cast<double>((this->rev_extrusion_speed + this->fwd_extrusion_speed) * s.nlefs) /
      static_cast<double>(this->average_lef_lifetime);

  // Local counter for the number of contacts generated by the current kernel instance
  auto n_contacts = 0UL;

  assert(avg_nlefs_to_sample <= static_cast<double>(s.nlefs));   // NOLINT
  assert(avg_nlefs_to_release <= static_cast<double>(s.nlefs));  // NOLINT
  assert((s.n_target_contacts != 0) ==                           // NOLINT
         (s.n_target_epochs == std::numeric_limits<std::size_t>::max()));

  // Declare the spans used by the burnin phase as well as the simulation itself
  auto lefs = absl::MakeSpan(s.lef_buff);
  auto lef_unloader_affinity = absl::MakeSpan(s.lef_unloader_affinity);
  const auto barriers = absl::MakeConstSpan(s.barriers);
  auto rev_lef_ranks = absl::MakeSpan(s.rank_buff1);
  auto fwd_lef_ranks = absl::MakeSpan(s.rank_buff2);
  auto rev_moves = absl::MakeSpan(s.moves_buff1);
  auto fwd_moves = absl::MakeSpan(s.moves_buff2);
  auto rev_collision_mask = absl::MakeSpan(s.idx_buff1);
  auto fwd_collision_mask = absl::MakeSpan(s.idx_buff2);

  // Generate initial extr. barrier states, so that they are already at or close to equilibrium
  for (auto i = 0UL; i < s.barrier_mask.size(); ++i) {
    s.barrier_mask[i] =
        std::bernoulli_distribution{this->probability_of_extrusion_barrier_block}(s.rand_eng);
  }

  std::size_t nlefs_to_release;

  // Start the burnin phase (followed by the actual simulation)
  for (auto epoch = 0UL; epoch < s.n_target_epochs; ++epoch) {
    if (!lef_initial_loading_epoch.empty()) {  // Execute this branch only during the burnin phase
      if (this->skip_burnin) {
        // The following will cause every LEF to be loaded in the current iteration
        lef_initial_loading_epoch.clear();
      }

      // Consume epochs for LEFs that are supposed to be loaded in the current epoch
      while (!lef_initial_loading_epoch.empty() && lef_initial_loading_epoch.back() == epoch) {
        lef_initial_loading_epoch.pop_back();
      }
      // Don't remove this assertion! It is very useful to prevent empty Spans (which can cause
      // weird issues later in the simulation) NOLINTNEXTLINE
      assert(lef_initial_loading_epoch.size() < s.nlefs);  // NOLINT

      // Compute the current number of active LEFs (i.e. LEFs that are either bound to DNA, or
      // that are valid candidates for re-binding). Grow spans accordingly
      const auto nlefs = s.nlefs - lef_initial_loading_epoch.size();

      lefs = absl::MakeSpan(s.lef_buff.data(), nlefs);
      lef_unloader_affinity = absl::MakeSpan(s.lef_unloader_affinity.data(), nlefs);
      rev_lef_ranks = absl::MakeSpan(s.rank_buff1.data(), nlefs);
      fwd_lef_ranks = absl::MakeSpan(s.rank_buff2.data(), nlefs);
      rev_moves = absl::MakeSpan(s.moves_buff1.data(), nlefs);
      fwd_moves = absl::MakeSpan(s.moves_buff2.data(), nlefs);
      rev_collision_mask = absl::MakeSpan(s.idx_buff1.data(), nlefs);
      fwd_collision_mask = absl::MakeSpan(s.idx_buff2.data(), nlefs);

      // Sample nlefs to be released while in burn-in phase
      nlefs_to_release = std::min(
          lefs.size(), std::poisson_distribution<std::size_t>{
                           static_cast<double>(
                               (this->rev_extrusion_speed + this->fwd_extrusion_speed) * s.nlefs) /
                           static_cast<double>(this->average_lef_lifetime)}(s.rand_eng));
    } else {  // Sample nlefs to be released after the burn-in phase has been completed
      nlefs_to_release = std::min(
          lefs.size(), std::poisson_distribution<std::size_t>{avg_nlefs_to_release}(s.rand_eng));
    }
    ////////////////////////
    // Simulate one epoch //
    ////////////////////////

    {  // Select and bind LEFs
      auto lef_mask = absl::MakeSpan(s.idx_buff1.data(), lefs.size());
      this->select_lefs_to_bind(lefs, lef_mask);
      this->bind_lefs(s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, lef_mask, s.rand_eng,
                      epoch == 0);
    }

    if (epoch > n_burnin_epochs) {              // Register contacts
      assert(fwd_lef_ranks.size() == s.nlefs);  // NOLINT

      auto nlefs_to_sample = std::min(
          lefs.size(), std::poisson_distribution<std::size_t>{avg_nlefs_to_sample}(s.rand_eng));

      if (s.n_target_contacts != 0) {  // When using the target contact density as stopping
                                       // criterion, don't overshoot the target number of contacts
        nlefs_to_sample = std::min(nlefs_to_sample, s.n_target_contacts - n_contacts);
      }

      // Select LEFs to be used for contact registration.
      // We are sampling from fwd ranks to avoid having to allocate a vector of indices just to do
      // this sampling
      auto lef_idx = absl::MakeSpan(s.idx_buff1.data(), nlefs_to_sample);
      std::sample(fwd_lef_ranks.begin(), fwd_lef_ranks.end(), lef_idx.begin(), lef_idx.size(),
                  s.rand_eng);
      n_contacts += this->register_contacts(s.chrom, lefs, lef_idx);

      if (s.n_target_contacts != 0 && n_contacts >= s.n_target_contacts) {
        return;  // Enough contact have been generated. Yay!
      }
    }

    this->generate_moves(s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                         s.rand_eng);

    this->generate_ctcf_states(barriers, s.barrier_mask, s.rand_eng);

    // Reset collision masks
    std::fill(rev_collision_mask.begin(), rev_collision_mask.end(), NO_COLLISION);
    std::fill(fwd_collision_mask.begin(), fwd_collision_mask.end(), NO_COLLISION);

    const auto& [nrev_units_at_5prime, nfwd_units_at_3prime] = Simulation::process_collisions(
        s.chrom, lefs, barriers, s.barrier_mask, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
        rev_collision_mask, fwd_collision_mask, s.rand_eng);

    this->extrude(s.chrom, lefs, rev_moves, fwd_moves, nrev_units_at_5prime, nfwd_units_at_3prime);

    this->generate_lef_unloader_affinities(lefs, barriers, absl::MakeConstSpan(rev_collision_mask),
                                           absl::MakeConstSpan(fwd_collision_mask),
                                           lef_unloader_affinity);

    // Reusing this buffer is ok, as at this point we don't need access to collision information
    auto lef_idx = absl::MakeSpan(s.idx_buff1.data(), nlefs_to_release);
    this->select_lefs_to_release(lef_idx, lef_unloader_affinity, s.rand_eng);
    this->release_lefs(lefs, lef_idx);
  }
}

template <typename MaskT>
void Simulation::bind_lefs(const Chromosome* chrom, absl::Span<Lef> lefs,
                           absl::Span<std::size_t> rev_lef_ranks,
                           absl::Span<std::size_t> fwd_lef_ranks, MaskT& mask,
                           modle::PRNG_t& rand_eng,
                           bool first_epoch) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_same_v<boost::dynamic_bitset<>, MaskT> ||
                    std::is_integral_v<std::decay_t<decltype(std::declval<MaskT&>().operator[](
                        std::declval<std::size_t>()))>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");

  assert(lefs.size() <= mask.size() || mask.empty());  // NOLINT
  assert(chrom);                                       // NOLINT
  chrom_pos_generator_t pos_generator{chrom->start_pos(), chrom->end_pos() - 1};
  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (mask.empty() || mask[i]) {  // Bind all LEFs when mask is empty
      lefs[i].bind_at_pos(pos_generator(rand_eng));
    }
  }

  this->rank_lefs(lefs, rev_lef_ranks, fwd_lef_ranks, !first_epoch);
}

void Simulation::generate_ctcf_states(absl::Span<const ExtrusionBarrier> extr_barriers,
                                      boost::dynamic_bitset<>& mask,
                                      modle::PRNG_t& rand_eng) noexcept(utils::ndebug_defined()) {
  assert(extr_barriers.size() == mask.size());  // NOLINT
  for (auto i = 0UL; i < extr_barriers.size(); ++i) {
    mask[i] = CTCF::next_state(mask[i] ? CTCF::OCCUPIED : CTCF::NOT_OCCUPIED,
                               extr_barriers[i].prob_occupied_to_occupied(),
                               extr_barriers[i].prob_not_occupied_to_not_occupied(), rand_eng);
  }
}

bp_t Simulation::generate_rev_move(const Chromosome* const chrom, const ExtrusionUnit& unit,
                                   modle::PRNG_t& rand_eng) {
  assert(unit.pos() >= chrom->start_pos());  // NOLINT
  if (this->rev_extrusion_speed_std == 0) {  // When std == 0 always return the avg. extrusion speed
    // (except when unit is close to chrom start pos.)
    return std::min(this->rev_extrusion_speed, unit.pos() - chrom->start_pos());
  }
  // Generate the move distance (and make sure it is not a negative distance)
  // NOTE: on my laptop generating doubles from a normal distribution, rounding them, then
  // casting double to uint is a lot faster (~4x) than drawing uints directly from a Poisson
  // distr.
  return std::clamp(static_cast<bp_t>(std::round(
                        lef_move_generator_t{static_cast<double>(this->rev_extrusion_speed),
                                             this->rev_extrusion_speed_std}(rand_eng))),
                    0UL, unit.pos() - chrom->start_pos());
}

bp_t Simulation::generate_fwd_move(const Chromosome* const chrom, const ExtrusionUnit& unit,
                                   modle::PRNG_t& rand_eng) {
  // See Simulation::generate_rev_move for comments
  assert(unit.pos() < chrom->end_pos());  // NOLINT
  if (this->fwd_extrusion_speed_std == 0) {
    return std::min(this->fwd_extrusion_speed, (chrom->end_pos() - 1) - unit.pos());
  }
  return std::clamp(static_cast<bp_t>(std::round(
                        lef_move_generator_t{static_cast<double>(this->fwd_extrusion_speed),
                                             this->fwd_extrusion_speed_std}(rand_eng))),
                    0UL, (chrom->end_pos() - 1) - unit.pos());
}

void Simulation::generate_moves(const Chromosome* const chrom, absl::Span<const Lef> lefs,
                                absl::Span<const std::size_t> rev_lef_ranks,
                                absl::Span<const std::size_t> fwd_lef_ranks,
                                absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
                                modle::PRNG_t& rand_eng,
                                bool adjust_moves_) noexcept(utils::ndebug_defined()) {
  assert(lefs.size() == fwd_lef_ranks.size());  // NOLINT
  assert(lefs.size() == rev_lef_ranks.size());  // NOLINT
  assert(lefs.size() == fwd_moves.size());      // NOLINT
  assert(lefs.size() == rev_moves.size());      // NOLINT

  // As long as a LEF is bound to DNA, always generate a move
  for (auto i = 0UL; i < lefs.size(); ++i) {
    rev_moves[i] = lefs[i].is_bound() ? generate_rev_move(chrom, lefs[i].rev_unit, rand_eng) : 0UL;
    fwd_moves[i] = lefs[i].is_bound() ? generate_fwd_move(chrom, lefs[i].fwd_unit, rand_eng) : 0UL;
  }

  if (adjust_moves_) {  // Adjust moves of consecutive extr. units to make LEF behavior more
                        // realistic See comments in adjust_moves for more details on what this
                        // entails
    this->adjust_moves(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves);
  }
}

void Simulation::adjust_moves(const Chromosome* chrom, absl::Span<const Lef> lefs,
                              absl::Span<const std::size_t> rev_lef_ranks,
                              absl::Span<const std::size_t> fwd_lef_ranks,
                              absl::Span<bp_t> rev_moves,
                              absl::Span<bp_t> fwd_moves) noexcept(utils::ndebug_defined()) {
  (void)chrom;

  // Loop over pairs of consecutive extr. units.
  // Extr. units moving in rev direction are processed in 3'-5' order, while units moving in fwd
  // direction are processed in 5'-3' direction
  const auto rev_offset = lefs.size() - 1;
  for (auto i = 0UL; i < lefs.size() - 1; ++i) {
    const auto& idx1 = rev_lef_ranks[rev_offset - 1 - i];
    const auto& idx2 = rev_lef_ranks[rev_offset - i];

    if (lefs[idx1].is_bound() && lefs[idx2].is_bound()) {
      assert(lefs[idx1].rev_unit.pos() >= chrom->start_pos() + rev_moves[idx1]);  // NOLINT
      assert(lefs[idx2].rev_unit.pos() >= chrom->start_pos() + rev_moves[idx2]);  // NOLINT

      const auto pos1 = lefs[idx1].rev_unit.pos() - rev_moves[idx1];
      const auto pos2 = lefs[idx2].rev_unit.pos() - rev_moves[idx2];

      // If moving extr. units 1 and 2 by their respective moves would cause unit 1 (which comes
      // first in 3'-5' order) to be surpassed by unit 2 (which comes second in 3'-5' order),
      // increment the move for unit 1 so that after moving both units, unit 1 will be located one
      // bp upstream of unit 2 in 5'-3' direction.
      // This mimics what would probably happen in a real system, where extr. unit 2 would most
      // likely push extr. unit 1, temporarily increasing extr. speed of unit 1.
      if (pos2 < pos1) {
        rev_moves[idx1] += pos1 - pos2;
      }
    }

    const auto& idx3 = fwd_lef_ranks[i];
    const auto& idx4 = fwd_lef_ranks[i + 1];

    // See above for detailed comments. The logic is the same used on rev units (but mirrored!)
    if (lefs[idx3].is_bound() && lefs[idx4].is_bound()) {
      assert(lefs[idx3].fwd_unit.pos() + fwd_moves[idx3] < chrom->end_pos());  // NOLINT
      assert(lefs[idx4].fwd_unit.pos() + fwd_moves[idx4] < chrom->end_pos());  // NOLINT

      const auto pos3 = lefs[idx3].fwd_unit.pos() + fwd_moves[idx3];
      const auto pos4 = lefs[idx4].fwd_unit.pos() + fwd_moves[idx4];

      if (pos3 > pos4) {
        fwd_moves[idx4] += pos3 - pos4;
      }
    }
  }
}

void Simulation::rank_lefs(absl::Span<const Lef> lefs, absl::Span<std::size_t> rev_lef_rank_buff,
                           absl::Span<std::size_t> fwd_lef_rank_buff,
                           bool ranks_are_partially_sorted,
                           bool init_buffers) noexcept(utils::ndebug_defined()) {
  assert(lefs.size() == fwd_lef_rank_buff.size());  // NOLINT
  assert(lefs.size() == rev_lef_rank_buff.size());  // NOLINT

  auto rev_comparator = [&](const auto r1, const auto r2) {
    const auto& p1 = lefs[r1].rev_unit.pos();
    const auto& p2 = lefs[r2].rev_unit.pos();
    if (p1 == p2) MODLE_UNLIKELY {
        // A LEF will have both extr. units bound at the same position if and only if the DNA
        // binding took place in the current iteration.
        // When this is the case, and there are multiple rev. units binding the same bp, then we
        // we want the LEF that bound last to be ranked last in the sequence of tied LEFs
        if (lefs[r1].fwd_unit.pos() == p1) {
          return false;
        }
        if (lefs[r2].fwd_unit.pos() == p2) {
          return true;
        }
        // This return ensures that in case of a tie between rev. units of two LEFs, neither of
        // which was bound the DNA in the current iteration, the existing order is maintained (which
        // means that the LEF that was the first to collide is positioned right next to the extr.
        // barrier, and any other LEFs follow based on their collision iteration)
        return r1 < r2;
      }
    return p1 < p2;
  };

  // See comments for rev_comparator.
  auto fwd_comparator = [&](const auto r1, const auto r2) {
    const auto& p1 = lefs[r1].fwd_unit.pos();
    const auto& p2 = lefs[r2].fwd_unit.pos();
    if (p1 == p2) MODLE_UNLIKELY {
        // Notice that the return value of the following two branches is the opposite of what's
        // found in rev_comparator. This is because in case of a tie we want the LEF that was bound
        // in the current iteration to be ranked lowest
        if (lefs[r1].rev_unit.pos() == p1) {
          return true;
        }
        if (lefs[r2].rev_unit.pos() == p2) {
          return false;
        }
        return r1 < r2;
      }
    return p1 < p2;
  };

  if (init_buffers) MODLE_UNLIKELY {  // Init rank buffers
      std::iota(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), 0);
      std::iota(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), 0);
    }

  if (lefs.size() < 20) MODLE_UNLIKELY {  // TODO Come up with a reasonable threshold
      cppsort::insertion_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), rev_comparator);
      cppsort::insertion_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), fwd_comparator);
    }
  else if (ranks_are_partially_sorted)
    MODLE_LIKELY {
      /*
      fmt::print(
          stderr, "fwd_inv={}\nrev_inv={}\n",
          static_cast<double>(cppsort::probe::inv(fwd_lef_rank_buff.begin(),
      fwd_lef_rank_buff.end(),
                                                  [&](const auto r1, const auto r2) {
                                                    return lefs[r1].fwd_unit.pos() <
                                                           lefs[r2].fwd_unit.pos();
                                                  })) /
              static_cast<double>(fwd_lef_rank_buff.size() * (fwd_lef_rank_buff.size() - 1) / 2),
          static_cast<double>(cppsort::probe::inv(rev_lef_rank_buff.begin(),
      rev_lef_rank_buff.end(),
                                                  [&](const auto r1, const auto r2) {
                                                    return lefs[r1].rev_unit.pos() <
                                                           lefs[r2].rev_unit.pos();
                                                  })) /
              static_cast<double>(rev_lef_rank_buff.size() * (rev_lef_rank_buff.size() - 1) / 2));
      */
      // Drop merge sort is a Rem-adaptive algorithm that it is particularly suitable for our
      // use-case, where after the initial binding, we expect a small fraction of extr. units to be
      // out of place when this function is called (from initial testing, it looks like most of the
      // times we have 2.5-5% inversions, and very rarely > 10%, which is where drop merge sort
      // shines).
      // https://github.com/Morwenn/cpp-sort/blob/develop/docs/Benchmarks.md#inv-adaptive-algorithms
      // https://github.com/Morwenn/cpp-sort/blob/develop/docs/Sorters.md#drop_merge_sorter
      // The algorithm itself is not stable, but the *_comparator should take care of this, making
      // sorting stable in practice
      cppsort::drop_merge_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), rev_comparator);
      cppsort::drop_merge_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), fwd_comparator);
    }
  else
    MODLE_UNLIKELY {
      // Fallback on pattern-defeating quicksort we have no information regarding the level of
      // pre-sortedness of LEFs
      cppsort::pdq_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), rev_comparator);
      cppsort::pdq_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), fwd_comparator);
    }
}

void Simulation::extrude(const Chromosome* chrom, absl::Span<Lef> lefs,
                         absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
                         std::size_t nrev_units_at_5prime,
                         std::size_t nfwd_units_at_3prime) noexcept(utils::ndebug_defined()) {
  assert(lefs.size() == rev_moves.size());      // NOLINT
  assert(lefs.size() == fwd_moves.size());      // NOLINT
  assert(lefs.size() >= nrev_units_at_5prime);  // NOLINT
  assert(lefs.size() >= nfwd_units_at_3prime);  // NOLINT
  (void)chrom;

  auto i1 = nrev_units_at_5prime == 0 ? 0UL : nrev_units_at_5prime - 1;
  const auto i2 = lefs.size() - nfwd_units_at_3prime;
  for (; i1 < i2; ++i1) {
    auto& lef = lefs[i1];
    if (!lef.is_bound()) {  // Do not process inactive LEFs
      continue;
    }
    assert(lef.rev_unit.pos() <= lef.fwd_unit.pos());  // NOLINT

    // Extrude rev unit
    assert(lef.rev_unit.pos() >= chrom->start_pos() + rev_moves[i1]);  // NOLINT
    lef.rev_unit._pos -= rev_moves[i1];  // Advance extr. unit in 3'-5' direction

    // Extrude fwd unit
    assert(lef.fwd_unit.pos() + fwd_moves[i1] <= chrom->end_pos() - 1);  // NOLINT
    lef.fwd_unit._pos += fwd_moves[i1];  // Advance extr. unit in 5'-3' direction

    assert(lef.rev_unit.pos() <= lef.fwd_unit.pos());  // NOLINT
  }
}

template <typename I, typename MaskT>
std::pair<std::size_t, std::size_t> Simulation::process_collisions(
    const Chromosome* chrom, absl::Span<const Lef> lefs,
    absl::Span<const ExtrusionBarrier> barriers, const MaskT& barrier_mask,
    absl::Span<const std::size_t> rev_lef_ranks, absl::Span<const std::size_t> fwd_lef_ranks,
    absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
    absl::Span<I> fwd_collisions, PRNG_t& rand_eng) noexcept(utils::ndebug_defined()) {
  // TODO Update this comment!
  // NOTE: The call order here is important!
  //       TLDR: Always call process_lef_bar* before process_lef_lef*
  //       process_lef_lef_collision assumes that the collision masks contain already LEF-BAR
  //       collisions.

  // process_*_collisions detect LEF-BAR/LEF collisions and update the moves so that after
  // extude() is called, stalled LEFs are located at their respective collision site
  const auto& [nrev_units_at_5prime, nfwd_units_at_3prime] =
      Simulation::detect_collisions_at_chrom_boundaries(chrom, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                        rev_moves, fwd_moves, rev_collisions,
                                                        fwd_collisions);

  this->detect_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                                  barriers, barrier_mask, rev_collisions, fwd_collisions, rand_eng,
                                  nrev_units_at_5prime, nfwd_units_at_3prime);

  this->detect_primary_lef_lef_collisions(chrom, lefs, barriers, rev_lef_ranks, fwd_lef_ranks,
                                          rev_moves, fwd_moves, rev_collisions, fwd_collisions,
                                          rand_eng, nrev_units_at_5prime, nfwd_units_at_3prime);
  this->adjust_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves, rev_collisions,
                                            fwd_collisions);

  this->adjust_moves_for_primary_lef_lef_collisions(lefs, barriers, rev_lef_ranks, fwd_lef_ranks,
                                                    rev_moves, fwd_moves, rev_collisions,
                                                    fwd_collisions);
  this->detect_secondary_lef_lef_collisions(
      chrom, lefs, barriers.size(), rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
      rev_collisions, fwd_collisions, rand_eng, nrev_units_at_5prime, nfwd_units_at_3prime);

  this->adjust_moves_for_secondary_lef_lef_collisions(lefs, barriers.size(), rev_lef_ranks,
                                                      fwd_lef_ranks, rev_moves, fwd_moves,
                                                      rev_collisions, fwd_collisions);

  return std::make_pair(nrev_units_at_5prime, nfwd_units_at_3prime);
}

template <typename I>
std::pair<std::size_t, std::size_t> Simulation::detect_collisions_at_chrom_boundaries(
    const Chromosome* chrom, absl::Span<const Lef> lefs,
    absl::Span<const std::size_t> rev_lef_ranks, absl::Span<const std::size_t> fwd_lef_ranks,
    absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
    absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");
  {
    assert(lefs.size() == fwd_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == rev_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == fwd_moves.size());                           // NOLINT
    assert(lefs.size() == rev_moves.size());                           // NOLINT
    assert(lefs.size() == fwd_collisions.size());                      // NOLINT
    assert(lefs.size() == rev_collisions.size());                      // NOLINT
    assert(std::is_sorted(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),  // NOLINT
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                          }));
    assert(std::is_sorted(rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                          }));
    assert(std::all_of(rev_collisions.begin(), rev_collisions.end(),  // NOLINT
                       [](const auto c) { return c == NO_COLLISION; }));
    assert(std::all_of(fwd_collisions.begin(), fwd_collisions.end(),  // NOLINT
                       [](const auto c) { return c == NO_COLLISION; }));
  }
  // Detect if the first rev unit or last fwd unit are about to fall off chrom. boundaries
  // Also detect extr. units that are already at chrom boundaries

  auto nrev_units_at_5prime = 0UL;
  auto nfwd_units_at_3prime = 0UL;

  assert(lefs[fwd_lef_ranks[0]].fwd_unit.pos() != std::numeric_limits<bp_t>::max());  // NOLINT
  const auto& first_active_fwd_unit = lefs[fwd_lef_ranks[0]].fwd_unit;
  const auto& last_active_rev_unit =
      lefs[*std::find_if(rev_lef_ranks.rbegin(), rev_lef_ranks.rend(), [&](const auto i) {
        return lefs[i].is_bound();
      })].rev_unit;

  for (auto i = 0UL; i < lefs.size(); ++i) {
    const auto& rev_idx = rev_lef_ranks[i];
    const auto& rev_unit = lefs[rev_idx].rev_unit;
    auto& rev_move = rev_moves[rev_idx];
    assert(lefs[rev_idx].is_bound());                         // NOLINT
    assert(chrom->start_pos() + rev_move <= rev_unit.pos());  // NOLINT
    if (rev_unit.pos() == chrom->start_pos()) {
      assert(rev_moves[rev_idx] == 0);  // NOLINT
      ++nrev_units_at_5prime;
      rev_collisions[rev_idx] = REACHED_CHROM_BOUNDARY;
    } else if (rev_unit.pos() > first_active_fwd_unit.pos()) {
      break;
    } else if (rev_unit.pos() - rev_move == chrom->start_pos()) {
      rev_collisions[rev_idx] = REACHED_CHROM_BOUNDARY;
      ++nrev_units_at_5prime;
      break;
    }
  }

  for (auto i = lefs.size() - 1; i > 0; --i) {
    const auto& fwd_idx = fwd_lef_ranks[i];
    const auto& fwd_unit = lefs[fwd_idx].fwd_unit;
    const auto& fwd_move = fwd_moves[fwd_idx];
    if (!lefs[fwd_idx].is_bound()) {
      ++nfwd_units_at_3prime;  // Inactive units are technically not at the 3'-end, but we count
      continue;                // them anyway so that we can shrink the spans on LEF-related
                               // buffers to avoid doing some work in later steps
    }
    assert(fwd_unit.pos() + fwd_move < chrom->end_pos());  // NOLINT
    if (fwd_unit.pos() == chrom->end_pos() - 1) {
      assert(fwd_moves[fwd_idx] == 0);  // NOLINT
      ++nfwd_units_at_3prime;
      fwd_collisions[fwd_idx] = REACHED_CHROM_BOUNDARY;
    } else if (fwd_unit.pos() < last_active_rev_unit.pos()) {
      break;
    } else if (fwd_unit.pos() + fwd_move == chrom->end_pos() - 1) {
      fwd_collisions[fwd_idx] = REACHED_CHROM_BOUNDARY;
      ++nfwd_units_at_3prime;
      break;
    }
  }

  return std::make_pair(nrev_units_at_5prime, nfwd_units_at_3prime);
}

template <typename I>
void Simulation::detect_lef_bar_collisions(
    absl::Span<const Lef> lefs, absl::Span<const std::size_t> rev_lef_ranks,
    absl::Span<const std::size_t> fwd_lef_ranks, absl::Span<const bp_t> rev_moves,
    absl::Span<const bp_t> fwd_moves, absl::Span<const ExtrusionBarrier> extr_barriers,
    const boost::dynamic_bitset<>& barrier_mask, absl::Span<I> rev_collisions,
    absl::Span<I> fwd_collisions, modle::PRNG_t& rand_eng, std::size_t nrev_units_at_5prime,
    std::size_t nfwd_units_at_3prime) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");
  {
    assert(lefs.size() == fwd_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == rev_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == fwd_moves.size());                           // NOLINT
    assert(lefs.size() == rev_moves.size());                           // NOLINT
    assert(lefs.size() == fwd_collisions.size());                      // NOLINT
    assert(lefs.size() == rev_collisions.size());                      // NOLINT
    assert(barrier_mask.size() == extr_barriers.size());               // NOLINT
    assert(lefs.size() >= nrev_units_at_5prime);                       // NOLINT
    assert(lefs.size() >= nfwd_units_at_3prime);                       // NOLINT
    assert(std::is_sorted(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),  // NOLINT
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                          }));
    assert(std::is_sorted(rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                          }));
  }
  // Loop over LEFs, using a procedure similar to merge in mergesort.
  // The idea here is that if we have a way to visit extr. units in 5'-3' order, then detecting
  // LEF-BAR collisions boils down to the following:
  //  - Loop over all extrusion barriers:
  //    - Find the first rev. unit located downstream of the current extr. barrier
  //    - Find the last fwd. unit located upstream of the current extr. barrier
  //    - If the distance between the extr. barrier and rev/fwd extr. units is less than the
  //      distance that will be covered by the rev/fwd extr. unit in the current iteration, then
  //      there may be a collision.
  //      To determine whether or not there will be a collision, run a bernoulli trial with a
  //      probability of success equal to the probability of block of the current barrier
  //
  // When a collision is detected the appropriate entry in rev/fwd_collision buffer is set to the
  // index of the extrusion barrier that caused the collision.
  // Furthermore the appropriate entry in the move vector is updated such that after calling
  // Simulation::extrude, the extr. unit that is being blocked will be located 1pb up/downstream
  // of the extr. barrier that caused the collision.

  // Init indices and position with the rev and fwd extr. units that are the closest to the 5'-end
  auto j1 = nrev_units_at_5prime == 0UL ? 0 : nrev_units_at_5prime - 1;
  std::size_t j2 = 0;
  const auto j2_end = lefs.size() - nfwd_units_at_3prime;

  auto rev_idx = rev_lef_ranks[j1];
  auto fwd_idx = fwd_lef_ranks[j2];
  auto rev_unit_pos = lefs[rev_idx].rev_unit.pos();
  auto fwd_unit_pos = lefs[fwd_idx].fwd_unit.pos();

  // Loop over extr. barriers and find the first, possibly colliding extr. unit
  for (auto i = 0UL; i < extr_barriers.size(); ++i) {
    if (barrier_mask[i] == CTCF::NOT_OCCUPIED) {  // Extrusion barriers that are not occupied are
      continue;                                   // transparent to extr. units
    }

    const auto& barrier = extr_barriers[i];

    if (j1 < lefs.size()) {  // Process rev unit
      // Probability of block is set based on the extr. barrier blocking direction
      const auto& pblock = barrier.blocking_direction_major() == dna::rev
                               ? this->lef_hard_collision_pblock
                               : this->lef_soft_collision_pblock;

      // Look for the first rev extr. unit that comes after the current barrier
      while (rev_unit_pos <= barrier.pos()) {
        if (++j1 >= lefs.size()) {  // All rev units have been processed
          goto process_fwd_unit;    // Move to the next section
        }

        // Update idx and position with those corresponding to the rev extr. unit that comes next
        rev_idx = rev_lef_ranks[j1];  // in 5'-3' order
        rev_unit_pos = lefs[rev_idx].rev_unit.pos();
      }

      if (lefs[rev_idx].is_bound()) {
        // We have a LEF-BAR collision event if the distance between the rev. unit and the extr.
        // barrier is less or equal than the distance that the rev extr. unit is set to move in
        // the current iteration. If pblock != 1, then we also require a successful bernoulli
        // trial before calling a collision
        const auto delta = rev_unit_pos - barrier.pos();
        if (delta > 0 && delta <= rev_moves[rev_idx] &&
            (pblock == 1.0 || std::bernoulli_distribution{pblock}(rand_eng))) {
          // Collision detected. Assign barrier idx to the respective entry in the collision mask
          rev_collisions[rev_idx] = i;
          // Move LEF close to the extr. barrier (i.e 1bp upstream of the extr. barrier)
          // rev_move = delta > 1 ? delta - 1 : 0;
        }
      }
    }

    // Look in the previous section for detailed comments
  process_fwd_unit:
    if (j2 < j2_end) {
      const auto& pblock = barrier.blocking_direction_major() == dna::fwd
                               ? this->lef_hard_collision_pblock
                               : this->lef_soft_collision_pblock;
      // Look for the next fwd unit that comes strictly before the current extr. barrier
      while (fwd_unit_pos < barrier.pos()) {
        if (++j2 >= j2_end) {
          goto end_of_loop;
        }

        fwd_idx = fwd_lef_ranks[j2];
        fwd_unit_pos = lefs[fwd_idx].fwd_unit.pos();
      }

      // Decrement j2 by one (if it is legal to do so), so that j2 corresponds to the index of the
      // fwd extr. unit that is located as close as possible to the extr. barrier that is being
      // processed
      fwd_idx = fwd_lef_ranks[j2 > 0 ? --j2 : 0];
      fwd_unit_pos = lefs[fwd_idx].fwd_unit.pos();

      if (lefs[fwd_idx].is_bound()) {
        const auto delta = barrier.pos() - fwd_unit_pos;
        if (delta > 0 && delta <= fwd_moves[fwd_idx] &&
            (pblock == 1.0 || std::bernoulli_distribution{pblock}(rand_eng))) {
          fwd_collisions[fwd_idx] = i;
          // fwd_move = delta > 1 ? delta - 1 : 0;
        }
      }
    }

  end_of_loop:
    // Return immediately if all extr. units have been processed (regardless of whether there are
    // still extr. barriers to be processed)
    if (j1 == lefs.size() && j2 == j2_end) {
      return;
    }
  }
}

template <typename I>  // TODO: make move spans const
void Simulation::detect_lef_lef_collisions(
    const Chromosome* chrom, absl::Span<const Lef> lefs,
    absl::Span<const ExtrusionBarrier> barriers, absl::Span<const std::size_t> rev_lef_ranks,
    absl::Span<const std::size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
    absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions,
    modle::PRNG_t& rand_eng, std::size_t nrev_units_at_5prime,
    std::size_t nfwd_units_at_3prime) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");

  // Process collisions between LEFs moving in opposite directions
  this->detect_primary_lef_lef_collisions(chrom, lefs, barriers, rev_lef_ranks, fwd_lef_ranks,
                                          rev_moves, fwd_moves, rev_collisions, fwd_collisions,
                                          rand_eng, nrev_units_at_5prime, nfwd_units_at_3prime);

  // Process collisions between LEFs moving in the same direction
  this->detect_secondary_lef_lef_collisions(
      chrom, lefs, barriers.size(), rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
      rev_collisions, fwd_collisions, rand_eng, nrev_units_at_5prime, nfwd_units_at_3prime);
}

template <typename I>  // TODO: make move spans const
void Simulation::detect_primary_lef_lef_collisions(
    const Chromosome* chrom, absl::Span<const Lef> lefs,
    absl::Span<const ExtrusionBarrier> barriers, absl::Span<const std::size_t> rev_lef_ranks,
    absl::Span<const std::size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
    absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions,
    PRNG_t& rand_eng, std::size_t nrev_units_at_5prime,
    std::size_t nfwd_units_at_3prime) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");
  {
    assert(lefs.size() == fwd_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == rev_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == fwd_moves.size());                           // NOLINT
    assert(lefs.size() == rev_moves.size());                           // NOLINT
    assert(lefs.size() == fwd_collisions.size());                      // NOLINT
    assert(lefs.size() == rev_collisions.size());                      // NOLINT
    assert(lefs.size() >= nrev_units_at_5prime);                       // NOLINT
    assert(lefs.size() >= nfwd_units_at_3prime);                       // NOLINT
    assert(std::is_sorted(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),  // NOLINT
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                          }));
    assert(std::is_sorted(rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                          }));
  }
  (void)chrom;  // TODO remove me
  // Loop over LEFs, using a procedure similar to merge in mergesort
  // The idea here is that if we have a way to visit extr. units in 5'-3' order, then detecting
  // LEF-LEF collisions boils down to:
  //  - Starting from the first fwd extr. unit in 5'-3' order:
  //    - Look for the first rev extr. unit that is located downstream of the fwd unit that is being
  //    processed
  //    - If the distance between the pair of extr. units identified in the previous step is below a
  //      certain threshold, then we may have a LEF-LEF collision.
  //      If the probability of unit bypass is 0, then we always predict a LEF-LEF collision. If
  //      this probability is larger than 0, then we predict a LEF-LEF collision based on the
  //      outcome of a bernoulli trial with probability of success equal to 1 - the prob. of bypass
  //    - If a LEF-LEF collision caused by two extr. unit moving in opposite direction is detected,
  //      encode the index of the extr. unit that caused the collision in the appropriate collision
  //      mask. This kind of collisions are encoded as offset + i, where offset = nbarriers and i =
  //      the index of the unit that is colliding.

  if (nrev_units_at_5prime == lefs.size() || nfwd_units_at_3prime == lefs.size()) MODLE_UNLIKELY {
      return;
    }

  //    Initialize indexes so that we skip over rev units at the 5' and fwd units at the 3' (if any)
  auto i1 = 0UL, j1 = nrev_units_at_5prime;
  const auto i2 = lefs.size() - nfwd_units_at_3prime, j2 = lefs.size();

  while (true) {
    auto rev_idx = rev_lef_ranks[j1];             // index of the jth rev unit in 5'-3' order
    auto rev_pos = lefs[rev_idx].rev_unit.pos();  // pos of the jth unit
    auto fwd_idx = fwd_lef_ranks[i1];             // index of the ith fwd unit in 5'-3' order
    auto fwd_pos = lefs[fwd_idx].fwd_unit.pos();  // pos of the ith unit

    // Find the first rev unit that comes right after the ith fwd unit
    while (rev_pos <= fwd_pos) {
      if (++j1 == j2) MODLE_UNLIKELY {
          return;  // all rev units have been processed
        }
      rev_idx = rev_lef_ranks[j1];
      rev_pos = lefs[rev_idx].rev_unit.pos();
    }

    // Find the last fwd unit that comes right before the jth rev unit
    // This is necessary in order to handle ties in fwd units ranking, as well as the case where
    // there are many fwd units between the pos of the jth-1 and jth rev extr. units
    while (fwd_pos < rev_pos) {
      if (++i1 == i2) MODLE_UNLIKELY {
          return;  // all fwd units have been processed
        }
      fwd_idx = fwd_lef_ranks[i1];             // index of the ith fwd unit in 5'-3' order
      fwd_pos = lefs[fwd_idx].fwd_unit.pos();  // pos of the ith fwd unit
    }

    // The previous while loop finds the first fwd unit that comes at or after the rev unit that
    // is being processed, but we are actually interested in the fwd unit before this one.
    // This is why we take the unit at index i1-1
    fwd_idx = fwd_lef_ranks[i1 - 1];
    fwd_pos = lefs[fwd_idx].fwd_unit.pos();

    // We have a LEF-LEF collision event if the distance between the two extr. units is less than
    // the sum of their moves.
    // If probability_of_extrusion_unit_bypass != 0, then we also require a successful bernoulli
    // trial before calling a collision
    if (const auto delta = rev_pos - fwd_pos;
        delta > 0 && delta < rev_moves[rev_idx] + fwd_moves[fwd_idx] &&
        (this->probability_of_extrusion_unit_bypass == 0 ||
         std::bernoulli_distribution{1.0 - this->probability_of_extrusion_unit_bypass}(rand_eng))) {
      // Declare few aliases to reduce code verbosity later on
      auto& rev_move = rev_moves[rev_idx];
      auto& fwd_move = fwd_moves[fwd_idx];
      auto& cause_of_collision_rev = rev_collisions[rev_idx];
      auto& cause_of_collision_fwd = fwd_collisions[fwd_idx];

      // Note that the rev collision pos will always be 1bp upstream of the fwd collision pos
      auto [collision_pos_rev, collision_pos_fwd] = compute_lef_lef_collision_pos(
          lefs[rev_idx].rev_unit, lefs[fwd_idx].fwd_unit, rev_move, fwd_move);

      if (cause_of_collision_rev == NO_COLLISION && cause_of_collision_fwd == NO_COLLISION) {
        // In the simplest case both units are free to move.
        // Thus we update their respective moves so that after calling Simulate::extrude, the two
        // units will be locate at their respective collision sites
        cause_of_collision_rev = barriers.size() + fwd_idx;
        cause_of_collision_fwd = barriers.size() + rev_idx;

      } else if (cause_of_collision_rev != NO_COLLISION && cause_of_collision_fwd == NO_COLLISION) {
        // In this case only the fwd unit is free to move.
        // This case is a bit more complicated than the first one, because we have to handle the
        // edge case where we have mistakenly predicted a LEF-BAR collision. This usually happens
        // when two extr. units are moving in opposite directions, and the unit extruding in rev
        // direction is located downstream of the unit extruding in fwd direction. If located
        // between the two units there's an extr. barrier, and one of the two units manages to
        // bypass the extr. barrier, then there's a chance that the LEF-LEF collision will take
        // place before the LEF-BAR collision (whether this happens or not depends on the distance
        // between the extr. units and the extr. barriers, as well as their respective extr. speed

        assert(cause_of_collision_rev < barriers.size() /* is LEF-BAR collision*/);  // NOLINT
        assert(collision_pos_rev != 0 && collision_pos_fwd != 0);                    // NOLINT
        const auto& barrier_pos = barriers[cause_of_collision_rev].pos();
        if (collision_pos_fwd > barrier_pos) {
          // Detected the mis-prediction mentioned above: make the LEF-BAR collision a LEF-LEF
          // collision
          cause_of_collision_rev = barriers.size() + fwd_idx;
          cause_of_collision_fwd = barriers.size() + rev_idx;
        } else {
          // fwd extr unit is being blocked by a rev unit that is itself being blocked by an extr.
          // barrier
          cause_of_collision_fwd = barriers.size() + rev_idx;
        }
        // This branch follows the same logic as the previous one. In this case the rev unit is free
        // to move, while the fwd unit has been predicted to be stalled by an extr. barrier
      } else if (cause_of_collision_rev == NO_COLLISION && cause_of_collision_fwd != NO_COLLISION) {
        assert(cause_of_collision_fwd < barriers.size() /* is LEF-BAR collision*/);  // NOLINT
        assert(collision_pos_rev != 0 && collision_pos_fwd != 0);                    // NOLINT
        const auto& barrier_pos = barriers[cause_of_collision_fwd].pos();
        if (collision_pos_rev < barrier_pos) {
          cause_of_collision_rev = barriers.size() + fwd_idx;
          cause_of_collision_fwd = barriers.size() + rev_idx;
        } else {
          cause_of_collision_rev = barriers.size() + fwd_idx;
        }
      }
    }
  }
}

template <typename I>  // TODO Make move spans const
void Simulation::detect_secondary_lef_lef_collisions(
    const Chromosome* chrom, absl::Span<const Lef> lefs, std::size_t nbarriers,
    absl::Span<const std::size_t> rev_lef_ranks, absl::Span<const std::size_t> fwd_lef_ranks,
    absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
    absl::Span<I> fwd_collisions, PRNG_t& rand_eng, std::size_t nrev_units_at_5prime,
    std::size_t nfwd_units_at_3prime) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");
  {
    assert(lefs.size() == fwd_lef_ranks.size());   // NOLINT
    assert(lefs.size() == rev_lef_ranks.size());   // NOLINT
    assert(lefs.size() == fwd_moves.size());       // NOLINT
    assert(lefs.size() == rev_moves.size());       // NOLINT
    assert(lefs.size() == fwd_collisions.size());  // NOLINT
    assert(lefs.size() == rev_collisions.size());  // NOLINT
    assert(lefs.size() >= nrev_units_at_5prime);   // NOLINT
    assert(lefs.size() >= nfwd_units_at_3prime);   // NOLINT
    (void)chrom;
  }

  // Loop over pairs of consecutive fwd units.
  // Throughout the comments for this function we will use the following notation:
  //  - U1 = unit at index i - 1 (i.e. the unit closer to the 5'-end);
  //  - U2 = unit at index i (i.e. the unit closer to the 3'-end);
  //
  // 0. Keep advancing the index i until we find the first U2 unit that is stalled. Then:
  //     1. Check if calling Simulation::extrude would cause U1 to end up at the same position or
  //        downstream of U2.
  //        If this is the case, then we may have a LEF-LEF collision.
  //          2. If the probability of bypass is 0, then we certainly have a collision, otherwise
  //             whether a collision will occur or not is decided by a bernoulli trial with
  //             probability of success equal to 1 - prob. of bypass.
  //          3. If the previous step determines that a LEF-LEF collision is about to take place,
  //             update the appropriate collision mask with the encoded index of the extr. unit that
  //             is causing the collision. Index i is encoded as nbarriers + nlefs + i.
  //          4. Keep backtracking repeating steps 1-4 until step 1 fails
  // 5. Continue iterating from where we left at step 0

  const auto offset = nbarriers + lefs.size();

  // TBC after fixing some bugs
  for (auto i = lefs.size() - nfwd_units_at_3prime - 1; i > 1; --i) {
    const auto& fwd_idx2 = fwd_lef_ranks[i];  // index of the ith fwd unit in 5'-3' order
    const auto& fwd_pos2 = lefs[fwd_idx2].fwd_unit.pos();  // pos of the ith unit
    if (fwd_collisions[fwd_idx2] == NO_COLLISION) {
      continue;
    }

    const auto& fwd_idx1 = fwd_lef_ranks[i - 1];  // index of the ith-1 fwd unit in 5'-3' order
    const auto& fwd_pos1 = lefs[fwd_idx1].fwd_unit.pos();  // pos of the ith-1 unit
    if (fwd_collisions[fwd_idx1] != NO_COLLISION) {
      continue;
    }

    auto& move1 = fwd_moves[fwd_idx1];
    const auto& move2 = fwd_moves[fwd_idx2];

    assert(fwd_pos2 >= fwd_pos1);                 // NOLINT
    assert(fwd_pos1 + move1 < chrom->end_pos());  // NOLINT
    assert(fwd_pos2 + move2 < chrom->end_pos());  // NOLINT

    if (fwd_pos1 + move1 >= fwd_pos2 + move2 &&
        (this->probability_of_extrusion_unit_bypass == 0 ||
         std::bernoulli_distribution{1.0 - this->probability_of_extrusion_unit_bypass}(rand_eng))) {
      fwd_collisions[fwd_idx1] = offset + fwd_idx2;
    }
  }

  for (auto i = std::max(1UL, nrev_units_at_5prime); i < lefs.size(); ++i) {
    const auto& rev_idx1 = rev_lef_ranks[i - 1];  // index of the ith rev unit in 5'-3' order
    const auto& rev_pos1 = lefs[rev_idx1].rev_unit.pos();  // pos of the ith unit
    if (rev_collisions[rev_idx1] == NO_COLLISION) {
      continue;
    }

    const auto& rev_idx2 = rev_lef_ranks[i];  // index of the ith-1 rev unit in 5'-3' order
    const auto& rev_pos2 = lefs[rev_idx2].rev_unit.pos();  // pos of the ith-1 unit
    if (rev_collisions[rev_idx2] != NO_COLLISION) {
      continue;
    }

    const auto& move1 = rev_moves[rev_idx1];
    auto& move2 = rev_moves[rev_idx2];

    assert(rev_pos2 >= rev_pos1);                    // NOLINT
    assert(rev_pos1 >= move1);                       // NOLINT
    assert(rev_pos2 >= move2);                       // NOLINT
    assert(rev_pos1 - move1 >= chrom->start_pos());  // NOLINT
    assert(rev_pos2 - move2 >= chrom->start_pos());  // NOLINT

    if (rev_pos2 - move2 <= rev_pos1 - move1 &&
        (this->probability_of_extrusion_unit_bypass == 0 ||
         std::bernoulli_distribution{1.0 - this->probability_of_extrusion_unit_bypass}(rand_eng))) {
      rev_collisions[rev_idx2] = offset + rev_idx1;
    }
  }
}

template <typename I>
void Simulation::adjust_moves_for_lef_bar_collisions(
    absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
    absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
    absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");

  const auto upper_bound = barriers.size();

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (rev_collisions[i] < upper_bound) MODLE_UNLIKELY {
        const auto& barrier_idx = rev_collisions[i];
        const auto& barrier = barriers[barrier_idx];
        assert(lefs[i].rev_unit.pos() > barrier.pos());  // NOLINT

        const auto delta = lefs[i].rev_unit.pos() - barrier.pos();
        rev_moves[i] = delta > 1UL ? delta - 1 : 0UL;
      }

    if (fwd_collisions[i] < upper_bound) MODLE_UNLIKELY {
        const auto& barrier_idx = fwd_collisions[i];
        const auto& barrier = barriers[barrier_idx];
        assert(lefs[i].fwd_unit.pos() < barrier.pos());  // NOLINT

        const auto delta = barrier.pos() - lefs[i].fwd_unit.pos();
        fwd_moves[i] = delta > 1UL ? delta - 1 : 0UL;
      }
  }
}

template <typename I>
inline void Simulation::adjust_moves_for_primary_lef_lef_collisions(
    absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
    absl::Span<const std::size_t> rev_ranks, absl::Span<const std::size_t> fwd_ranks,
    absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
    absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");

  (void)fwd_ranks;  // TODO Remove me
  const auto lower_bound = barriers.size();
  const auto upper_bound = lower_bound + lefs.size();

  auto is_lef_lef_primary_collision = [&](const auto i) constexpr {
    return i >= lower_bound && i < upper_bound;
  };

  auto is_lef_bar_collision = [&](const auto i) constexpr { return i < lower_bound; };

  for (auto rev_idx : rev_ranks) {
    if (auto rev_collision = rev_collisions[rev_idx]; is_lef_lef_primary_collision(rev_collision))
      MODLE_UNLIKELY {
        rev_collision -= lower_bound;
        const auto& fwd_idx = rev_collision;
        if (auto fwd_collision = fwd_collisions[fwd_idx];
            is_lef_lef_primary_collision(fwd_collision)) {
          const auto& rev_unit = lefs[rev_idx].rev_unit;
          const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

          auto& rev_move = rev_moves[rev_idx];
          auto& fwd_move = fwd_moves[fwd_idx];

          const auto [p1, p2] =
              compute_lef_lef_collision_pos(rev_unit, fwd_unit, rev_move, fwd_move);
          assert(rev_unit.pos() >= p1);  // NOLINT
          assert(fwd_unit.pos() <= p2);  // NOLINT

          rev_move = rev_unit.pos() - p1;
          fwd_move = p2 - fwd_unit.pos();
        } else if (is_lef_bar_collision(fwd_collision)) {
          const auto& rev_unit = lefs[rev_idx].rev_unit;
          const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

          auto& rev_move = rev_moves[rev_idx];
          const auto& fwd_move = fwd_moves[fwd_idx];

          assert(rev_unit.pos() >= fwd_unit.pos() + fwd_move);  // NOLINT
          rev_move = rev_unit.pos() - (fwd_unit.pos() + fwd_move);
        }
      }
  }

  for (auto fwd_idx : fwd_ranks) {
    if (auto fwd_collision = fwd_collisions[fwd_idx]; is_lef_lef_primary_collision(fwd_collision))
      MODLE_UNLIKELY {
        fwd_collision -= lower_bound;
        const auto& rev_idx = fwd_collision;
        if (const auto& rev_collision = rev_collisions[rev_idx];
            is_lef_bar_collision(rev_collision)) {
          const auto& rev_unit = lefs[rev_idx].rev_unit;
          const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

          const auto& rev_move = rev_moves[rev_idx];
          auto& fwd_move = fwd_moves[fwd_idx];

          assert(rev_unit.pos() >= fwd_unit.pos() + rev_move);  // NOLINT
          fwd_move = (rev_unit.pos() - rev_move) - fwd_unit.pos();
        }
      }
  }
}

template <typename I>
inline void Simulation::adjust_moves_for_secondary_lef_lef_collisions(
    absl::Span<const Lef> lefs, std::size_t nbarriers, absl::Span<const std::size_t> rev_ranks,
    absl::Span<const std::size_t> fwd_ranks, absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
    absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");

  const auto lower_bound = nbarriers + lefs.size();
  const auto upper_bound = lower_bound + lefs.size();

  for (const auto rev_idx2 : rev_ranks) {  // TODO skip the first rank
    if (auto rev_idx1 = rev_collisions[rev_idx2]; rev_idx1 >= lower_bound && rev_idx1 < upper_bound)
      MODLE_UNLIKELY {
        rev_idx1 -= lower_bound;
        const auto& rev_unit1 = lefs[rev_idx1].rev_unit;
        const auto& rev_unit2 = lefs[rev_idx2].rev_unit;

        assert(rev_unit1.pos() >= rev_moves[rev_idx1]);  // NOLINT
        assert(rev_unit2.pos() >= rev_moves[rev_idx2]);  // NOLINT
        assert(rev_unit1.pos() - rev_moves[rev_idx1] >=
               rev_unit2.pos() - rev_moves[rev_idx2]);  // NOLINT

        const auto move = rev_unit2.pos() - (rev_unit1.pos() - rev_moves[rev_idx1]);
        rev_moves[rev_idx2] = move > 0UL ? move - 1 : 0UL;
      }
  }

  for (const auto fwd_idx1 : boost::adaptors::reverse(fwd_ranks)) {  // TODO skip the last rank
    if (auto fwd_idx2 = fwd_collisions[fwd_idx1]; fwd_idx2 >= lower_bound && fwd_idx2 < upper_bound)
      MODLE_UNLIKELY {
        fwd_idx2 -= lower_bound;
        const auto& fwd_unit1 = lefs[fwd_idx1].fwd_unit;
        const auto& fwd_unit2 = lefs[fwd_idx2].fwd_unit;

        assert(fwd_unit1.pos() <= fwd_unit2.pos());  // NOLINT
        assert(fwd_unit1.pos() + fwd_moves[fwd_idx1] >=
               fwd_unit2.pos() + fwd_moves[fwd_idx2]);  // NOLINT

        const auto delta = (fwd_unit2.pos() + fwd_moves[fwd_idx2]) - fwd_unit1.pos();
        fwd_moves[fwd_idx1] = delta > 0UL ? delta - 1 : 0UL;
      }
  }
}

template <typename I>
void Simulation::adjust_moves_based_on_collisions(
    absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
    absl::Span<const std::size_t> rev_ranks, absl::Span<const std::size_t> fwd_ranks,
    absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
    absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined()) {
  Simulation::adjust_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                  rev_collisions, fwd_collisions);
  Simulation::adjust_moves_for_primary_lef_lef_collisions(
      lefs, barriers, rev_ranks, fwd_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);
  Simulation::adjust_moves_for_secondary_lef_lef_collisions(lefs, barriers.size(), rev_ranks,
                                                            fwd_ranks, rev_moves, fwd_moves,
                                                            rev_collisions, fwd_collisions);
}

std::pair<bp_t, bp_t> Simulation::compute_lef_lef_collision_pos(const ExtrusionUnit& rev_unit,
                                                                const ExtrusionUnit& fwd_unit,
                                                                bp_t rev_move, bp_t fwd_move) {
  const auto& rev_speed = rev_move;
  const auto& fwd_speed = fwd_move;
  const auto& rev_pos = rev_unit.pos();
  const auto& fwd_pos = fwd_unit.pos();

  const auto relative_speed = rev_speed + fwd_speed;
  const auto time_to_collision =
      static_cast<double>(rev_pos - fwd_pos) / static_cast<double>(relative_speed);

  const auto collision_pos =
      fwd_pos + static_cast<bp_t>(std::round((static_cast<double>(fwd_speed) * time_to_collision)));
  assert(collision_pos <= rev_pos);  // NOLINT
#ifndef NDEBUG
  const auto collision_pos_ =
      static_cast<double>(rev_pos) - static_cast<double>(rev_speed) * time_to_collision;
  assert(abs(static_cast<double>(collision_pos) - collision_pos_) < 1.0);  // NOLINT
#endif
  if (collision_pos == fwd_pos) {
    assert(collision_pos >= fwd_pos);      // NOLINT
    assert(collision_pos + 1 <= rev_pos);  // NOLINT
    return std::make_pair(collision_pos + 1, collision_pos);
  }
  assert(collision_pos > 0);             // NOLINT
  assert(collision_pos - 1 >= fwd_pos);  // NOLINT
  return std::make_pair(collision_pos, collision_pos - 1);
}

std::size_t Simulation::register_contacts(
    Chromosome* chrom, absl::Span<const Lef> lefs,
    absl::Span<const std::size_t> selected_lef_idx) noexcept(utils::ndebug_defined()) {
  // Register contacts for the selected LEFs (excluding LEFs that have one of their units at the
  // beginning/end of a chromosome)
  std::size_t new_contacts = 0;
  for (const auto i : selected_lef_idx) {
    assert(i < lefs.size());  // NOLINT
    const auto& lef = lefs[i];
    if (lef.is_bound() && lef.rev_unit.pos() > chrom->start_pos() &&
        lef.rev_unit.pos() < chrom->end_pos() && lef.fwd_unit.pos() > chrom->start_pos() &&
        lef.fwd_unit.pos() < chrom->end_pos() - 1)
      MODLE_LIKELY {
        chrom->increment_contacts(lef.rev_unit.pos(), lef.fwd_unit.pos(), this->bin_size);
        ++new_contacts;
      }
  }
  return new_contacts;
}

template <typename MaskT>
void Simulation::select_lefs_to_bind(absl::Span<const Lef> lefs,
                                     MaskT& mask) noexcept(utils::ndebug_defined()) {
  static_assert(
      std::is_integral_v<
          std::decay_t<decltype(std::declval<MaskT&>().operator[](std::declval<std::size_t>()))>>,
      "mask should be a vector of integral numbers or a boost::dynamic_bitset.");
  assert(lefs.size() == mask.size());  // NOLINT
  std::transform(lefs.begin(), lefs.end(), mask.begin(),
                 [](const auto& lef) { return !lef.is_bound(); });
}

void Simulation::generate_lef_unloader_affinities(
    absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
    absl::Span<const collision_t> rev_collisions, absl::Span<const collision_t> fwd_collisions,
    absl::Span<double> lef_unloader_affinity) noexcept(utils::ndebug_defined()) {
  assert(lefs.size() == rev_collisions.size());         // NOLINT
  assert(lefs.size() == fwd_collisions.size());         // NOLINT
  assert(lefs.size() == lef_unloader_affinity.size());  // NOLINT

  auto is_lef_bar_collision = [&](const auto i) { return i < barriers.size(); };

  for (auto i = 0UL; i < lefs.size(); ++i) {
    const auto& lef = lefs[i];
    if (!lef.is_bound()) {
      lef_unloader_affinity[i] = 0.0;
    } else if (!is_lef_bar_collision(rev_collisions[i]) || !is_lef_bar_collision(fwd_collisions[i]))
      MODLE_LIKELY { lef_unloader_affinity[i] = 1.0; }
    else {
      const auto& rev_barrier = barriers[rev_collisions[i]];
      const auto& fwd_barrier = barriers[fwd_collisions[i]];

      if (rev_barrier.blocking_direction_major() == dna::rev &&
          fwd_barrier.blocking_direction_major() == dna::fwd)
        MODLE_UNLIKELY { lef_unloader_affinity[i] = 1.0 / this->hard_stall_multiplier; }
      else
        MODLE_LIKELY { lef_unloader_affinity[i] = 1.0; }
    }
  }
}

void Simulation::select_lefs_to_release(absl::Span<std::size_t> lef_idx,
                                        absl::Span<const double> lef_unloader_affinity,
                                        modle::PRNG_t& rand_eng) noexcept(utils::ndebug_defined()) {
  std::discrete_distribution<std::size_t> idx_gen(lef_unloader_affinity.begin(),
                                                  lef_unloader_affinity.end());
  std::generate(lef_idx.begin(), lef_idx.end(), [&]() { return idx_gen(rand_eng); });
}

void Simulation::release_lefs(absl::Span<Lef> lefs,
                              absl::Span<const std::size_t> lef_idx) noexcept {
  for (const auto i : lef_idx) {
    assert(i < lefs.size());  // NOLINT
    lefs[i].release();
  }
}

boost::asio::thread_pool Simulation::instantiate_thread_pool() const {
  return boost::asio::thread_pool(this->nthreads);
}

template <typename I>
boost::asio::thread_pool Simulation::instantiate_thread_pool(I nthreads, bool clamp_nthreads) {
  static_assert(std::is_integral_v<I>, "nthreads should have an integral type.");
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  if (clamp_nthreads) {
    return boost::asio::thread_pool(
        std::min(std::thread::hardware_concurrency(), static_cast<unsigned int>(nthreads)));
  }
  assert(nthreads > 0);
  return boost::asio::thread_pool(static_cast<unsigned int>(nthreads));
  DISABLE_WARNING_POP
}

Simulation::State& Simulation::State::operator=(const Task& task) {
  this->id = task.id;
  this->chrom = task.chrom;
  this->cell_id = task.cell_id;
  this->n_target_epochs = task.n_target_epochs;
  this->n_target_contacts = task.n_target_contacts;
  this->nlefs = task.nlefs;
  this->barriers = task.barriers;
  return *this;
}

void Simulation::State::resize(std::size_t size) {
  if (size == std::numeric_limits<std::size_t>::max()) {
    size = this->nlefs;
  }
  lef_buff.resize(size);
  lef_unloader_affinity.resize(size);
  rank_buff1.resize(size);
  rank_buff2.resize(size);
  moves_buff1.resize(size);
  moves_buff2.resize(size);
  idx_buff1.resize(size);
  idx_buff2.resize(size);
  epoch_buff.resize(size);
}

void Simulation::State::reset() {  // TODO figure out which resets are redundant
  std::for_each(lef_buff.begin(), lef_buff.end(), [](auto& lef) { lef.reset(); });
  std::fill(lef_unloader_affinity.begin(), lef_unloader_affinity.end(), 0.0);
  std::iota(rank_buff1.begin(), rank_buff1.end(), 0);
  std::copy(rank_buff1.begin(), rank_buff1.end(), rank_buff2.begin());
  std::fill(moves_buff1.begin(), moves_buff1.end(), 0);
  std::fill(moves_buff2.begin(), moves_buff2.end(), 0);
  std::fill(idx_buff1.begin(), idx_buff1.end(), 0);
  std::fill(idx_buff2.begin(), idx_buff2.end(), 0);
  std::fill(epoch_buff.begin(), epoch_buff.end(), 0);
}

}  // namespace modle
