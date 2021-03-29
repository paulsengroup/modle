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
#include "modle/chr_sizes.hpp"                   // for ChrSize, Parser
#include "modle/common.hpp"                      // for Bp, fwd, rev
#include "modle/config.hpp"                      // for Config
#include "modle/contacts.hpp"                    // for ContactMatrix
#include "modle/cooler.hpp"                      // for Cooler, Cooler::...
#include "modle/dna.hpp"                         // for Chromosome
#include "modle/extrusion_barriers.hpp"          // for ExtrusionBarrier
#include "modle/extrusion_factors.hpp"           // for Lef, ExtrusionUnit
#include "modle/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP
#include "modle/utils.hpp"                       // for traced

#ifndef BOOST_STACKTRACE_USE_NOOP
#include <boost/exception/get_error_info.hpp>  // for get_error_info
#include <boost/stacktrace/stacktrace.hpp>     // for operator<<
#include <ios>                                 // IWYU pragma: keep for streamsize
#endif

#ifdef USE_XOSHIRO
#include <Xoshiro-cpp/XoshiroCpp.hpp>  // for XoshiroCpp::Xoshiro256PlusPlus XoshiroCpp::SplitMix64
#endif

namespace modle {

#ifdef USE_XOSHIRO
using PRNG = XoshiroCpp::Xoshiro256PlusPlus;
using seeder = XoshiroCpp::SplitMix64;
#else
using PRNG = std::mt19937_64;
using seeder = std::seed_seq;
#endif

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
  assert(!path_to_chrom_sizes.empty());
  assert(!path_to_extr_barriers.empty());
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
  for (auto&& chrom : modle::chr_sizes::Parser(path_to_chrom_sizes).parse_all()) {
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
    if (b.strand == '+' || b.strand == '-') {
      const auto pos = (b.chrom_start + b.chrom_end + 1) / 2;
      // TODO figure out how to read both transition probabilities
      // const auto pblock = b.score != 0 ? b.score : this->probability_of_extrusion_barrier_block;
      barriers.emplace_back(pos, this->ctcf_occupied_self_prob, this->ctcf_not_occupied_self_prob,
                            b.strand);
    } else {
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

  // TODO with the current implementation using a thread_pool is probably overkill
  auto tpool = Simulation::instantiate_thread_pool(this->nthreads + 1, false);

  /* These are the threads spawned by simulate_extrusion:
   * - 1 thread to write contacts to disk. This thread pops Chromosome* from a std::deque once the
   *   simulation on the Chromosome* has been ultimated. This thread is also responsible of
   * freeing memory as soon as it is not needed by any other thread.
   * - The main thread (i.e. this thread), which loops over chromosomes, instantiates data
   *   structures that are shared across simulation threads (e.g. the vector of extr. barriers),
   * then it submits simulation tasks to a concurrent queue
   * - N simulation threads. These threads pop tasks from the concurrent queue fed by the main
   *   thread and carry out the actual simulation
   *
   * Switching to this thread-layout allowed us to significantly reduce memory allocations, as by
   * running few simulation threads for the whole simulation allows us to recycle buffers more
   * efficiently. Most of the time, the only data moved are tasks, which are structs consisting of
   * few numbers, a pointer and a Span (i.e. ptr + size)
   */

  std::atomic<bool> end_of_simulation = false;
  std::mutex progress_mutex;
  std::deque<std::pair<Chromosome*, std::size_t>> progress_queue;

  // Barriers are shared across threads that are simulating loop-extrusion on a given chromosome.
  // Extrusion barriers are placed on this std::deque, and are passed to simulation threads
  // through a Span of const ExtrusionBarriers
  std::mutex barrier_mutex;
  std::deque<std::unique_ptr<std::vector<ExtrusionBarrier>>> barriers;

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
        std::scoped_lock l(progress_mutex);
        if (auto& [chrom, count] = progress_queue.front(); chrom == nullptr) {
          end_of_simulation = true;
          return;
        } else if (count == ncells) {
          chrom_to_be_written = chrom;
          progress_queue.pop_front();
        } else {
          assert(count < ncells);
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
    boost::asio::post(tpool, [&, ctok = std::move(ctok),
                              progress = progress_queue.begin()]() mutable {
      fmt::print(stderr, FMT_STRING("Spawning simulation thread {}...\n"),
                 std::this_thread::get_id());
      try {
        std::array<Task, task_batch_size> task_buff;  // Tasks are dequeue in batch.
                                                      // This is to reduce contention when
                                                      // accessing the queue

        // Buffers declared here are recycled by every simulation kernel launched from this
        // thread Performance gain mostly come from reusing lef_buff. The impact of recycling
        // the other buffers is barely noticeable
        Simulation::State state_buff;

        while (true) {  // Try to dequeue a batch of tasks
          const auto avail_tasks = task_queue.wait_dequeue_bulk_timed(
              ctok, task_buff.begin(), task_buff.size(), std::chrono::milliseconds(10));
          if (avail_tasks == 0) {     // Dequeue timed-out before any task became available
            if (end_of_simulation) {  // Check whether this is because all tasks have been
                                      // processed
              return;
            }
            continue;
          }

          // Loop over new tasks
          auto tasks = absl::MakeSpan(task_buff.data(), avail_tasks);
          for (auto& task : tasks) {
            // Resize and reset buffers
            state_buff = task;
            state_buff.resize();  // resize to nlefs for the current task
            state_buff.barrier_mask.resize(task.barriers.size());
            state_buff.reset();

            if (task.cell_id == 0) {
              // Print a status update when we are processing cell #0 for a given chromosome
              const auto target_iterations = static_cast<std::size_t>(
                  std::max(1.0, std::round((this->target_contact_density *
                                            static_cast<double>(task.chrom->contacts().npixels())) /
                                           (this->lef_fraction_contact_sampling *
                                            static_cast<double>(this->ncells * task.nlefs)))));
              fmt::print(stderr, FMT_STRING("Simulating ~{} epochs for '{}' across {} cells...\n"),
                         target_iterations, task.chrom->name(), ncells);
            }

            // Start the simulation kernel
            Simulation::simulate_extrusion_kernel(state_buff);

            // Update progress
            std::scoped_lock l1(progress_mutex);
            if (progress == progress_queue.end() || progress->first != task.chrom) {
              progress = std::find_if(progress_queue.begin(), progress_queue.end(),
                                      [&](const auto& p) { return task.chrom == p.first; });
              assert(progress != progress_queue.end());
            }

            // We are done simulating loop-extrusion on task.chrom
            if (++progress->second == ncells) {
              // Print a status update and deallocate extr. barriers
              fmt::print(stderr, "Simulation for '{}' successfully completed.\n",
                         task.chrom->name());
              std::scoped_lock l2(barrier_mutex);
              barriers.pop_front();
            }
          }
        }
      } catch (const std::exception& err) {
        // This is needed, as exceptions don't seem to always propagate to the main()
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

  /*
   * The remaining code submits simulation tasks to the queue. Then it waits until all the tasks
   * have been completed, and contacts have been written to disk
   */

  // Create a vector of pointers to chroms, and sort the pointers by chrom size
  std::vector<Chromosome*> chroms(static_cast<std::size_t>(this->_chromosomes.size()));
  std::transform(this->_chromosomes.begin(), this->_chromosomes.end(), chroms.begin(),
                 [](auto& chrom) { return &chrom; });
  std::sort(chroms.begin(), chroms.end(),
            [](const auto* const c1, const auto* const c2) { return c1->size() > c2->size(); });

  absl::Span<const ExtrusionBarrier> extr_barriers_buff{};
  std::vector<Task> tasks(ncells);
  std::size_t taskid = 0;

  // Loop over chromosomes
  for (auto* chrom : chroms) {
    // Don't simulate KO chroms (but write them to disk if the user specified so)
    if (!chrom->ok()) {
      fmt::print(stderr, "SKIPPING '{}'...\n", chrom->name());
      if (this->write_contacts_for_ko_chroms) {
        chrom->allocate_contacts(this->bin_size, this->diagonal_width);
        std::scoped_lock l(barrier_mutex, progress_mutex);
        progress_queue.emplace_back(chrom, ncells);
        barriers.emplace_back(nullptr);
      }
      continue;
    }

    // Allocate the contact matrix. Once it's not needed anymore, the contact matrix will be
    // de-allocated by the thread that is writing contacts to disk
    chrom->allocate_contacts(this->bin_size, this->diagonal_width);
    {
      std::scoped_lock l(barrier_mutex, progress_mutex);
      progress_queue.emplace_back(chrom, 0UL);  // Signal that we have started processing chrom

      // Allocate extr. barriers mapping on the chromosome that is being simulated
      barriers.emplace_back(
          std::make_unique<std::vector<ExtrusionBarrier>>(Simulation::allocate_barriers(chrom)));
      extr_barriers_buff = absl::MakeConstSpan(barriers.back()->data(), barriers.back()->size());
    }

    // Compute # of LEFs to be simulated based on chrom. sizes
    const auto nlefs = static_cast<std::size_t>(std::round(
        this->number_of_lefs_per_mbp * (static_cast<double>(chrom->simulated_size()) / 1.0e6)));

    tasks.resize(ncells);  // This resize is not needed at the moment. It's there as a safe-guard,
    // just in case we decide to allow the number of simulated cells to
    // change across different chromosomes

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
    std::generate(tasks.begin(), tasks.end(), [&, i = 0UL]() mutable {
      return Task{taskid++, chrom, i++, target_epochs, target_contacts, nlefs, extr_barriers_buff};
    });

    while (  // Enqueue the tasks
        !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), tasks.size())) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  }

  // Signal end of simulation to the thread that is writing contacts to disk
  progress_queue.emplace_back(nullptr, 0UL);
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

  s.seed = this->seed + std::hash<std::string_view>{}(s.chrom->name()) +
           std::hash<std::size_t>{}(s.chrom->size()) + std::hash<std::size_t>{}(s.cell_id);
  {
#ifdef USE_XOSHIRO
    modle::seeder seeder_(s.seed);
    s.rand_eng = modle::PRNG(seeder_.generateSeedSequence<4>());
#else
    modle::seeder seeder_{s.seed};
    s.rand_eng = modle::PRNG(seeder_);
#endif
  }

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

  const auto n_burnin_epochs = this->skip_burnin ? 0 : lef_initial_loading_epoch.front();
  s.n_target_epochs += s.n_target_contacts != 0 ? 0 : n_burnin_epochs;

  // Compute the avg. # of LEFs to use to sample contact every iterations
  const auto avg_nlefs_to_sample =
      static_cast<double>(s.nlefs) * this->lef_fraction_contact_sampling;

  // Compute the avg. # of LEFs to be released every iterations
  const auto avg_nlefs_to_release =
      static_cast<double>((this->rev_extrusion_speed + this->fwd_extrusion_speed) * s.nlefs) /
      static_cast<double>(this->average_lef_lifetime);

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

  // Generate initial extr. barrier states
  for (auto i = 0UL; i < s.barrier_mask.size(); ++i) {
    s.barrier_mask[i] =
        std::bernoulli_distribution{this->probability_of_extrusion_barrier_block}(s.rand_eng);
  }

  std::size_t nlefs_to_sample;
  std::size_t nlefs_to_release;

  // Start the burnin phase (followed by the actual simulation)
  for (auto epoch = 0UL; epoch < s.n_target_epochs; ++epoch) {
    if (!lef_initial_loading_epoch.empty()) {  // Execute this branch only during the burnin phase
      if (this->skip_burnin) {
        lef_initial_loading_epoch.clear();
      }
      while (!lef_initial_loading_epoch.empty() && lef_initial_loading_epoch.back() == epoch) {
        // Consume epochs for LEFs that are supposed to be loaded in the current epoch
        lef_initial_loading_epoch.pop_back();
      }
      // Don't remove this assertion. It is very useful to prevent empty Spans (which cause weird
      // issues later in the simulation) NOLINTNEXTLINE
      assert(lef_initial_loading_epoch.size() < s.nlefs);

      // Compute the current number of active LEFs (i.e. LEFs that are either bound to DNA, or
      // that are valid candidates for re-binding) and grow spans accordingly
      const auto nlefs = s.nlefs - lef_initial_loading_epoch.size();

      lefs = absl::MakeSpan(s.lef_buff.data(), nlefs);
      lef_unloader_affinity = absl::MakeSpan(s.lef_unloader_affinity.data(), nlefs);
      rev_lef_ranks = absl::MakeSpan(s.rank_buff1.data(), nlefs);
      fwd_lef_ranks = absl::MakeSpan(s.rank_buff2.data(), nlefs);
      rev_moves = absl::MakeSpan(s.moves_buff1.data(), nlefs);
      fwd_moves = absl::MakeSpan(s.moves_buff2.data(), nlefs);
      rev_collision_mask = absl::MakeSpan(s.idx_buff1.data(), nlefs);
      fwd_collision_mask = absl::MakeSpan(s.idx_buff2.data(), nlefs);
      nlefs_to_release = std::poisson_distribution<std::size_t>{
          static_cast<double>((this->rev_extrusion_speed + this->fwd_extrusion_speed) * s.nlefs) /
          static_cast<double>(this->average_lef_lifetime)}(s.rand_eng);
    } else {
      nlefs_to_release = std::poisson_distribution<std::size_t>{avg_nlefs_to_release}(s.rand_eng);
    }

    // Simulate one epoch
    {
      auto lef_mask = absl::MakeSpan(s.idx_buff1.data(), lefs.size());
      this->select_lefs_to_bind(lefs, lef_mask);
      this->bind_lefs(s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, lef_mask, s.rand_eng);
    }

    if (epoch > n_burnin_epochs) {
      assert(fwd_lef_ranks.size() == s.nlefs);
      nlefs_to_sample = std::poisson_distribution<std::size_t>{avg_nlefs_to_sample}(s.rand_eng);
      if (s.n_target_contacts != 0) {
        nlefs_to_sample = std::min(nlefs_to_sample, s.n_target_contacts - n_contacts);
      }
      auto lef_idx = absl::MakeSpan(s.idx_buff1.data(), nlefs_to_sample);
      std::sample(fwd_lef_ranks.begin(), fwd_lef_ranks.end(), lef_idx.begin(), lef_idx.size(),
                  s.rand_eng);
      n_contacts += this->register_contacts(s.chrom, lefs, lef_idx);
      if (s.n_target_contacts != 0 && n_contacts >= s.n_target_contacts) {
        return;
      }
    }

    this->generate_moves(s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                         s.rand_eng);

    this->generate_ctcf_states(barriers, s.barrier_mask, s.rand_eng);
    std::fill(rev_collision_mask.begin(), rev_collision_mask.end(), NO_COLLISION);
    std::fill(fwd_collision_mask.begin(), fwd_collision_mask.end(), NO_COLLISION);

    this->check_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                                   barriers, s.barrier_mask, rev_collision_mask, fwd_collision_mask,
                                   s.rand_eng);

    this->check_lef_lef_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                                   rev_collision_mask, fwd_collision_mask, s.rand_eng);

    this->extrude(s.chrom, lefs, rev_moves, fwd_moves);

    this->generate_lef_unloader_affinities(lefs, barriers, rev_collision_mask, fwd_collision_mask,
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
                           modle::PRNG& rand_eng) {
  static_assert(std::is_same_v<boost::dynamic_bitset<>, MaskT> ||
                    std::is_integral_v<std::decay_t<decltype(std::declval<MaskT&>().operator[](
                        std::declval<std::size_t>()))>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");

  assert(lefs.size() <= mask.size() || mask.empty());
  assert(chrom);
  chrom_pos_generator_t pos_generator{chrom->start_pos(), chrom->end_pos()};
  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (mask.empty() || mask[i]) {
      lefs[i].bind_at_pos(pos_generator(rand_eng));
    }
  }
  this->rank_lefs(lefs, rev_lef_ranks, fwd_lef_ranks);
}

inline void Simulation::generate_ctcf_states(absl::Span<const ExtrusionBarrier> extr_barriers,
                                             boost::dynamic_bitset<>& mask, modle::PRNG& rand_eng) {
  assert(extr_barriers.size() == mask.size());
  for (auto i = 0UL; i < extr_barriers.size(); ++i) {
    mask[i] = CTCF::next_state(mask[i] ? CTCF::OCCUPIED : CTCF::NOT_OCCUPIED,
                               this->ctcf_occupied_self_prob, this->ctcf_not_occupied_self_prob,
                               rand_eng);
  }
}

void Simulation::generate_moves(const Chromosome* const chrom, absl::Span<const Lef> lefs,
                                absl::Span<const std::size_t> rev_lef_ranks,
                                absl::Span<const std::size_t> fwd_lef_ranks,
                                absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
                                modle::PRNG& rand_eng) {
  assert(lefs.size() == fwd_lef_ranks.size());  // NOLINT
  assert(lefs.size() == rev_lef_ranks.size());  // NOLINT
  assert(lefs.size() == fwd_moves.size());      // NOLINT
  assert(lefs.size() == rev_moves.size());      // NOLINT
  auto generate_rev_move = [&](std::size_t i) {
    const auto& rev_unit = lefs[i].rev_unit;
    assert(rev_unit.pos() >= chrom->start_pos());  // NOLINT
    if (this->rev_extrusion_speed_std == 0) {      // When std == 0 always return the mean
      return std::min(this->rev_extrusion_speed, rev_unit.pos() - chrom->start_pos());
    }
    // Generate the move distance (and make sure it is not a negative distance)
    // NOTE: on my laptop generating doubles from a normal distribution, rounding them, then
    // casting double to uint is a lot faster (~4x) than drawing uints directly from a Poisson
    // distr.
    return std::min(static_cast<bp_t>(std::round(
                        lef_move_generator_t{static_cast<double>(this->rev_extrusion_speed),
                                             this->rev_extrusion_speed_std}(rand_eng))),
                    rev_unit.pos() - chrom->start_pos());
  };

  auto generate_fwd_move =  // Generate the move distance (and clamp it, so that the move does not
      [&](std::size_t i) {  // result in moving the extrusion unit past the chrom. end)
        const auto& fwd_unit = lefs[i].fwd_unit;
        assert(fwd_unit.pos() < chrom->end_pos());  // NOLINT
        if (this->fwd_extrusion_speed_std == 0) {   // When std == 0 always return the mean
          return std::min(this->fwd_extrusion_speed, (chrom->end_pos() - 1) - fwd_unit.pos());
        }
        // NOTE: on my laptop generating doubles from a normal
        // distribution, rounding them, then casting double to uint is a lot faster (~4x) than
        // drawing uints directly from a Poisson distr.
        return std::min(static_cast<bp_t>(std::round(
                            lef_move_generator_t{static_cast<double>(this->fwd_extrusion_speed),
                                                 this->fwd_extrusion_speed_std}(rand_eng))),
                        (chrom->end_pos() - 1) - fwd_unit.pos());
      };

  // As long as a LEF is bound to DNA, always generate a move, even in case a LEF has one or more
  // of its extr. unit stalled. In the latter case, the number of moves will be used to decrease
  // the number of stalls
  for (auto i = 0UL; i < lefs.size(); ++i) {
    rev_moves[i] = lefs[i].is_bound() ? generate_rev_move(i) : 0UL;
    fwd_moves[i] = lefs[i].is_bound() ? generate_fwd_move(i) : 0UL;
  }
  this->adjust_moves(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves);
}

void Simulation::adjust_moves(const Chromosome* chrom, absl::Span<const Lef> lefs,
                              absl::Span<const std::size_t> rev_lef_ranks,
                              absl::Span<const std::size_t> fwd_lef_ranks,
                              absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves) {
  // Make sure that consecutive extr. units that are moving in the same direction do not bypass
  // each other. Consider the following example: LEF1.fwd_unit.pos() = 0; LEF2.fwd_unit.pos() =
  // 10; fwd_extrusion_speed = 1000; fwd_extrusion_speed_std > 0; the fwd_unit of LEF1 is set to
  // move 1050 bp in the current iteration, while LEF2 is set to move only 950 bp. In this
  // scenario, the fwd_unit of LEF1 would bypass the fwd_unit of LEF2. In a real system, what
  // would likely happen, is that the fwd_unit of LEF1 would push the fwd_unit of LEF2,
  // temporarily increasing the fwd extr. speed of LEF2. The following code does exactly what was
  // described for the real-system scenario.
  (void)chrom;
  // Loop over pairs of consecutive LEFs in 3'-5' direction
  const auto rev_offset = lefs.size() - 1;
  for (auto i = 0UL; i < lefs.size() - 1; ++i) {
    if (!lefs[i].is_bound()) {
      continue;
    }
    const auto& idx1 = rev_lef_ranks[rev_offset - 1 - i];
    const auto& idx2 = rev_lef_ranks[rev_offset - i];

    assert(lefs[idx1].rev_unit.pos() >= chrom->start_pos() + rev_moves[idx1]);
    assert(lefs[idx2].rev_unit.pos() >= chrom->start_pos() + rev_moves[idx2]);
    const auto pos1 = lefs[idx1].rev_unit.pos() - rev_moves[idx1];
    const auto pos2 = lefs[idx2].rev_unit.pos() - rev_moves[idx2];
    if (pos2 < pos1) {
      rev_moves[idx1] += pos1 - pos2;
    }
    const auto& idx3 = fwd_lef_ranks[i];
    const auto& idx4 = fwd_lef_ranks[i + 1];

    assert(lefs[idx3].fwd_unit.pos() + fwd_moves[idx3] < chrom->end_pos());
    assert(lefs[idx4].fwd_unit.pos() + fwd_moves[idx4] < chrom->end_pos());
    const auto pos3 = lefs[idx3].fwd_unit.pos() + fwd_moves[idx3];
    const auto pos4 = lefs[idx4].fwd_unit.pos() + fwd_moves[idx4];
    if (pos3 > pos4) {
      fwd_moves[idx4] += pos3 - pos4;
    }
  }
}

void Simulation::rank_lefs(absl::Span<const Lef> lefs, absl::Span<std::size_t> rev_lef_rank_buff,
                           absl::Span<std::size_t> fwd_lef_rank_buff, bool init_buffers) {
  if (init_buffers) {  // Init rank buffers
    std::iota(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), 0);
    std::copy(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), rev_lef_rank_buff.begin());
  }
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());

  if (lefs.size() < 20) {  // TODO Come up with a reasonable threshold
    cppsort::insertion_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                            [&](const auto r1, const auto r2) {
                              return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                            });

    cppsort::insertion_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                            [&](const auto r1, const auto r2) {
                              return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                            });
  } else {
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
    // times we have 2.5-5% inversions, and very rarely > 10%, which is drop merge sort shines).
    // https://github.com/Morwenn/cpp-sort/blob/develop/docs/Benchmarks.md#inv-adaptive-algorithms
    // https://github.com/Morwenn/cpp-sort/blob/develop/docs/Sorters.md#drop_merge_sorter
    cppsort::drop_merge_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                             [&](const auto r1, const auto r2) {
                               return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                             });

    cppsort::drop_merge_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                             [&](const auto r1, const auto r2) {
                               return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                             });
  }
}

void Simulation::extrude(const Chromosome* chrom, absl::Span<Lef> lefs,
                         absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves) {
  assert(lefs.size() == rev_moves.size());  // NOLINT
  assert(lefs.size() == fwd_moves.size());  // NOLINT
  (void)chrom;

  for (auto i = 0UL; i < lefs.size(); ++i) {
    auto& lef = lefs[i];
    if (!lef.is_bound()) {  // LEF does not currently bind DNA
      continue;
    }

    // Extrude rev unit
    assert(lef.rev_unit.pos() >= chrom->start_pos() + rev_moves[i]);  // NOLINT
    lef.rev_unit._pos -= rev_moves[i];  // Advance extr. unit in 3'-5' direction

    // Extrude fwd unit
    assert(lef.fwd_unit.pos() + fwd_moves[i] <= chrom->end_pos() - 1);
    lef.fwd_unit._pos += fwd_moves[i];  // Advance extr. unit in 5'-3' direction
  }
}

void Simulation::check_lef_bar_collisions(
    absl::Span<const Lef> lefs, absl::Span<const std::size_t> rev_lef_rank_buff,
    absl::Span<const std::size_t> fwd_lef_rank_buff, absl::Span<bp_t> rev_move_buff,
    absl::Span<bp_t> fwd_move_buff, absl::Span<const ExtrusionBarrier> extr_barriers,
    const boost::dynamic_bitset<>& barrier_mask, absl::Span<std::size_t> rev_collisions,
    absl::Span<std::size_t> fwd_collisions, modle::PRNG& rand_eng) {
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());
  assert(lefs.size() == fwd_move_buff.size());
  assert(lefs.size() == rev_move_buff.size());
  assert(lefs.size() == fwd_collisions.size());
  assert(lefs.size() == rev_collisions.size());
  assert(barrier_mask.size() == extr_barriers.size());
  assert(std::is_sorted(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                        }));
  assert(std::is_sorted(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                        }));
  assert(std::all_of(rev_collisions.begin(), rev_collisions.end(),
                     [](const auto c) { return c == NO_COLLISION; }));
  assert(std::all_of(fwd_collisions.begin(), fwd_collisions.end(),
                     [](const auto c) { return c == NO_COLLISION; }));

  /* Loop over lefs, using a procedure similar to merge in mergesort
   * The idea here is that if we have a way to visit extr. units in 5'-3' order, then detecting
   * LEF-LEF collisions boils down to:
   *  - Starting from the first fwd extr. unit in 5'-3' order:
   *    - Look for rev extr. units that are located downstream of the fwd unit that is being
   *      processed
   *    - Continue looking for the next rev unit, until the distance between the fwd and rev units
   *      being considered is less than a certain threshold (e.g. the bin size)
   *    - While doing so, increase the number of LEF-LEF collisions for the fwd/rev unit that are
   *      being processed
   * The fwd and rev_collision_buff will contain the index corresponding to the ext. barrier that
   * caused the collision, or -1 in case of no collisions
   */
  auto rev_idx = rev_lef_rank_buff.front();
  auto fwd_idx = fwd_lef_rank_buff.front();
  auto rev_unit_pos = lefs[rev_idx].rev_unit.pos();
  auto fwd_unit_pos = lefs[fwd_idx].fwd_unit.pos();
  std::size_t j1 = 0, j2 = 0;

  // Loop over extr. barriers and find the first, possibly colliding extr. unit
  for (auto i = 0UL; i < extr_barriers.size(); ++i) {
    if (barrier_mask[i] == CTCF::NOT_OCCUPIED) {
      continue;
    }
    const auto& barrier = extr_barriers[i];

    if (j1 < lefs.size()) {  // Process rev unit
      const auto& pblock = barrier.blocking_direction() == dna::rev
                               ? this->lef_hard_collision_pblock
                               : this->lef_soft_collision_pblock;
      // Look for the first rev extr. unit that comes after the current barrier
      while (rev_unit_pos <= barrier.pos()) {
        if (++j1 == lefs.size()) {  // All rev units have been processed
          goto process_fwd_unit;
        }
        rev_idx = rev_lef_rank_buff[j1];
        rev_unit_pos = lefs[rev_idx].rev_unit.pos();
      }

      auto& rev_move = rev_move_buff[rev_idx];
      if (const auto delta = rev_unit_pos - barrier.pos();
          delta < rev_move && std::bernoulli_distribution{pblock}(rand_eng)) {
        // Collision detected. Assign barrier idx to the respective entry in the collision mask
        rev_collisions[rev_idx] = i;
        // Move LEF close to the extr. barrier
        rev_move = delta > 1 ? delta - 1 : 0;
      }
    }

  process_fwd_unit:
    if (j2 < lefs.size()) {  // Process fwd unit
      const auto& pblock = barrier.blocking_direction() == dna::fwd
                               ? this->lef_hard_collision_pblock
                               : this->lef_soft_collision_pblock;
      // Look for the first fwd extr. unit that comes before the current barrier
      while (fwd_unit_pos <= barrier.pos()) {
        if (++j2 == lefs.size()) {  // All fwd units have been processed
          goto end_of_loop;
        }
        fwd_idx = fwd_lef_rank_buff[j2];
        fwd_unit_pos = lefs[fwd_idx].fwd_unit.pos();
      }
      fwd_idx = fwd_lef_rank_buff[j2 > 0 ? j2 - 1 : 0];
      fwd_unit_pos = lefs[fwd_idx].fwd_unit.pos();

      auto& fwd_move = fwd_move_buff[fwd_idx];
      if (const auto delta = barrier.pos() - fwd_unit_pos;
          delta < fwd_move && std::bernoulli_distribution{pblock}(rand_eng)) {
        // Collision detected. Assign barrier idx to the respective entry in the collision mask
        fwd_collisions[fwd_idx] = i;
        // Move LEF close to the extr. barrier
        fwd_move = delta > 1 ? delta - 1 : 0;
      }
    }

  end_of_loop:
    // Return if there are no more extr.units to be processed
    if (j1 == extr_barriers.size() && j2 == 0) {
      return;
    }
  }
}

void Simulation::check_lef_lef_collisions(
    absl::Span<const Lef> lefs, absl::Span<const std::size_t> rev_lef_rank_buff,
    absl::Span<const std::size_t> fwd_lef_rank_buff, absl::Span<bp_t> rev_move_buff,
    absl::Span<bp_t> fwd_move_buff, absl::Span<std::size_t> rev_collision_mask,
    absl::Span<std::size_t> fwd_collision_mask, modle::PRNG& rand_eng) {
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());
  assert(lefs.size() == fwd_collision_mask.size());
  assert(lefs.size() == rev_collision_mask.size());
  assert(std::is_sorted(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                        }));
  assert(std::is_sorted(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                        }));

  /* Loop over lefs, using a procedure similar to merge in mergesort
   * The idea here is that if we have a way to visit extr. units in 5'-3' order, then detecting
   * LEF-LEF collisions boils down to:
   *  - Starting from the first fwd extr. unit in 5'-3' order:
   *    - Look for rev extr. units that are located downstream of the fwd unit that is being
   *      processed
   *    - Continue looking for the next rev unit, until the distance between the fwd and rev units
   *      being considered is less than a certain threshold (e.g. the bin size)
   *    - While doing so, increase the number of LEF-LEF collisions for the fwd/rev unit that are
   *      being processed
   */
  for (auto i = 0UL, j = 0UL; i < lefs.size() && j < lefs.size(); ++i) {
    auto rev_idx = rev_lef_rank_buff[j];          // index of the ith rev unit in 5'-3' order
    auto rev_pos = lefs[rev_idx].rev_unit.pos();  // pos of the ith unit
    auto fwd_idx = fwd_lef_rank_buff[i];          // index of the jth fwd unit in 5'-3' order
    auto fwd_pos = lefs[fwd_idx].fwd_unit.pos();  // pos of the jth unit

    while (rev_pos <= fwd_pos) {
      if (++j == lefs.size()) {
        goto process_fwd_unit;
      }
      rev_idx = rev_lef_rank_buff[j];          // index of the jth fwd unit in 5'-3' order
      rev_pos = lefs[rev_idx].rev_unit.pos();  // pos of the jth unit
    }

    if (const auto delta = rev_pos - fwd_pos;
        delta < rev_move_buff[rev_idx] + fwd_move_buff[fwd_idx] &&
        std::bernoulli_distribution{1.0 - this->probability_of_extrusion_unit_bypass}(rand_eng)) {
      if (delta > 1) {
        if (rev_collision_mask[rev_idx] == NO_COLLISION &&
            fwd_collision_mask[fwd_idx] == NO_COLLISION) {
          const auto collision_pos = (fwd_pos + rev_pos + 1) / 2;
          rev_move_buff[rev_idx] = rev_pos - collision_pos;
          fwd_move_buff[fwd_idx] = (collision_pos - fwd_pos) - 1;
        } else if (rev_collision_mask[rev_idx] != NO_COLLISION &&
                   fwd_collision_mask[fwd_idx] == NO_COLLISION) {
          rev_move_buff[rev_idx] = 0;
          fwd_move_buff[fwd_idx] = (rev_pos - fwd_pos) - 1;
        } else if (rev_collision_mask[rev_idx] == NO_COLLISION &&
                   fwd_collision_mask[fwd_idx] != NO_COLLISION) {
          rev_move_buff[rev_idx] = (rev_pos - fwd_pos) - 1;
          fwd_move_buff[fwd_idx] = 0;
        }
      } else {
        rev_move_buff[rev_idx] = 0;
        fwd_move_buff[fwd_idx] = 0;
      }
      rev_collision_mask[rev_idx] = LEF_LEF_COLLISION;
      fwd_collision_mask[fwd_idx] = LEF_LEF_COLLISION;
    }
  }

process_fwd_unit:
  for (auto i = 1UL; i < lefs.size(); ++i) {
    const auto& fwd_idx2 = fwd_lef_rank_buff[i];  // index of the ith fwd unit in 5'-3' order
    const auto& fwd_pos2 = lefs[fwd_idx2].fwd_unit.pos();  // pos of the ith unit
    if (fwd_collision_mask[fwd_idx2] == NO_COLLISION) {
      continue;
    }

    for (auto j = i - 1; j > 0; --j) {
      const auto& fwd_idx1 = fwd_lef_rank_buff[j];  // index of the ith-1 fwd unit in 5'-3' order
      const auto& fwd_pos1 = lefs[fwd_idx1].fwd_unit.pos();  // pos of the ith-1 unit
      auto& move1 = fwd_move_buff[fwd_idx1];
      const auto& move2 = fwd_move_buff[fwd_idx2];

      if (fwd_pos2 - fwd_pos1 <= move1 + move2 &&
          std::bernoulli_distribution{1.0 - this->probability_of_extrusion_unit_bypass}(rand_eng)) {
        fwd_collision_mask[fwd_idx1] = LEF_LEF_COLLISION;
        move1 = (fwd_pos2 + move2) - fwd_pos1;
        move1 -= move1 > 0 ? 1 : 0;
      } else {
        break;
      }
    }
  }

  for (auto i = 0UL; i < lefs.size() - 1; ++i) {
    const auto& rev_idx1 = rev_lef_rank_buff[i];  // index of the ith rev unit in 5'-3' order
    const auto& rev_pos1 = lefs[rev_idx1].rev_unit.pos();  // pos of the ith unit
    if (rev_collision_mask[rev_idx1] == NO_COLLISION) {
      continue;
    }

    const auto& rev_idx2 = rev_lef_rank_buff[i + 1];  // index of the ith-1 rev unit in 5'-3' order
    const auto& rev_pos2 = lefs[rev_idx2].rev_unit.pos();  // pos of the ith-1 unit

    const auto& move1 = rev_move_buff[rev_idx1];
    auto& move2 = rev_move_buff[rev_idx2];
    if (rev_pos2 - rev_pos1 <= move1 + move2 &&
        std::bernoulli_distribution{1.0 - this->probability_of_extrusion_unit_bypass}(rand_eng)) {
      rev_collision_mask[rev_idx2] = LEF_LEF_COLLISION;
      move2 = rev_pos2 - (rev_pos1 - move1);
      move2 -= move2 > 0 ? 1 : 0;
    }
  }
}

std::size_t Simulation::register_contacts(Chromosome* chrom, absl::Span<const Lef> lefs,
                                          absl::Span<const std::size_t> selected_lef_idx) {
  // Register contacts for the selected LEFs (excluding LEFs that have one of their units at the
  // beginning/end of a chromosome)
  std::size_t new_contacts = 0;
  for (const auto i : selected_lef_idx) {
    assert(i < lefs.size());
    const auto& lef = lefs[i];
    if (lef.is_bound() && lef.rev_unit.pos() > chrom->start_pos() &&
        lef.rev_unit.pos() < chrom->end_pos() && lef.fwd_unit.pos() > chrom->start_pos() &&
        lef.fwd_unit.pos() < chrom->end_pos()) {
      chrom->increment_contacts(lef.rev_unit.pos(), lef.fwd_unit.pos(), this->bin_size);
      ++new_contacts;
    }
  }
  return new_contacts;
}

template <typename MaskT>
void Simulation::select_lefs_to_bind(absl::Span<const Lef> lefs, MaskT& mask) {
  static_assert(std::is_same_v<boost::dynamic_bitset<>, MaskT> ||
                    std::is_integral_v<std::decay_t<decltype(std::declval<MaskT&>().operator[](
                        std::declval<std::size_t>()))>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");
  assert(lefs.size() == mask.size());
  std::transform(lefs.begin(), lefs.end(), mask.begin(),
                 [](const auto& lef) { return !lef.is_bound(); });
}

void Simulation::generate_lef_unloader_affinities(absl::Span<const Lef> lefs,
                                                  absl::Span<const ExtrusionBarrier> barriers,
                                                  absl::Span<const std::size_t> rev_collisions,
                                                  absl::Span<const std::size_t> fwd_collisions,
                                                  absl::Span<double> lef_unloader_affinity) {
  assert(lefs.size() == rev_collisions.size());         // NOLINT
  assert(lefs.size() == fwd_collisions.size());         // NOLINT
  assert(lefs.size() == lef_unloader_affinity.size());  // NOLINT

  for (auto i = 0UL; i < lefs.size(); ++i) {
    const auto& lef = lefs[i];
    if (!lef.is_bound()) {
      lef_unloader_affinity[i] = 0.0;
    } else if (rev_collisions[i] == NO_COLLISION || fwd_collisions[i] == NO_COLLISION ||
               rev_collisions[i] == LEF_LEF_COLLISION || fwd_collisions[i] == LEF_LEF_COLLISION) {
      lef_unloader_affinity[i] = 1.0;
    } else {
      const auto& rev_barrier = barriers[rev_collisions[i]];
      const auto& fwd_barrier = barriers[fwd_collisions[i]];

      if (rev_barrier.blocking_direction() == dna::rev &&
          fwd_barrier.blocking_direction() == dna::fwd) {
        lef_unloader_affinity[i] = 1.0 / this->hard_stall_multiplier;
      } else {
        lef_unloader_affinity[i] = 1.0;
      }
    }
  }
}

void Simulation::select_lefs_to_release(absl::Span<std::size_t> lef_idx,
                                        absl::Span<const double> lef_unloader_affinity,
                                        modle::PRNG& rand_eng) {
  std::discrete_distribution<std::size_t> idx_gen(lef_unloader_affinity.begin(),
                                                  lef_unloader_affinity.end());
  std::generate(lef_idx.begin(), lef_idx.end(), [&]() { return idx_gen(rand_eng); });
}

void Simulation::release_lefs(absl::Span<Lef> lefs, absl::Span<const std::size_t> lef_idx) {
  for (const auto& i : lef_idx) {
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

void Simulation::State::operator=(const Task& task) {
  this->id = task.id;
  this->chrom = task.chrom;
  this->cell_id = task.cell_id;
  this->n_target_epochs = task.n_target_epochs;
  this->n_target_contacts = task.n_target_contacts;
  this->nlefs = task.nlefs;
  this->barriers = task.barriers;
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
