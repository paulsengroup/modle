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

#ifndef BOOST_STACKTRACE_USE_NOOP
#include <boost/exception/get_error_info.hpp>  // for get_error_info
#include <boost/stacktrace/stacktrace.hpp>     // for operator<<
#include <ios>                                 // IWYU pragma: keep for streamsize
#endif

#ifdef USE_XOSHIRO
#include <XoshiroCpp.hpp>
#endif

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
      const auto pos =
          static_cast<Bp>(std::round(static_cast<double>(b.chrom_start + b.chrom_end) / 2.0));
      const auto pblock = b.score != 0 ? b.score : this->probability_of_extrusion_barrier_block;
      barriers.emplace_back(pos, pblock, b.strand);
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

        // Buffers declared here are recycled by every simulation kernel launched from this thread
        // Performance gain mostly come from reusing lef_buff. The impact of recycling the other
        // buffers is barely noticeable
        std::vector<Lef> lef_buff;
        std::vector<std::size_t> rev_lef_rank_buff;
        std::vector<std::size_t> fwd_lef_rank_buff;
        boost::dynamic_bitset<> mask;
        std::vector<collision_t> rev_lef_collision_buff;
        std::vector<collision_t> fwd_lef_collision_buff;

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
            // Resize buffers
            lef_buff.resize(task.nlefs);
            rev_lef_rank_buff.resize(task.nlefs);
            fwd_lef_rank_buff.resize(task.nlefs);
            mask.resize(task.nlefs);
            rev_lef_collision_buff.resize(task.nlefs);
            fwd_lef_collision_buff.resize(task.nlefs);

            // Reset buffers
            std::for_each(lef_buff.begin(), lef_buff.end(), [](auto& lef) { lef.reset(); });
            std::iota(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), 0);
            std::iota(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), 0);
            mask.reset();
            std::fill(rev_lef_collision_buff.begin(), rev_lef_collision_buff.end(), 0);
            std::fill(fwd_lef_collision_buff.begin(), fwd_lef_collision_buff.end(), 0);

            if (task.cell_id == 0) {
              // Print a status update when we are processing cell #0 for a given chromosome
              fmt::print(stderr, FMT_STRING("Simulating {} epochs for '{}' across {} cells...\n"),
                         task.nrounds, task.chrom->name(), ncells);
            }

            // Start the simulation kernel
            Simulation::simulate_extrusion_kernel(
                task.chrom, task.cell_id, task.nrounds, lef_buff, task.barriers,
                absl::MakeSpan(rev_lef_rank_buff), absl::MakeSpan(fwd_lef_rank_buff), mask,
                absl::MakeSpan(rev_lef_collision_buff), absl::MakeSpan(fwd_lef_collision_buff));

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

    if (this->target_contact_density != 0) {  // Compute the number of simulation rounds required
      // to reach the target contact density
      this->simulation_iterations = static_cast<std::size_t>(
          std::max(1.0, std::round((this->target_contact_density *
                                    static_cast<double>(this->contact_sampling_interval *
                                                        chrom->contacts().npixels())) /
                                   static_cast<double>(ncells * nlefs))));
    }

    // Generate a batch of tasks for all the simulations involving the current chrom
    std::generate(tasks.begin(), tasks.end(), [&, i = 0UL]() mutable {
      return Task{taskid++, chrom, i++, this->simulation_iterations, nlefs, extr_barriers_buff};
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

void Simulation::simulate_extrusion_kernel(
    Chromosome* chrom, std::size_t cell_id, std::size_t n_target_epochs, std::vector<Lef> lef_buff,
    const absl::Span<const ExtrusionBarrier> extr_barrier_buff,
    absl::Span<std::size_t> rev_lef_rank_buff, absl::Span<std::size_t> fwd_lef_rank_buff,
    boost::dynamic_bitset<>& mask, absl::Span<collision_t> rev_lef_collision_buff,
    absl::Span<collision_t> fwd_lef_collision_buff) {
  const auto seed_ = this->seed + std::hash<std::string_view>{}(chrom->name()) +
                     std::hash<std::size_t>{}(chrom->size()) + std::hash<std::size_t>{}(cell_id);
#ifdef USE_XOSHIRO
  modle::seeder seeder_(seed_);
  modle::PRNG rand_eng(seeder_.generateSeedSequence<4>());
#else
  modle::seeder seeder_{seed_};
  modle::PRNG rand_eng(seeder_);
#endif

  // Generate the epoch at which each LEF is supposed to be initially loaded
  std::vector<std::size_t> lef_initial_loading_epoch(lef_buff.size());

  // TODO Consider using a Poisson process instead of sampling from an uniform distribution
  std::uniform_int_distribution<std::size_t> round_gen{
      0, (4 * this->average_lef_lifetime) / this->bin_size};
  std::generate(lef_initial_loading_epoch.begin(), lef_initial_loading_epoch.end(),
                [&]() { return round_gen(rand_eng); });

  // Sort epochs in descending order
  if (round_gen.max() > 2048) {
    // Counting sort uses n + r space in memory, where r is the number of unique values in the
    // range to be sorted. For this reason it is not a good idea to use it when the sampling
    // interval is relatively large. Whether 2048 is a reasonable threshold has yet to be tested
    cppsort::ska_sort(lef_initial_loading_epoch.rbegin(), lef_initial_loading_epoch.rend());
  } else {
    cppsort::counting_sort(lef_initial_loading_epoch.rbegin(), lef_initial_loading_epoch.rend());
  }

  // Shift epochs so that the first epoch == 0
  if (const auto offset = lef_initial_loading_epoch.back(); offset != 0) {
    std::for_each(lef_initial_loading_epoch.begin(), lef_initial_loading_epoch.end(),
                  [&](auto& n) { n -= offset; });
  }

  const auto n_burnin_epochs = lef_initial_loading_epoch.front();
  n_target_epochs += n_burnin_epochs;

  std::bernoulli_distribution sample_contacts{1.0 / this->contact_sampling_interval};

  // Declare the spans used by the burnin phase as well as the simulation itself
  auto lefs = absl::MakeSpan(lef_buff.data(), lef_buff.size());
  const auto barriers = absl::MakeConstSpan(extr_barrier_buff);
  auto rev_lef_ranks = absl::MakeSpan(rev_lef_rank_buff);
  auto fwd_lef_ranks = absl::MakeSpan(fwd_lef_rank_buff);
  auto rev_lef_collisions = absl::MakeSpan(rev_lef_collision_buff);
  auto fwd_lef_collisions = absl::MakeSpan(fwd_lef_collision_buff);

  // Start the burnin phase (followed by the actual simulation)
  for (auto epoch = 0UL; epoch < n_target_epochs; ++epoch) {
    const auto register_contacts =
        epoch > n_burnin_epochs &&
        (this->randomize_contact_sampling_interval ? sample_contacts(rand_eng)
                                                   : epoch % this->contact_sampling_interval == 0);

    if (!lef_initial_loading_epoch.empty()) {  // Execute this branch only during the burnin phase
      while (!lef_initial_loading_epoch.empty() && lef_initial_loading_epoch.back() == epoch) {
        // Consume epochs for LEFs that are supposed to be loaded in the current epoch
        lef_initial_loading_epoch.pop_back();
      }
      // Don't remove this assertion. It is very useful to prevent empty Spans (which cause weird
      // issues later in the simulation) NOLINTNEXTLINE
      assert(lef_initial_loading_epoch.size() < lef_buff.size());

      // Compute the current number of active LEFs (i.e. LEFs that are either bound to DNA, or
      // that are valid candidates for re-binding) and grow spans accordingly
      const auto nlefs = lef_buff.size() - lef_initial_loading_epoch.size();

      lefs = absl::MakeSpan(lef_buff.data(), nlefs);
      rev_lef_ranks = absl::MakeSpan(rev_lef_rank_buff.data(), nlefs);
      fwd_lef_ranks = absl::MakeSpan(fwd_lef_rank_buff.data(), nlefs);
      rev_lef_collisions = absl::MakeSpan(rev_lef_collision_buff.data(), nlefs);
      fwd_lef_collisions = absl::MakeSpan(fwd_lef_collision_buff.data(), nlefs);
    }

    // Simulate one epoch
    this->select_lefs_to_bind(lefs, mask);
    this->bind_lefs(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rand_eng, mask);

    this->check_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, barriers, rev_lef_collisions,
                                   fwd_lef_collisions);
    this->apply_lef_bar_stalls(lefs, rev_lef_collisions, fwd_lef_collisions, barriers, rand_eng);

    this->check_lef_lef_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_lef_collisions,
                                   fwd_lef_collisions);
    this->apply_lef_lef_stalls(lefs, rev_lef_collisions, fwd_lef_collisions, rev_lef_ranks,
                               fwd_lef_ranks, rand_eng);

    if (register_contacts) {
      this->register_contacts(chrom, lefs);
    }

    this->extrude(chrom, lefs);
  }
}

template <typename MaskT>
void Simulation::bind_lefs(const Chromosome* const chrom, absl::Span<Lef> lefs,
                           absl::Span<std::size_t> rev_lef_rank_buff,
                           absl::Span<std::size_t> fwd_lef_rank_buff, modle::PRNG& rand_eng,
                           MaskT& mask) {
  static_assert(std::is_same_v<boost::dynamic_bitset<>, MaskT> ||
                    std::is_integral_v<std::decay_t<decltype(std::declval<MaskT&>().operator[](
                        std::declval<std::size_t>()))>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");

  assert(lefs.size() <= mask.size() || mask.empty());
  assert(chrom);
  chrom_pos_generator_t pos_generator{chrom->start_pos(), chrom->end_pos()};
  lef_lifetime_generator_t lifetime_generator{1.0 /
                                              static_cast<double>(this->average_lef_lifetime)};
  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (mask.empty() || mask[i]) {
      lefs[i].bind_at_pos(pos_generator(rand_eng), lifetime_generator(rand_eng));
    }
  }
  this->rank_lefs(lefs, rev_lef_rank_buff, fwd_lef_rank_buff);
}

void Simulation::bind_all_lefs(const Chromosome* chrom, absl::Span<Lef> lefs,
                               absl::Span<std::size_t> rev_lef_rank_buff,
                               absl::Span<std::size_t> fwd_lef_rank_buff, PRNG& rand_eng) {
  chrom_pos_generator_t chrom_pos_generator{chrom->start_pos(), chrom->end_pos()};
  lef_lifetime_generator_t lef_lifetime_generator{1.0 /
                                                  static_cast<double>(this->average_lef_lifetime)};
  absl::btree_multiset<Bp> positions;
  std::generate_n(std::inserter(positions, positions.begin()), lefs.size(),
                  [&]() { return chrom_pos_generator(rand_eng); });
  std::size_t i = 0;
  for (auto& pos : positions) {
    lefs[i].bind_at_pos(pos, lef_lifetime_generator(rand_eng));
    rev_lef_rank_buff[i] = i;
    fwd_lef_rank_buff[i] = i;
    ++i;
  }
}

void Simulation::rank_lefs(absl::Span<Lef> lefs, absl::Span<std::size_t> rev_lef_rank_buff,
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

void Simulation::extrude(const Chromosome* chrom, absl::Span<Lef> lefs) {
  for (auto& lef : lefs) {
    if (lef.lifetime == 0) {  // LEF does not currently bind DNA
      continue;
    }

    // Extrude rev unit
    if (lef.rev_unit.stalled()) {  // Unit is stalled
      lef.rev_unit.decrement_stalls(this->bin_size);
    } else if (lef.rev_unit._pos < chrom->start_pos() + this->bin_size) {  // Unit reached chrom.
                                                                           // boundaries
      lef.rev_unit._pos = chrom->start_pos();
      lef.rev_unit._nstalls_lef_lef = std::numeric_limits<Bp>::max();
    } else {
      lef.rev_unit._pos -= this->bin_size;  // Advance extr. unit in 3'-5' direction
    }

    // Extrude fwd unit
    if (lef.fwd_unit.stalled()) {  // Unit is stalled
      lef.fwd_unit.decrement_stalls(this->bin_size);
    } else if (lef.fwd_unit._pos + this->bin_size > chrom->end_pos() - 1) {  // Unit reached
                                                                             // chrom. boundaries
      lef.fwd_unit._pos = chrom->end_pos() - 1;
      lef.fwd_unit._nstalls_lef_lef = std::numeric_limits<Bp>::max();
    } else {
      lef.fwd_unit._pos += this->bin_size;  // Advance extr. unit in 5'-3' direction
    }

    // Reduce lifetime and release LEFs when appropriate (i.e. LEF lifetime reaches 0)
    if ((lef.lifetime -= std::min(lef.lifetime, 2 * this->bin_size)) == 0) {
      lef.release();
    }
  }
}

void Simulation::check_lef_lef_collisions(absl::Span<const Lef> lefs,
                                          absl::Span<const std::size_t> rev_lef_rank_buff,
                                          absl::Span<const std::size_t> fwd_lef_rank_buff,
                                          absl::Span<collision_t> rev_collision_buff,
                                          absl::Span<collision_t> fwd_collision_buff) {
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());
  assert(lefs.size() == fwd_collision_buff.size());
  assert(lefs.size() == rev_collision_buff.size());
  assert(std::is_sorted(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                        }));
  assert(std::is_sorted(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                        }));

  // Clear buffers
  std::fill(fwd_collision_buff.begin(), fwd_collision_buff.end(), 0);
  std::fill(rev_collision_buff.begin(), rev_collision_buff.end(), 0);

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
  std::size_t i = 0, j = 0;
  while (i < lefs.size() && j < lefs.size()) {
    const auto& fwd_idx = fwd_lef_rank_buff[i];      // index of the ith fwd unit in 5'-3' order
    const auto& fwd = lefs[fwd_idx].fwd_unit.pos();  // pos of the ith unit

    for (auto k = j; k < lefs.size(); ++k) {
      const auto& rev_idx = rev_lef_rank_buff[k];      // index of the kth rev unit in 5'-3' order
      const auto& rev = lefs[rev_idx].rev_unit.pos();  // pos of the kth unit
      if (rev < fwd) {  // Look for the first rev unit that comes after the current fwd unit
        ++j;  // NOTE k is increased by the for loop. Increasing j does not affect the current for
              // loop
        continue;
      }

      // Note: this delta between uints is safe to do, because at this point rev_pos >= fwd_pos
      const auto delta = rev - fwd;

      // Increment the # of collisions if delta is less than the threshold, and if the two ext.
      // units do not belong to the same LEF
      fwd_collision_buff[fwd_idx] += delta <= this->bin_size && rev_idx != fwd_idx;
      rev_collision_buff[rev_idx] += delta <= this->bin_size && rev_idx != fwd_idx;

      /*
      fmt::print(stderr,
                 "i={}; j={}; fwd_idx={}; rev_idx={}; rev={}; fwd={}; delta={}; fwd_mask={}; "
                 "rev_mask={};\n",
                 i, j, fwd_idx, rev_idx, rev, fwd, delta, fwd_collision_buff[fwd_idx],
                 rev_collision_buff[rev_idx]);
      */
      // Break out of the loop if the units being processed are too far from each other, or if k
      // points to the last rev unit
      if (delta > this->bin_size || k == lefs.size() - 1) {
        // We always move to the next fwd unit, because when this branch is executed, we have
        // already detected all possible collisions for the current fwd unit
        ++i;

        // If we already know that rev == fwd, we can spare one loop iteration by increasing j
        // directly here
        j += rev == fwd;
        break;
      }
    }
  }
}

void Simulation::apply_lef_lef_stalls(absl::Span<Lef> lefs,
                                      absl::Span<const collision_t> rev_collision_buff,
                                      absl::Span<const collision_t> fwd_collision_buff,
                                      absl::Span<const std::size_t> rev_lef_rank_buff,
                                      absl::Span<const std::size_t> fwd_lef_rank_buff,
                                      PRNG& rand_eng) {
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());
  assert(lefs.size() == fwd_collision_buff.size());
  assert(lefs.size() == rev_collision_buff.size());
  assert(std::is_sorted(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                        }));
  assert(std::is_sorted(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                        }));
  assert(std::accumulate(fwd_collision_buff.begin(), fwd_collision_buff.end(), 0UL) ==
         std::accumulate(rev_collision_buff.begin(), rev_collision_buff.end(), 0UL));

  lef_lef_stall_generator_t rng(this->probability_of_extrusion_unit_bypass);
  std::size_t i = 0, j = 0;
  auto fwd_idx = fwd_lef_rank_buff[i];
  auto rev_idx = rev_lef_rank_buff[j];

  std::size_t collisions_fwd = fwd_collision_buff[fwd_idx];
  std::size_t collisions_rev = rev_collision_buff[rev_idx];

  /* Loop over the two buffers storing the number of collision for each LEF. We use the fwd/rev
   * ranks to iterate through collision events in 5'-3' direction. Draw and apply the appropriate
   * number of stalls
   */
  while (i < lefs.size() && j < lefs.size()) {
    while (collisions_rev == 0) {  // Find the first collision event for rev units
      if (++j == lefs.size()) {
        return;
      }
      rev_idx = rev_lef_rank_buff[j];
      collisions_rev = rev_collision_buff[rev_idx];
    }

    while (collisions_fwd == 0) {  // Find the first collision event for fwd units
      if (++i == lefs.size()) {
        return;
      }
      fwd_idx = fwd_lef_rank_buff[i];
      collisions_fwd = fwd_collision_buff[fwd_idx];
    }

    // Each iteration consumes a pair of collision events. When fwd or rev collision events (or
    // both) reach 0, this while loop is terminated, and we go back to the outer while loop (where
    // we look for the next non-zero collision event)
    while (collisions_fwd > 0 && collisions_rev > 0) {
      const auto nstalls =
          std::min({rng(rand_eng), lefs[rev_idx].lifetime, lefs[fwd_idx].lifetime});
      lefs[fwd_idx].fwd_unit._nstalls_lef_lef += nstalls;
      lefs[rev_idx].rev_unit._nstalls_lef_lef += nstalls;
      --collisions_fwd;
      --collisions_rev;
    }
  }
}

void Simulation::check_lef_bar_collisions(absl::Span<const Lef> lefs,
                                          absl::Span<const std::size_t> rev_lef_rank_buff,
                                          absl::Span<const std::size_t> fwd_lef_rank_buff,
                                          absl::Span<const ExtrusionBarrier> extr_barriers,
                                          absl::Span<collision_t> rev_collision_buff,
                                          absl::Span<collision_t> fwd_collision_buff) {
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());
  assert(lefs.size() == fwd_collision_buff.size());
  assert(lefs.size() == rev_collision_buff.size());
  assert(std::is_sorted(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                        }));
  assert(std::is_sorted(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                        }));

  // Clear buffers (-1 indicates no collisions)
  std::fill(fwd_collision_buff.begin(), fwd_collision_buff.end(),
            std::numeric_limits<std::size_t>::max());
  std::fill(rev_collision_buff.begin(), rev_collision_buff.end(),
            std::numeric_limits<std::size_t>::max());

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
  std::size_t j = 0;  // Process fwd extrusion units
  for (auto i = 0UL; i < lefs.size(); ++i) {
    const auto& fwd_idx = fwd_lef_rank_buff[i];  // index of the ith fwd unit in 5'-3' order
    if (lefs[fwd_idx].fwd_unit.stalled()) {
      continue;
    }
    const auto& fwd_pos = lefs[fwd_idx].fwd_unit.pos();  // pos of the ith fwd unit
    auto extr_barr_pos = extr_barriers[j].pos();  // pos of the next (possibly blocking) barrier

    // Look for the first extr. barrier downstream of the current fwd extrusion unit
    while (extr_barr_pos < fwd_pos) {
      if (++j == extr_barriers.size()) {
        goto endloop1;  // Reached the last barrier at the 5'-3' direction
      }
      extr_barr_pos = extr_barriers[j].pos();
    }

    // Detect collisions between extr. units/barriers that are less then delta bp away from each
    // other
    const auto delta = extr_barr_pos - fwd_pos;
    fwd_collision_buff[fwd_idx] = delta <= this->bin_size ? j : -1UL;
  }
endloop1:

  j = extr_barriers.size() - 1;  // Process rev extrusion units
  for (auto i = lefs.size() - 1; i > 0; --i) {
    const auto& rev_idx = rev_lef_rank_buff[i];  // index of the ith fwd unit in 5'-3' order
    // fmt::print(stderr, "i={}/{}; rev_idx={}/{}\n", i, lefs.size(), rev_idx,
    // rev_lef_rank_buff.size());
    if (lefs[rev_idx].rev_unit.stalled()) {
      continue;
    }
    const auto& rev_pos = lefs[rev_idx].rev_unit.pos();  // pos of the ith rev unit
    auto extr_barr_pos = extr_barriers[j].pos();  // pos of the next (possibly blocking) barrier

    // Look for the first extr. barrier upstream of the current rev extrusion unit
    while (extr_barr_pos > rev_pos) {
      if (j-- == 0) {
        return;  // reached the last barrier in 3'-5' direction
      }
      extr_barr_pos = extr_barriers[j].pos();
    }

    // Detect collisions between extr. units/barriers that are less then delta bp away from each
    // other
    const auto delta = rev_pos - extr_barr_pos;
    rev_collision_buff[rev_idx] = delta <= this->bin_size ? j : -1UL;
  }
}

void Simulation::apply_lef_bar_stalls(absl::Span<Lef> lefs,
                                      absl::Span<const collision_t> rev_collision_buff,
                                      absl::Span<const collision_t> fwd_collision_buff,
                                      absl::Span<const ExtrusionBarrier> extr_barriers,
                                      PRNG& rand_eng) {
  assert(lefs.size() == fwd_collision_buff.size());
  assert(lefs.size() == rev_collision_buff.size());

  for (auto i = 0UL; i < lefs.size(); ++i) {
    // Skip over LEFs where no collision occurred or where both units are currently stalled
    if ((rev_collision_buff[i] == NO_COLLISION && fwd_collision_buff[i] == NO_COLLISION) ||
        (lefs[i].rev_unit.stalled() && lefs[i].fwd_unit.stalled())) {
      continue;
    }
    auto& lef = lefs[i];

    // Process rev collisions
    if (const auto barr_idx = rev_collision_buff[i]; barr_idx != NO_COLLISION) {
      assert(!lef.rev_unit.stalled());
      const auto& barrier = extr_barriers[barr_idx];
      auto nstalls = lef_bar_stall_generator_t{1.0 - barrier.pblock()}(rand_eng);
      if (barrier.blocking_direction() != dna::rev) {  // Reduce nstalls in case of soft-stalls
        assert(barrier.blocking_direction() == dna::fwd);
        nstalls = static_cast<Bp>(std::round(soft_stall_multiplier * static_cast<double>(nstalls)));
      }
      lef.rev_unit._nstalls_lef_bar = nstalls;
    }

    // Process fwd collisions
    if (const auto barr_idx = fwd_collision_buff[i]; barr_idx != NO_COLLISION) {
      assert(!lef.fwd_unit.stalled());
      const auto& barrier = extr_barriers[barr_idx];
      auto nstalls = lef_bar_stall_generator_t{1.0 - barrier.pblock()}(rand_eng);
      if (barrier.blocking_direction() != dna::fwd) {  // Reduce nstalls in case of soft-stalls
        assert(barrier.blocking_direction() == dna::rev);
        nstalls = static_cast<Bp>(std::round(soft_stall_multiplier * static_cast<double>(nstalls)));
      }
      lef.fwd_unit._nstalls_lef_bar = nstalls;
    }

    const auto rev_barr_idx = rev_collision_buff[i];
    const auto fwd_barr_idx = fwd_collision_buff[i];
    // Detect and process hard-stalls
    if (lef.rev_unit.stalled() && lef.fwd_unit.stalled() && rev_barr_idx != NO_COLLISION &&
        fwd_barr_idx != NO_COLLISION &&
        extr_barriers[rev_barr_idx].blocking_direction() == dna::rev &&
        extr_barriers[fwd_barr_idx].blocking_direction() == dna::fwd) {
      const auto nstalls = static_cast<Bp>(std::round(
          hard_stall_multiplier * static_cast<double>(std::min(lef.rev_unit._nstalls_lef_bar,
                                                               lef.fwd_unit._nstalls_lef_bar))));
      lef.rev_unit._nstalls_lef_bar += nstalls;
      lef.fwd_unit._nstalls_lef_bar += nstalls;

      if (this->allow_lef_lifetime_extension) {
        assert(std::numeric_limits<Bp>::max() - nstalls > lef.lifetime);  // NOLINT
        lef.lifetime += nstalls;
      }
    }
  }
}

void Simulation::register_contacts(Chromosome* chrom, absl::Span<const Lef> lefs) {
  for (const auto& lef : lefs) {  // Register contacts at LEF sites (excluding LEFs that have one
                                  // of their units at the beginning/end of a chromosome)
    if (lef.is_bound() && lef.rev_unit.pos() != chrom->start_pos() &&
        lef.rev_unit.pos() != chrom->end_pos() - 1 && lef.fwd_unit.pos() != chrom->start_pos() &&
        lef.fwd_unit.pos() != chrom->end_pos() - 1) {
      chrom->increment_contacts(lef.rev_unit.pos(), lef.fwd_unit.pos(), this->bin_size);
    }
  }
}

template <typename MaskT>
void Simulation::select_lefs_to_bind(absl::Span<const Lef> lefs, MaskT& mask) {
  static_assert(std::is_same_v<boost::dynamic_bitset<>, MaskT> ||
                    std::is_integral_v<std::decay_t<decltype(std::declval<MaskT&>().operator[](
                        std::declval<std::size_t>()))>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");
  assert(lefs.size() <= mask.size());
  for (auto i = 0UL; i < lefs.size(); ++i) {
    mask[i] = lefs[i].lifetime == 0;
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

}  // namespace modle
