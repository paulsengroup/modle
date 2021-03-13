#pragma once

#include <absl/container/flat_hash_map.h>
#include <absl/strings/str_cat.h>
#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration, Time
#include <absl/types/span.h>
#include <fmt/format.h>  // for FMT_STRING
#ifdef USE_XOSHIRO
#include <XoshiroCpp.hpp>
#endif

#include <cpp-sort/sorters/drop_merge_sorter.h>
#include <cpp-sort/sorters/insertion_sorter.h>

#include <algorithm>  // for max, min, sort, transform, all_of, for_each, remove_if
#include <atomic>     // for atomic, memory_order_relaxed
#include <boost/asio/thread_pool.hpp>
#include <boost/dynamic_bitset.hpp>
#include <cassert>
#include <chrono>  // for seconds
#include <cmath>   // for llround, floor, lround, sqrt
#include <condition_variable>
#include <cstdint>     // for uint*_t, UINT*_MAX
#include <cstdio>      // for stderr
#include <filesystem>  // for create_directories
#include <functional>  // for greater, hash
#include <iosfwd>      // for size_t
#include <mutex>       // for mutex, unique_lock, scoped_lock
#include <numeric>     // for accumulate, partial_sum
#include <random>  // for mt19937,  bernoulli_distribution, seed_seq, discrete_distribution, uniform_int_distribution
#include <stdexcept>  // for runtime_error
#include <thread>
#include <type_traits>  // for declval

#include "modle/bed.hpp"        // for BED, Parser, BED::BED6, BED::Standard
#include "modle/chr_sizes.hpp"  // for ChrSize, Parser
#include "modle/common.hpp"
#include "modle/config.hpp"    // for config
#include "modle/contacts.hpp"  // for ContactMatrix
#include "modle/cooler.hpp"
#include "modle/extrusion_barriers.hpp"
#include "modle/extrusion_factors.hpp"
#include "modle/suppress_compiler_warnings.hpp"

namespace modle {

#ifdef USE_XOSHIRO
using PRNG = XoshiroCpp::Xoshiro256PlusPlus;
using seeder = XoshiroCpp::SplitMix64;
#else
using PRNG = std::mt19937_64;
using seeder = std::seed_seq;
#endif

void Genome::test() {
  for (const auto& chrom : this->_chromosomes) {
    fmt::print(stderr, "{}: {}\n", chrom.name(), chrom.nbarriers());

    for (const auto& b : chrom.get_barriers()) {
      fmt::print(stderr, "\t{}[{}-{}]", b.chrom, b.chrom_start, b.chrom_end);
    }
    fmt::print(stderr, "\n");
    if (chrom.nbarriers() > 0) {
      this->simulate_extrusion(1);
    }
  }
}

Genome::Genome(const config& c, bool import_chroms)
    : _path_to_chrom_sizes(c.path_to_chr_sizes),
      _path_to_chrom_subranges(c.path_to_chr_subranges),
      _path_to_extr_barriers(c.path_to_extr_barriers_bed),
      _bin_size(c.bin_size),
      _diagonal_width(c.diagonal_width),
      _avg_lef_lifetime(c.average_lef_lifetime),
      _nlefs_per_mbp(c.number_of_lefs_per_mbp),
      _probability_of_barrier_block(c.probability_of_extrusion_barrier_block),
      _probability_of_lef_rebind(c.probability_of_lef_rebind),
      _probability_of_extr_unit_bypass(c.probability_of_extrusion_unit_bypass),
      _soft_stall_multiplier(c.soft_stall_multiplier),
      _hard_stall_multiplier(c.hard_stall_multiplier),
      _allow_lef_lifetime_extension(c.allow_lef_lifetime_extension),
      _sampling_interval(c.contact_sampling_interval),
      _randomize_contact_sampling(c.randomize_contact_sampling_interval),
      _nthreads(std::min(std::thread::hardware_concurrency(), c.nthreads)),
      _seed(c.seed),
      _chromosomes(import_chroms ? import_chromosomes(_path_to_chrom_sizes, _path_to_extr_barriers,
                                                      _path_to_chrom_subranges)
                                 : Chromosomes{}) {}

std::size_t Genome::size() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
      [](auto accumulator, const auto& chrom) { return accumulator + chrom.size(); });
}

std::size_t Genome::simulated_size() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
                         [](auto accumulator, const auto& chrom) {
                           return accumulator + (chrom.end_pos() - chrom.start_pos());
                         });
}

// TODO Add flag to import skip chrom without barriers
Genome::Chromosomes Genome::import_chromosomes(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_extr_barriers,
    const std::filesystem::path& path_to_chrom_subranges) {
  assert(!path_to_chrom_sizes.empty());
  assert(!path_to_extr_barriers.empty());
  Chromosomes chroms;

  {
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
      } else if (chrom_ranges.empty()) {
        chroms.emplace(std::move(chrom));
      }
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
    if (auto match = chroms.find(record.chrom); match != chroms.end()) {
      fmt::print(stderr, "Adding barrier for {}[{}-{}]...\n", record.chrom, record.chrom_start,
                 record.chrom_end);
      match->add_extrusion_barrier(record);
    }
  }
  return chroms;
}

/*
void Genome::simulate_extrusion() { this->simulate_extrusion(1, 0); }
void Genome::simulate_extrusion(uint32_t iterations) { this->simulate_extrusion(iterations, 0); }
void Genome::simulate_extrusion(double target_contact_density) {
  this->simulate_extrusion(0, target_contact_density);
}

void Genome::simulate_extrusion(uint32_t iterations, double target_contact_density) {
  auto get_status = [](const Chromosome& c) constexpr->std::string_view {
    if (c.ok()) {
      return "OK!";
    }
    if (c.barriers.empty()) {
      return "KO! Chromosome won't be simulated. Reason: chromosome has 0 extrusion barriers.";
    }
    utils::throw_with_trace(std::logic_error("Unreachable code"));
  };

  fmt::print(stderr, FMT_STRING("Chromosome status report:\n"));

  for (const auto& chr : this->_chromosomes) {
    fmt::print(stderr, FMT_STRING("'{}' status: {}\n"), chr.name, get_status(chr));
  }
  fmt::print(stderr, FMT_STRING("Simulating loop extrusion on {}/{} chromosomes...\n"),
             this->get_n_ok_chromosomes(), this->get_nchromosomes());

  // If the simulation is set to stop when a target contact density is reached, set the number of
  // iterations to a very large number (2^32)
  if (target_contact_density != 0.0) {
    iterations = UINT32_MAX;
  }

  auto tpool = this->instantiate_thread_pool();

  // Initialize variables for simulation progress tracking
  std::atomic<uint64_t> ticks_done{0};
  std::atomic<uint64_t> extrusion_events{0};
  std::atomic<uint64_t> chromosomes_completed{0};
  std::atomic<bool> simulation_completed{false};
  std::mutex m;  // This mutex is supposed to protect simulation_complete. As of C++17 we also
need
  // to acquire a lock when modifying the condition variable simulation_completed_cv
  // (otherwise the change might not be communicated to waiting threads)
  std::condition_variable simulation_completed_cv;

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_CONVERSION
  DISABLE_WARNING_DOUBLE_PROMOTION
  std::thread progress_tracker([&]() {  // This thread is used to periodically print the
    // simulation progress to stderr
    // The total number of ticks tot_ticks is set depending on whether or not the simulation is
    // set to run for a fixed number of iterations:
    //  - When we know in advance the number of iterations (e.g. because it was specified
    //    through --number-of-iterations), then the total number of ticks is calculated as:
    //       n. of iterations * n. of chromosomes
    //  - When the target number of iterations is unknown (e.g. because we are targeting the
    //    contact density specified through --target-contact-density), then the total ticks
    //    number is given by the sum of the target number of contacts for each of the
    //    chromosomes simulated, where the target number of contacts is defined as:
    //       target contact density * (matrix columns * matrix rows)
    const long double tot_ticks =
        target_contact_density != 0.0
            ? std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0.0L,
                              [&](auto accumulator, const auto& chr) {
                                return accumulator + (target_contact_density *
                                                      chr.contacts.ncols() *
chr.contacts.nrows());
                              })
            : iterations * this->get_n_ok_chromosomes();

    const auto t0 = absl::Now();
    while (!simulation_completed) {
      {  // Wait for 5 seconds or until the simulation terminates
        std::unique_lock<std::mutex> lk(m);
        simulation_completed_cv.wait_for(lk, std::chrono::seconds(5));  // NOLINT
        if (simulation_completed) {  // Stop printing updates once there are no chromosomes left,
          // and return immediately, so that the main thread can join this thread
          return;
        }
      }
      if (extrusion_events > 0) {  // Print progress
        const auto progress = ticks_done / tot_ticks;
        const auto throughput = static_cast<double>(extrusion_events) / 5.0e6; // 5s * 1M NOLINT
        const auto delta_t = absl::ToDoubleSeconds(absl::Now() - t0);
        const auto eta = absl::Seconds(
            (std::min(this->_nthreads, this->get_n_ok_chromosomes()) /
             std::min(static_cast<double>(this->_nthreads),
                      static_cast<double>(this->get_n_ok_chromosomes() - chromosomes_completed)))
* (delta_t / std::max(1.0e-06L, progress) - delta_t)); if (eta > absl::ZeroDuration()) {  // Avoid
printing updates when simulation is about to end fmt::print(stderr, FMT_STRING("### ~{:.2f}% ###
{:.2f}M extr/sec - Simulation completed for "
                                "{}/{} chromosomes - ETA {}.\n"),
                     100.0 * progress, throughput, chromosomes_completed,
                     this->get_n_ok_chromosomes(), absl::FormatDuration(eta));
          extrusion_events = 0;
        }
      }
    }
  });
  DISABLE_WARNING_POP

  // Loop extrusion is simulated using boost::thread_pool, where each thread simulates loop
  // extrusion at the chromosome level. This constrains the level of parallelism to the number of
  // chromosomes that are being simulated. We can certainly do better.
  // This is just a way to get significant speedups with very little effort
  for (auto nchr = 0U; nchr < this->_chromosomes.size(); ++nchr) {
    // Simulating loop extrusion on chromosomes without extr. barrier does not make sense
    if (!this->_chromosomes[nchr].ok()) {
      continue;
    }
    boost::asio::post(tpool, [&, nchr]() {
      const auto t0 = absl::Now();
      auto& chr = this->_chromosomes[nchr];  // Alias for the chromosome that is being simulated
      // This random number generator is used to randomize sampling interval
      // (e.g. when --randomize-contact-sampling-interval is specified)
      std::bernoulli_distribution sample_contacts{1.0 / this->_sampling_interval};

      // Calculate the number of contacts after which the simulation for the current chromosome is
      // stopped. This is set to a very large number (2^64) when the simulation is set to run for
      // a fixed number of iterations
      const auto target_n_of_contacts =
          target_contact_density != 0.0
              ? target_contact_density *
                    static_cast<double>(chr.contacts.nrows() * chr.contacts.ncols())
              : std::numeric_limits<double>::max();

      // Variables to track the simulation progress for the current chromosome
      uint64_t local_extr_events_counter = 0;
      uint64_t ticks_local = 0;

      // This for loop is where the simulation actually takes place.
      // We basically keep iterating until one of the stopping condition is met
      // TODO: Instead of a raw loop, we could use std::transform and exec. policies and achieve a
      // better level of parallelism To do this, we probably need to have a lock on lefs and
      // process all left units first, then right all the right ones
      for (auto i = 1UL; i <= iterations; ++i) {
        // Determine whether we will register contacts produced during the current iteration
        bool register_contacts = this->_randomize_contact_sampling
                                     ? sample_contacts(chr._rand_eng)
                                     : i % this->_sampling_interval == 0;

        // Loop over the LEFs belonging to the chromosome that is being simulated and
        // extrude/register contacts when appropriate
        for (auto& lef : chr.lefs) {
          if (lef->is_bound()) {  // Attempt to register a contact and move the LEF only if the
            // latter is bound to the chromosome
            if (register_contacts) {
              lef->register_contact();
            }
            lef->try_extrude();
            ++local_extr_events_counter;
          }
        }

        // Once we are done extruding for the current iteration, check if the constrains are
        // satisfied and apply a stall or rebind free LEFs when appropriate
        for (auto& lef : chr.lefs) {
          if (lef->is_bound()) {
            lef->check_constraints(chr._rand_eng);
          } else {
            lef->try_rebind(chr._rand_eng, this->_probability_of_lef_rebind, register_contacts);
          }
        }

        // Add the number of extr. events to the global counter. This is used to calculate the n.
        // of extr. events per seconds, which is a good way to assess the simulation throughput
        extrusion_events.fetch_add(local_extr_events_counter, std::memory_order_relaxed);
        local_extr_events_counter = 0;

        if (register_contacts) {
          // Propagate local progress to the global counters
          if (target_contact_density !=
              0.0) {  // Simulation is set to stop when the target contact density is reached
            assert(chr.contacts.get_tot_contacts() >= ticks_local);
            ticks_done.fetch_add(chr.contacts.get_tot_contacts() - ticks_local,
                                 std::memory_order_relaxed);
            ticks_local = chr.contacts.get_tot_contacts();
          } else {  // Simulation is set to stop at a fixed number of iterations
            assert(i >= ticks_local);
            ticks_done.fetch_add(i - ticks_local, std::memory_order_relaxed);
            ticks_local = i;
          }

          // Alt the simulation when we have reached the target number of contacts.
          // This will never happen when the simulation is set to run for a fixed number of
          // iterations, as when this is the case, target_n_of_contacts == 2^64
          if (static_cast<double>(chr.contacts.get_tot_contacts()) >= target_n_of_contacts) {
            break;
          }
        }
      }
      if (++chromosomes_completed == this->get_n_ok_chromosomes()) {
        {  // Notify the thread that is tracking the simulation progress that we are done
          // simulating
          std::scoped_lock<std::mutex> lk(m);
          simulation_completed = true;
        }
        simulation_completed_cv.notify_all();
      }
      fmt::print(stderr, FMT_STRING("DONE simulating loop extrusion on '{}'! Simulation took
{}\n"), chr.name, absl::FormatDuration(absl::Now() - t0));

      // Free up resources
      // chr = Chromosome{"", 0, 0, 0};
    });
  }
  // This blocks until all the tasks posted to the thread_pool have been completed
  tpool.join();
  progress_tracker.join();
}
*/

void Genome::simulate_extrusion(std::size_t nthreads) {
  std::size_t ncells = 5000;
  std::size_t simulation_rounds = 200;
  auto tpool = this->instantiate_thread_pool(nthreads);

  std::vector<ExtrusionBarrier> extr_barriers_buff;
  std::vector<Lef> lef_buff;

  std::mutex cooler_mutex;
  cooler::Cooler c(std::filesystem::path("/tmp/matrix.cool"), cooler::Cooler::WRITE_ONLY,
                   this->_bin_size);
  std::mutex progress_mutex;
  absl::flat_hash_map<const Chromosome* const, std::size_t> progress;
  std::transform(this->_chromosomes.begin(), this->_chromosomes.end(),
                 std::inserter(progress, progress.begin()),
                 [](const auto& chr) { return std::make_pair(&chr, 0UL); });

  for (auto& chrom : this->_chromosomes) {
    if (!chrom.ok()) {
      progress.erase(&chrom);
      continue;
    }
    chrom.allocate_contacts(this->_bin_size, this->_diagonal_width);
    const auto nlefs = static_cast<std::size_t>(
        std::round(this->_nlefs_per_mbp * (static_cast<double>(chrom.simulated_size()) / 1.0e6)));
    // Clearing buffers
    extr_barriers_buff.clear();
    lef_buff.resize(nlefs);

    for (const auto& b : chrom.get_barriers()) {
      if (b.strand == '+' || b.strand == '-') {
        const auto pos =
            static_cast<Bp>(std::round(static_cast<double>(b.chrom_start + b.chrom_end) / 2.0));
        const auto pblock = b.score != 0 ? b.score : this->_probability_of_barrier_block;
        extr_barriers_buff.emplace_back(pos, pblock, b.strand);
      }
    }

    for (auto i = 0UL; i < ncells; ++i) {
      boost::asio::post(tpool, [&, i]() {
        try {
          Genome::simulate_extrusion_kernel(&chrom, i, simulation_rounds, lef_buff,
                                            extr_barriers_buff);
          fmt::print(stderr, "Done simulating cell #{}!\n", i + 1);
          bool write_matrix = false;
          {
            std::scoped_lock l(progress_mutex);
            progress[&chrom] += 1;
            if (progress[&chrom] == ncells) {
              write_matrix = true;
            }
          }
          if (write_matrix) {
            std::scoped_lock l(cooler_mutex);
            c.write_or_append_cmatrix_to_file(chrom.contacts(), chrom.name(), chrom.start_pos(),
                                              chrom.end_pos(), chrom.size());
          }

        } catch (const std::exception& err) {
          fmt::print(stderr, FMT_STRING("{}\n"), err.what());
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
  }
  tpool.join();
  c.write_metadata();
}

void Genome::simulate_extrusion_kernel(Chromosome* chrom, std::size_t cell_id,
                                       std::size_t simulation_rounds, std::vector<Lef> lef_buff,
                                       std::vector<ExtrusionBarrier> extr_barrier_buff) {
  // Bind all lefs
  const auto seed = this->_seed + std::hash<std::string_view>{}(chrom->name()) +
                    std::hash<std::size_t>{}(chrom->size()) + std::hash<std::size_t>{}(cell_id);
#ifdef USE_XOSHIRO
  modle::seeder seeder_(seed);
  modle::PRNG rand_eng(seeder_.generateSeedSequence<4>());
#else
  modle::seeder seeder_{seed};
  modle::PRNG rand_eng(seeder_);
#endif

  std::vector<std::size_t> lef_initial_loading_round(lef_buff.size());
  {
    absl::btree_multiset<std::size_t> sorted_rounds;
    std::uniform_int_distribution<std::size_t> round_gen{
        0, (4 * this->_avg_lef_lifetime) / this->_bin_size};

    std::generate_n(std::inserter(sorted_rounds, sorted_rounds.begin()), lef_buff.size(),
                    [&]() { return round_gen(rand_eng); });
    std::transform(sorted_rounds.begin(), sorted_rounds.end(), lef_initial_loading_round.rbegin(),
                   [&](const auto n) { return n - *sorted_rounds.begin(); });
  }
  const auto burnin_rounds = lef_initial_loading_round.front();
  simulation_rounds += burnin_rounds;

  std::vector<std::size_t> rev_lef_rank_buff(lef_buff.size());
  std::iota(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), 0);
  std::vector<std::size_t> fwd_lef_rank_buff{rev_lef_rank_buff};
  boost::dynamic_bitset<> mask{lef_buff.size()};

  // this->bind_all_lefs(chrom, lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rand_eng);

  std::vector<collision_t> rev_lef_collision_buff(lef_buff.size());
  std::vector<collision_t> fwd_lef_collision_buff(lef_buff.size());
  std::bernoulli_distribution sample_contacts{1.0 / this->_sampling_interval};

  auto lefs = absl::MakeSpan(lef_buff);
  const auto barriers = absl::MakeConstSpan(extr_barrier_buff);
  auto rev_lef_ranks = absl::MakeSpan(rev_lef_rank_buff);
  auto fwd_lef_ranks = absl::MakeSpan(fwd_lef_rank_buff);
  auto rev_lef_collisions = absl::MakeSpan(rev_lef_collision_buff);
  auto fwd_lef_collisions = absl::MakeSpan(fwd_lef_collision_buff);

  for (auto round = 0UL; round < simulation_rounds; ++round) {
    const auto register_contacts =
        round > burnin_rounds &&
        (this->_randomize_contact_sampling ? sample_contacts(rand_eng)
                                           : round % this->_sampling_interval == 0);

    if (!lef_initial_loading_round.empty()) {
      while (!lef_initial_loading_round.empty() && lef_initial_loading_round.back() == round) {
        lef_initial_loading_round.pop_back();
      }
      const auto nlefs = lef_buff.size() - lef_initial_loading_round.size();

      lefs = absl::MakeSpan(lef_buff.data(), nlefs);
      rev_lef_ranks = absl::MakeSpan(rev_lef_rank_buff.data(), nlefs);
      fwd_lef_ranks = absl::MakeSpan(fwd_lef_rank_buff.data(), nlefs);
      rev_lef_collisions = absl::MakeSpan(rev_lef_collision_buff.data(), nlefs);
      fwd_lef_collisions = absl::MakeSpan(fwd_lef_collision_buff.data(), nlefs);
    }

    this->select_lefs_to_bind(lefs, mask);
    this->bind_lefs(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rand_eng, mask);

    this->check_lef_lef_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_lef_collisions,
                                   fwd_lef_collisions);
    this->apply_lef_lef_stalls(lefs, rev_lef_collisions, fwd_lef_collisions, rev_lef_ranks,
                               fwd_lef_ranks, rand_eng, this->_probability_of_extr_unit_bypass);

    this->check_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, barriers, rev_lef_collisions,
                                   fwd_lef_collisions);
    this->apply_lef_bar_stalls(lefs, rev_lef_collisions, fwd_lef_collisions, barriers,
                               rev_lef_ranks, fwd_lef_ranks, rand_eng, this->_soft_stall_multiplier,
                               this->_hard_stall_multiplier);

    if (register_contacts) {
      this->register_contacts(chrom, lefs);
    }

    this->extrude(chrom, lefs);
    this->rank_lefs(lefs, fwd_lef_ranks, rev_lef_ranks);
  }
}

template <typename MaskT>
void Genome::bind_lefs(const Chromosome* const chrom, absl::Span<Lef> lefs,
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
  lef_lifetime_generator_t lifetime_generator{1.0 / static_cast<double>(this->_avg_lef_lifetime)};
  for (auto i = 0UL; i < lefs.size(); ++i) {
    // Here we only consider even indexes, as binding a single extr. unit is not
    // allowed: either we bind all extr. units belonging to a LEF, or we bind none
    if (mask.empty() || mask[i]) {
      lefs[i].bind_at_pos(pos_generator(rand_eng), lifetime_generator(rand_eng));
    }
  }
  this->rank_lefs(lefs, rev_lef_rank_buff, fwd_lef_rank_buff);
}

void Genome::bind_all_lefs(const Chromosome* chrom, absl::Span<Lef> lefs,
                           absl::Span<std::size_t> rev_lef_rank_buff,
                           absl::Span<std::size_t> fwd_lef_rank_buff, PRNG& rand_eng) {
  chrom_pos_generator_t chrom_pos_generator{chrom->start_pos(), chrom->end_pos()};
  lef_lifetime_generator_t lef_lifetime_generator{1.0 /
                                                  static_cast<double>(this->_avg_lef_lifetime)};
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

void Genome::rank_lefs(absl::Span<Lef> lefs, absl::Span<std::size_t> rev_lef_rank_buff,
                       absl::Span<std::size_t> fwd_lef_rank_buff, bool init_buffers) {
  if (init_buffers) {
    std::iota(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), 0);
    std::copy(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), rev_lef_rank_buff.begin());
  }
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());

#if 0  // 55.77s
  std::sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
  });

  std::sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
  });
#else  // 35.50s
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
    cppsort::drop_merge_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                             [&](const auto r1, const auto r2) {
                               return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                             });

    cppsort::drop_merge_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                             [&](const auto r1, const auto r2) {
                               return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                             });
  }
#endif
}

void Genome::extrude(const Chromosome* chrom, absl::Span<Lef> lefs) {
  for (auto& lef : lefs) {
    assert(lef.lifetime > 0);
    if (lef.rev_unit.stalled()) {
      lef.rev_unit.decrement_stalls(this->_bin_size);
    } else if (lef.rev_unit._pos < chrom->start_pos() + this->_bin_size) {
      lef.rev_unit._pos = chrom->start_pos();
      lef.rev_unit._nstalls_lef_lef = std::numeric_limits<Bp>::max();
    } else {
      lef.rev_unit._pos -= this->_bin_size;
    }

    if (lef.fwd_unit.stalled()) {
      lef.fwd_unit.decrement_stalls(this->_bin_size);
    } else if (lef.fwd_unit._pos + this->_bin_size > chrom->end_pos() - 1) {
      lef.fwd_unit._pos = chrom->end_pos() - 1;
      lef.fwd_unit._nstalls_lef_lef = std::numeric_limits<Bp>::max();
    } else {
      lef.fwd_unit._pos += this->_bin_size;
    }

    if ((lef.lifetime -= std::min(lef.lifetime, 2 * this->_bin_size)) == 0) {
      lef.release();
    }
  }
}

void Genome::check_lef_lef_collisions(absl::Span<const Lef> lefs,
                                      absl::Span<const std::size_t> rev_lef_rank_buff,
                                      absl::Span<const std::size_t> fwd_lef_rank_buff,
                                      absl::Span<collision_t> rev_collision_buff,
                                      absl::Span<collision_t> fwd_collision_buff,
                                      Bp dist_threshold) {
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

      // Note: this delta between uint is safe to do, because at this point rev_pos >= fwd_pos
      const auto delta = rev - fwd;

      // Increment the # of collisions if delta is less than the threshold, and if the two ext.
      // units do not belong to the same LEF
      fwd_collision_buff[fwd_idx] += delta <= dist_threshold && rev_idx != fwd_idx;
      rev_collision_buff[rev_idx] += delta <= dist_threshold && rev_idx != fwd_idx;

      /*
      fmt::print(stderr,
                 "i={}; j={}; fwd_idx={}; rev_idx={}; rev={}; fwd={}; delta={}; fwd_mask={}; "
                 "rev_mask={};\n",
                 i, j, fwd_idx, rev_idx, rev, fwd, delta, fwd_collision_buff[fwd_idx],
                 rev_collision_buff[rev_idx]);
      */
      // Break out of the loop if the units being processed are too far from each other, or if k
      // points to the last rev unit
      if (delta > dist_threshold || k == lefs.size() - 1) {
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

void Genome::check_lef_lef_collisions(absl::Span<const Lef> lefs,
                                      absl::Span<const std::size_t> rev_lef_rank_buff,
                                      absl::Span<const std::size_t> fwd_lef_rank_buff,
                                      absl::Span<collision_t> rev_collision_buff,
                                      absl::Span<collision_t> fwd_collision_buff) {
  this->check_lef_lef_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_collision_buff,
                                 fwd_collision_buff, this->_bin_size);
}

void Genome::apply_lef_lef_stalls(absl::Span<Lef> lefs,
                                  absl::Span<const collision_t> rev_collision_buff,
                                  absl::Span<const collision_t> fwd_collision_buff,
                                  absl::Span<const std::size_t> rev_lef_rank_buff,
                                  absl::Span<const std::size_t> fwd_lef_rank_buff, PRNG& rand_eng,
                                  double prob_of_bypass) {
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

  lef_lef_stall_generator_t rng(prob_of_bypass);
  std::size_t i = 0, j = 0;
  auto fwd_idx = fwd_lef_rank_buff[i];
  auto rev_idx = rev_lef_rank_buff[j];

  std::size_t collisions_fwd = fwd_collision_buff[fwd_idx];
  std::size_t collisions_rev = rev_collision_buff[rev_idx];

  /* Loop over the two buffers storing the number of collision for each LEF. We use the fwd/rev
   * ranks to iterate through collision events in 5'-3' direction. Draw and apply the appropriate
   * number of stalls
   */
  while (true) {
    while (collisions_fwd == 0) {  // Find the first collision event for fwd units
      if (++i >= lefs.size()) {
        assert(collisions_rev == 0);
        return;  // All collisions have been processed
      }
      fwd_idx = fwd_lef_rank_buff[i];
      collisions_fwd = fwd_collision_buff[fwd_idx];
      // Reset the number of stalls if the number of collisions for the current extr. unit is 0
      lefs[fwd_idx].fwd_unit._nstalls_lef_lef *= collisions_fwd != 0;
    }

    while (collisions_rev == 0) {  // Find the first collision event for rev units
      if (++j >= lefs.size()) {
        assert(collisions_fwd == 0);
        return;  // All collisions have been processed
      }
      rev_idx = rev_lef_rank_buff[j];
      collisions_rev = rev_collision_buff[rev_idx];
      // Reset the number of stalls if the number of collisions for the current extr. unit is 0
      lefs[rev_idx].rev_unit._nstalls_lef_lef *= collisions_rev != 0;
    }

    // Each iteration consumes a pair of collision events. When fwd or rev collision events (or
    // both) reach 0, this while loop is terminated, and we go back to the outer while loop (where
    // we look for the next non-zero collision event)
    while (collisions_fwd > 0 && collisions_rev > 0) {
      const auto nstalls = rng(rand_eng);
      lefs[fwd_idx].fwd_unit._nstalls_lef_lef += nstalls;
      lefs[rev_idx].rev_unit._nstalls_lef_lef += nstalls;
      --collisions_fwd;
      --collisions_rev;
    }
  }
}

void Genome::check_lef_bar_collisions(absl::Span<const Lef> lefs,
                                      absl::Span<const std::size_t> rev_lef_rank_buff,
                                      absl::Span<const std::size_t> fwd_lef_rank_buff,
                                      absl::Span<const ExtrusionBarrier> extr_barriers,
                                      absl::Span<collision_t> rev_collision_buff,
                                      absl::Span<collision_t> fwd_collision_buff,
                                      Bp dist_threshold) {
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
  std::size_t j = 0;
  for (auto i = 0UL; i < lefs.size(); ++i) {
    const auto& fwd_idx = fwd_lef_rank_buff[i];  // index of the ith fwd unit in 5'-3' order
    if (lefs[fwd_idx].fwd_unit.stalled()) {
      continue;
    }
    const auto& fwd_pos = lefs[fwd_idx].fwd_unit.pos();  // pos of the ith unit
    auto extr_barr_pos = extr_barriers[j].pos();

    while (extr_barr_pos < fwd_pos) {
      if (++j == extr_barriers.size()) {
        goto endloop1;
      }
      extr_barr_pos = extr_barriers[j].pos();
    }

    const auto delta = extr_barr_pos - fwd_pos;
    fwd_collision_buff[fwd_idx] = delta <= dist_threshold ? j : -1UL;
  }
endloop1:

  j = extr_barriers.size() - 1;
  for (auto i = lefs.size() - 1; i > 0; --i) {
    const auto& rev_idx = rev_lef_rank_buff[i];  // index of the ith fwd unit in 5'-3' order
    if (lefs[rev_idx].rev_unit.stalled()) {
      continue;
    }
    const auto& rev_pos = lefs[rev_idx].rev_unit.pos();  // pos of the ith unit
    auto extr_barr_pos = extr_barriers[j].pos();

    while (extr_barr_pos > rev_pos) {
      if (j-- == 0) {
        return;
      }
      extr_barr_pos = extr_barriers[j].pos();
    }

    const auto delta = rev_pos - extr_barr_pos;
    rev_collision_buff[rev_idx] = delta <= dist_threshold ? j : -1UL;
  }
}

void Genome::check_lef_bar_collisions(absl::Span<const Lef> lefs,
                                      absl::Span<const std::size_t> rev_lef_rank_buff,
                                      absl::Span<const std::size_t> fwd_lef_rank_buff,
                                      absl::Span<const ExtrusionBarrier> extr_barriers,
                                      absl::Span<collision_t> rev_collision_buff,
                                      absl::Span<collision_t> fwd_collision_buff) {
  this->check_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, extr_barriers,
                                 rev_collision_buff, fwd_collision_buff, this->_bin_size);
}

void Genome::apply_lef_bar_stalls(absl::Span<Lef> lefs,
                                  absl::Span<const collision_t> rev_collision_buff,
                                  absl::Span<const collision_t> fwd_collision_buff,
                                  absl::Span<const ExtrusionBarrier> extr_barriers,
                                  absl::Span<const std::size_t> rev_lef_rank_buff,
                                  absl::Span<const std::size_t> fwd_lef_rank_buff, PRNG& rand_eng,
                                  double soft_stall_multiplier, double hard_stall_multiplier) {
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

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (rev_collision_buff[i] == -1 && fwd_collision_buff[i] == -1 ||
        (lefs[i].rev_unit.stalled() && lefs[i].fwd_unit.stalled())) {
      continue;
    }
    auto& lef = lefs[i];

    if (const auto barr_idx = rev_collision_buff[i]; barr_idx != -1) {
      assert(!lef.rev_unit.stalled());
      const auto& barrier = extr_barriers[barr_idx];
      auto nstalls = lef_bar_stall_generator_t{1.0 - barrier.pblock()}(rand_eng);
      if (barrier.blocking_direction() != dna::rev) {
        assert(barrier.blocking_direction() == dna::fwd);
        nstalls = static_cast<Bp>(std::round(soft_stall_multiplier * static_cast<double>(nstalls)));
      }
      lef.rev_unit._nstalls_lef_bar = nstalls;
    }

    if (const auto barr_idx = fwd_collision_buff[i]; barr_idx != -1) {
      assert(!lef.fwd_unit.stalled());
      const auto& barrier = extr_barriers[barr_idx];
      auto nstalls = lef_bar_stall_generator_t{1.0 - barrier.pblock()}(rand_eng);
      if (barrier.blocking_direction() != dna::fwd) {
        assert(barrier.blocking_direction() == dna::rev);
        nstalls = static_cast<Bp>(std::round(soft_stall_multiplier * static_cast<double>(nstalls)));
      }
      lef.fwd_unit._nstalls_lef_bar = nstalls;
    }

    const auto rev_barr_idx = rev_collision_buff[i];
    const auto fwd_barr_idx = fwd_collision_buff[i];

    if (lef.rev_unit.stalled() && lef.fwd_unit.stalled() && rev_barr_idx != -1 &&
        fwd_barr_idx != -1 && extr_barriers[rev_barr_idx].blocking_direction() == dna::rev &&
        extr_barriers[fwd_barr_idx].blocking_direction() == dna::fwd) {
      const auto nstalls = static_cast<Bp>(std::round(
          hard_stall_multiplier * static_cast<double>(std::min(lef.rev_unit._nstalls_lef_bar,
                                                               lef.fwd_unit._nstalls_lef_bar))));
      lef.rev_unit._nstalls_lef_bar += nstalls;
      lef.fwd_unit._nstalls_lef_bar += nstalls;

      if (this->_allow_lef_lifetime_extension) {
        assert(std::numeric_limits<Bp>::max() - nstalls > lef.lifetime);  // NOLINT
        lef.lifetime += nstalls;
      }
    }
  }
}

void Genome::register_contacts(Chromosome* chrom, absl::Span<const Lef> lefs) {
  for (const auto& lef : lefs) {
    if (lef.is_bound() && lef.rev_unit.pos() != chrom->start_pos() &&
        lef.rev_unit.pos() != chrom->end_pos() - 1 && lef.fwd_unit.pos() != chrom->start_pos() &&
        lef.fwd_unit.pos() != chrom->end_pos() - 1) {
      chrom->increment_contacts(lef.rev_unit.pos(), lef.fwd_unit.pos(), this->_bin_size);
    }
  }
}

template <typename MaskT>
void Genome::select_lefs_to_bind(absl::Span<const Lef> lefs, MaskT& mask) {
  static_assert(std::is_same_v<boost::dynamic_bitset<>, MaskT> ||
                    std::is_integral_v<std::decay_t<decltype(std::declval<MaskT&>().operator[](
                        std::declval<std::size_t>()))>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");
  assert(lefs.size() <= mask.size());
  for (auto i = 0UL; i < lefs.size(); ++i) {
    mask[i] = lefs[i].lifetime == 0;
  }
}

boost::asio::thread_pool Genome::instantiate_thread_pool() const {
  return boost::asio::thread_pool(this->_nthreads);
}

template <typename I>
boost::asio::thread_pool Genome::instantiate_thread_pool(I nthreads) {
  static_assert(std::is_integral_v<I>, "nthreads should have an integral type.");
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  return boost::asio::thread_pool(
      std::min(std::thread::hardware_concurrency(), static_cast<unsigned int>(nthreads)));
  DISABLE_WARNING_POP
}

}  // namespace modle
