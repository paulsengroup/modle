#pragma once

#include <absl/container/flat_hash_map.h>
#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration, Time
#include <fmt/format.h>       // for FMT_STRING

#include <algorithm>  // for max, min, sort, transform, all_of, for_each, remove_if
#include <atomic>     // for atomic, memory_order_relaxed
#include <boost/asio/thread_pool.hpp>
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
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/view/chunk.hpp>  // for chunk_view, chunk
#include <stdexcept>                // for runtime_error
#include <thread>
#include <type_traits>  // for declval

#include "modle/bed.hpp"        // for BED, Parser, BED::BED6, BED::Standard
#include "modle/chr_sizes.hpp"  // for ChrSize, Parser
#include "modle/common.hpp"
#include "modle/config.hpp"    // for config
#include "modle/contacts.hpp"  // for ContactMatrix
#include "modle/cooler.hpp"
#include "modle/extr_barrier.hpp"  // IWYU pragma: keep
#include "modle/genome.hpp"
#include "modle/suppress_compiler_warnings.hpp"

namespace modle {

Genome::Genome(const config& c)
    : _seed(c.seed),
      _path_to_chrom_sizes_file(c.path_to_chr_sizes),
      _path_to_chr_subranges_file(c.path_to_chr_subranges),
      _bin_size(c.bin_size),
      _avg_lef_processivity(c.average_lef_processivity),
      _probability_of_barrier_block(c.probability_of_extrusion_barrier_block),
      _probability_of_lef_rebind(c.probability_of_lef_rebind),
      _probability_of_extr_unit_bypass(c.probability_of_extrusion_unit_bypass),
      _lef_unloader_strength_coeff(c.lef_unloader_strength),
      _lefs(generate_lefs(c.number_of_lefs)),
      _chromosomes(init_chromosomes_from_file(c.diagonal_width)),
      _sampling_interval(c.contact_sampling_interval),
      _randomize_contact_sampling(c.randomize_contact_sampling_interval),
      _nthreads(std::min(std::thread::hardware_concurrency(), c.nthreads)) {}

std::vector<uint64_t> Genome::get_chromosome_lengths() const {
  std::vector<uint64_t> lengths(this->get_n_chromosomes());
  std::transform(this->_chromosomes.begin(), this->_chromosomes.end(), lengths.begin(),
                 [](const auto& chr) { return chr.simulated_length(); });
  return lengths;
}

std::vector<double> Genome::get_chromosome_lef_affinities() const {
  std::vector<double> lef_affinities(this->get_n_chromosomes());
  std::transform(this->_chromosomes.begin(), this->_chromosomes.end(), lef_affinities.begin(),
                 [](const auto& chr) { return chr.get_total_lef_affinity(); });
  return lef_affinities;
}

std::vector<std::string_view> Genome::get_chromosome_names() const {
  std::vector<std::string_view> names(this->get_n_chromosomes());
  std::transform(this->_chromosomes.begin(), this->_chromosomes.end(), names.begin(),
                 [](const auto& chr) { return chr.name; });
  return names;
}

uint64_t Genome::get_n_lefs() const { return this->_lefs.size(); }

uint64_t Genome::get_n_of_free_lefs() const {
  return std::accumulate(
      this->_lefs.begin(), this->_lefs.end(), 0UL,
      // NOLINTNEXTLINE(readability-implicit-bool-conversion)
      [](uint64_t accumulator, const auto& lef) { return accumulator + !lef.is_bound(); });
}

uint64_t Genome::get_n_of_busy_lefs() const {
  assert(this->get_n_lefs() >= this->get_n_of_free_lefs());  // NOLINT
  return this->get_n_lefs() - this->get_n_of_free_lefs();
}

void Genome::write_contacts_to_file(std::string_view output_file) {
  std::vector<ContactMatrix<uint32_t>*> cmatrices(this->get_n_chromosomes());
  std::vector<std::string> chr_names(this->get_n_chromosomes());
  std::vector<uint64_t> chr_starts(this->get_n_chromosomes());
  std::vector<uint64_t> chr_ends(this->get_n_chromosomes());
  std::vector<uint64_t> chr_sizes(this->get_n_chromosomes());

  for (auto i = 0UL; i < this->get_n_chromosomes(); ++i) {
    auto& chr = this->_chromosomes[i];
    cmatrices[i] = &chr.contacts;
    chr_names[i] = chr.name;
    chr_starts[i] = chr.start;
    chr_ends[i] = chr.end;
    chr_sizes[i] = chr.total_length;
  }
  assert(!std::filesystem::exists(output_file));  // NOLINT
  std::filesystem::create_directories(std::filesystem::path(output_file).parent_path());
  cooler::Cooler c(output_file, cooler::Cooler::WRITE_ONLY, this->_bin_size);
  c.write_cmatrix_to_file(cmatrices, chr_names, chr_starts, chr_ends, chr_sizes);
}

void Genome::write_extrusion_barriers_to_file(std::string_view output_dir,
                                              bool force_overwrite) const {
  fmt::fprintf(stderr, "Writing extrusion barriers for %lu chromosomes in folder '%s'...",
               this->get_n_chromosomes(), output_dir);
  std::filesystem::create_directories(output_dir);
  auto t0 = absl::Now();
  std::for_each(this->_chromosomes.begin(), this->_chromosomes.end(), [&](const Chromosome& chr) {
    chr.write_barriers_to_tsv(output_dir, force_overwrite);
  });
  fmt::fprintf(stderr, "DONE! Written extrusion barrier coordinates for %lu chromosomes in %s\n",
               this->get_n_chromosomes(), absl::FormatDuration(absl::Now() - t0));
}

std::vector<Chromosome> Genome::init_chromosomes_from_file(uint32_t diagonal_width) const {
  // We are reading all the chrom_sizes at once to detect and deal with duplicate chrom_sizes
  auto chrom_sizes = modle::chr_sizes::Parser(this->_path_to_chrom_sizes_file).parse_all();
  absl::flat_hash_map<std::string, std::pair<uint64_t, uint64_t>> chrom_ranges;
  if (!this->_path_to_chr_subranges_file.empty()) {
    for (const auto& record : modle::bed::Parser(this->_path_to_chr_subranges_file).parse_all()) {
      chrom_ranges.emplace(record.chrom, std::make_pair(record.chrom_start, record.chrom_end));
    }
  }
  std::vector<Chromosome> chromosomes;
  chromosomes.reserve(chrom_sizes.size());
  std::transform(
      chrom_sizes.begin(), chrom_sizes.end(), std::back_inserter(chromosomes),
      [&](const auto& chr) {
        uint64_t start = chr.start;
        uint64_t end = chr.end;
        uint64_t length = end;

        if (chrom_ranges.contains(chr.name)) {
          const auto& match = chrom_ranges.at(chr.name);
          if (match.first < chr.start || match.second > chr.end) {
            throw std::runtime_error(fmt::format(
                FMT_STRING(
                    "According to the chrom.sizes file '{}', chromosome '{}' has a size of '{}', "
                    "but the "
                    "subrange specified through BED file '{}' extends past this region: range "
                    "{}:{}-{} does not fit in range {}:{}-{}"),
                this->_path_to_chrom_sizes_file, chr.name, chr.end,
                this->_path_to_chr_subranges_file, chr.name, match.first, match.second, chr.name,
                chr.start, chr.end));
          }
          start = match.first;
          end = match.second;
        }

        return Chromosome{chr.name, start, end, length, this->_bin_size, diagonal_width};
      });

  return chromosomes;
}

std::vector<Lef> Genome::generate_lefs(uint32_t n) {
  // TODO: Figure out how to clean up this function
  std::vector<Lef> v;
  v.reserve(n);
  for (auto i = 0U; i < n; ++i) {
    v.emplace_back(this->_bin_size, this->_avg_lef_processivity,
                   this->_probability_of_extr_unit_bypass, this->_lef_unloader_strength_coeff);
  }
  return v;
}

uint64_t Genome::n50() const {
  auto chr_lengths = this->get_chromosome_lengths();
  std::sort(chr_lengths.rbegin(), chr_lengths.rend(), std::greater<>());
  uint64_t n50_thresh = this->size() / 2;
  uint64_t tot = 0;
  return *ranges::find_if(chr_lengths,
                          [&](const auto len) { return ((tot += len) >= n50_thresh); });
}

std::pair<uint64_t, uint64_t> Genome::import_extrusion_barriers_from_bed(
    std::string_view path_to_bed, double probability_of_block) {
  auto p = modle::bed::Parser(path_to_bed, bed::BED::Standard::BED6);
  uint64_t nrecords = 0;
  uint64_t nrecords_ignored = 0;
  absl::flat_hash_map<std::string_view, Chromosome*> chromosomes;
  chromosomes.reserve(this->get_n_chromosomes());
  for (auto& chr : this->get_chromosomes()) {
    chromosomes.emplace(chr.name, &chr);
  }
  for (auto& record : p.parse_all()) {
    ++nrecords;
    if (!chromosomes.contains(record.chrom) || record.strand == '.' ||
        record.chrom_start < chromosomes[record.chrom]->get_start_pos() ||
        record.chrom_end > chromosomes[record.chrom]->get_end_pos()) {
      ++nrecords_ignored;
      continue;
    }
    if (probability_of_block != 0) {
      record.score = probability_of_block;
    }
    if (record.score < 0 || record.score > 1) {
      throw std::runtime_error(
          fmt::format("Invalid score field detected for record %s[%lu-%lu]: expected a score "
                      "between 0 and 1, got %.4g.",
                      record.name, record.chrom_start, record.chrom_end, record.score));
    }
    record.chrom_start -= chromosomes[record.chrom]->get_start_pos();
    record.chrom_end -= chromosomes[record.chrom]->get_start_pos();
    chromosomes[record.chrom]->dna.add_extr_barrier(record);
  }

  for (auto& chr : this->_chromosomes) {
    chr.sort_barriers_by_pos();
  }

  return {nrecords, nrecords_ignored};
}

const std::vector<Chromosome>& Genome::get_chromosomes() const { return this->_chromosomes; }
std::vector<Chromosome>& Genome::get_chromosomes() { return this->_chromosomes; }

std::vector<Lef>& Genome::get_lefs() { return this->_lefs; }
const std::vector<Lef>& Genome::get_lefs() const { return this->_lefs; }

uint32_t Genome::get_n_chromosomes() const {
  assert(this->_chromosomes.size() <= UINT32_MAX);  // NOLINT
  return static_cast<uint32_t>(this->_chromosomes.size());
}

uint64_t Genome::size() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0ULL,
      [&](uint64_t accumulator, const auto& chr) { return accumulator + chr.simulated_length(); });
}

uint64_t Genome::get_n_bins() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0ULL,
      [](uint64_t accumulator, const auto& chr) { return accumulator + chr.get_n_bins(); });
}

uint64_t Genome::get_n_barriers() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0ULL,
      [](uint64_t accumulator, const auto& chr) { return accumulator + chr.get_n_barriers(); });
}

void Genome::assign_lefs(bool bind_lefs_after_assignment) {
  absl::flat_hash_map<Chromosome*, double> chr_lef_affinities;
  double tot_affinity = 0;
  chr_lef_affinities.reserve(this->get_n_chromosomes());
  for (auto& chr : this->_chromosomes) {  // Build a map of the per-chromosome lef affinity
    const auto& [node, _] = chr_lef_affinities.emplace(&chr, chr.get_total_lef_affinity());
    tot_affinity += node->second;
  }
  std::vector<std::pair<Chromosome*, std::size_t>> chr_sorted_by_affinity(
      this->get_n_chromosomes());
  std::transform(this->_chromosomes.begin(), this->_chromosomes.end(),
                 chr_sorted_by_affinity.begin(),
                 [](auto& chr) { return std::make_pair(&chr, 0U); });
  std::sort(chr_sorted_by_affinity.begin(), chr_sorted_by_affinity.end(),
            [&](const auto& chr_pair1, const auto& chr_pair2) {
              return chr_lef_affinities[chr_pair1.first] > chr_lef_affinities[chr_pair2.first];
            });
  std::size_t lefs_assigned = 0;
  for (auto& [chr, nlefs] : chr_sorted_by_affinity) {
    nlefs = static_cast<std::size_t>(std::floor((chr_lef_affinities[chr] / tot_affinity) *
                                                static_cast<double>(this->get_n_lefs())));
    lefs_assigned += nlefs;
  }
  while (lefs_assigned < this->get_n_lefs()) {
    for (auto& [chr, nlefs] : chr_sorted_by_affinity) {
      if (lefs_assigned++ < this->get_n_lefs()) {
        ++nlefs;
      } else {
        break;
      }
    }
  }

  std::size_t i = 0;
  for (auto& [chr, nlefs] : chr_sorted_by_affinity) {
    chr->lefs.reserve(nlefs);
    for (auto j = 0U; j < nlefs; ++j) {
      auto* last_lef = chr->lefs.emplace_back(&(this->_lefs[i++]));
      last_lef->assign_to_chr(chr);
      if (bind_lefs_after_assignment) {
        last_lef->randomly_bind_to_chr(chr, chr->_rand_eng);
      }
    }
  }
}

uint64_t Genome::remove_chromosomes_wo_extr_barriers() {
  const auto n_chromosomes = this->get_n_chromosomes();
  this->_chromosomes.erase(
      std::remove_if(this->_chromosomes.begin(), this->_chromosomes.end(),
                     [](const auto& chr) { return chr.get_n_barriers() == 0; }),
      this->_chromosomes.end());
  return n_chromosomes - this->get_n_chromosomes();
}

std::pair<double, double> Genome::run_burnin(double prob_of_rebinding,
                                             uint32_t target_n_of_unload_events,
                                             uint64_t min_extr_rounds) {
  std::vector<double> burnin_rounds_performed(this->get_n_chromosomes());
  boost::asio::thread_pool tpool(std::min(this->get_n_chromosomes(), this->_nthreads));
  for (auto nchr = 0U; nchr < this->get_n_chromosomes(); ++nchr) {
    boost::asio::post(tpool, [&, nchr]() {
      const auto& chr = this->_chromosomes[nchr];
      std::seed_seq seed{this->_seed + std::hash<std::string>{}(chr.name) +
                         std::hash<uint64_t>{}(chr.simulated_length())};
      std::mt19937 rand_eng{seed};

      double avg_num_of_extr_events_per_bind =
          (static_cast<double>(this->_avg_lef_processivity) / chr.get_bin_size()) /
          2 /* N of active extr. unit */;
      const auto n_of_lefs_to_bind_each_round = std::lround(
          std::max(static_cast<double>(chr.get_n_lefs()) / avg_num_of_extr_events_per_bind, 1.0));

      std::vector<uint16_t> unload_events(chr.get_n_lefs(), 0);

      auto chunks = ranges::views::chunk(chr.lefs, n_of_lefs_to_bind_each_round);
      std::size_t lefs_loaded = 0;
      for (auto rounds = 0U;; ++rounds) {
        if (rounds < chunks.size()) {
          for (const auto& lef : chunks[rounds]) {
            lef->try_rebind(rand_eng);
            ++lefs_loaded;
          }
        }
        for (auto i = 0U; i < lefs_loaded; ++i) {
          auto& lef = *chr.lefs[i];
          if (lef.is_bound()) {
            lef.extrude(rand_eng);
            /* clang-format off
             * TODO: Figure out why GCC 8.3 is complaining about this line. The warning is:
             * warning: conversion from ‘int’ to ‘__gnu_cxx::__alloc_traits<std::allocator<short unsigned int>, short unsigned int>::value_type’ {aka ‘short unsigned int’} may change value [-Wconversion]
             *          unload_events[i] += static_cast<uint16_t>(!lef.is_bound());  // NOLINT(readability-implicit-bool-conversion)
             *                                                                   ^
             * clang-format on
             */
            unload_events[i] += static_cast<uint16_t>(!lef.is_bound());
          }
        }

        for (auto i = 0U; i < lefs_loaded; ++i) {
          if (auto& lef = *chr.lefs[i]; lef.is_bound()) {
            lef.check_constraints(rand_eng);
          } else {
            lef.try_rebind(rand_eng, prob_of_rebinding, false);
          }
        }
        if (rounds >= min_extr_rounds &&
            std::all_of(unload_events.begin(), unload_events.end(),
                        [&](uint16_t n) { return n >= target_n_of_unload_events; })) {
          burnin_rounds_performed[nchr] = static_cast<double>(rounds);
          return;
        }
      }
    });
  }
  tpool.join();

  if (burnin_rounds_performed.size() == 1) {
    return std::make_pair(burnin_rounds_performed[0], 0.0);
  }

  const auto tot_burnin_rounds =
      std::accumulate(burnin_rounds_performed.begin(), burnin_rounds_performed.end(), 0.0);
  const auto avg_burnin_rounds = tot_burnin_rounds / static_cast<double>(this->get_n_chromosomes());
  const auto burnin_rounds_stdev =
      std::sqrt(tot_burnin_rounds / static_cast<double>(burnin_rounds_performed.size() - 1));

  return std::make_pair(avg_burnin_rounds, burnin_rounds_stdev);
}

void Genome::simulate_extrusion(uint32_t iterations, double target_contact_density) {
  {
    std::vector<std::string_view> chr_names(this->get_n_chromosomes());
    std::transform(this->_chromosomes.begin(), this->_chromosomes.end(), chr_names.begin(),
                   [](const auto& c) {
                     return std::string_view{c.name.data(), c.name.size()};
                   });
    fmt::print(stderr, FMT_STRING("Simulating loop extrusion on the following chromosomes: '{}'\n"),
               absl::StrJoin(chr_names, "', '"));
  }

  // If the simulation is set to stop when a target contact density is reached, set the number of
  // iterations to a very large number (2^32)
  if (target_contact_density != 0.0) {
    iterations = UINT32_MAX;
  }

  // Instantiate the thread pool to be used for simulation and IO operations
  boost::asio::thread_pool tpool(std::min(this->_nthreads, this->get_n_chromosomes()));

  // Initialize variables for simulation progress tracking
  std::atomic<uint64_t> ticks_done{0};
  std::atomic<uint64_t> extrusion_events{0};
  std::atomic<uint64_t> chromosomes_completed{0};
  std::atomic<bool> simulation_completed{false};
  std::mutex m;  // This mutex is supposed to protect simulation_complete. As of C++17 we also need
  // to acquire a lock when modifying the condition variable simulation_completed_cv
  // (otherwise the change might not be communicated to waiting threads)
  std::condition_variable simulation_completed_cv;

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_CONVERSION
  DISABLE_WARNING_DOUBLE_PROMOTION
  std::thread progress_tracker([&]() {  // This thread is used to periodically print the simulation
    // progress to stderr
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
            ? std::accumulate(
                  this->_chromosomes.begin(), this->_chromosomes.end(), 0.0L,
                  [&](auto accumulator, const auto& chr) {
                    return accumulator +
                           (target_contact_density * chr.contacts.n_cols() * chr.contacts.n_rows());
                  })
            : iterations * this->get_n_chromosomes();

    const auto t0 = absl::Now();
    while (true) {
      {  // Wait for 5 seconds or until the simulation terminates
        std::unique_lock<std::mutex> lk(m);
        simulation_completed_cv.wait_for(lk, std::chrono::seconds(5));  // NOLINT
        if (simulation_completed) {  // Stop printing updates once there are no chromosomes left,
          // and return immediately, so that the main thread can join
          // this thread
          // Unlocking here is probably unnecessary, as there's nothing waiting on this mutex
          lk.unlock();
          return;
        }
      }
      if (extrusion_events > 0) {  // Print progress
        const auto progress = ticks_done / tot_ticks;
        const auto throughput = static_cast<double>(extrusion_events) / 5.0e6; /* 5s * 1M NOLINT */
        const auto delta_t = absl::ToDoubleSeconds(absl::Now() - t0);
        fmt::print(
            stderr,
            FMT_STRING("### ~{:.2f}% ###   {:.2f}M extr/sec - Simulation completed for "
                       "{}/{} chromosomes - ETA {}.\n"),
            100.0 * progress, throughput, chromosomes_completed, this->get_n_chromosomes(),
            absl::FormatDuration(absl::Seconds(
                (std::min(this->_nthreads, this->get_n_chromosomes()) /
                 std::min(static_cast<double>(this->_nthreads),
                          static_cast<double>(this->get_n_chromosomes() - chromosomes_completed))) *
                (delta_t / std::max(1.0e-06L, progress) - delta_t))));
        extrusion_events = 0;
      }
    }
  });
  DISABLE_WARNING_POP

  // Loop extrusion is simulated using boost::thread_pool, where each thread simulates loop
  // extrusion at the chromosome level. This constrains the level of parallelism to the number of
  // chromosomes that are being simulated. We can certainly do better.
  // This is just a way to get significant speedups with very little effort
  for (auto nchr = 0U; nchr < this->_chromosomes.size(); ++nchr) {
    boost::asio::post(tpool, [&, nchr]() {
      const auto t0 = absl::Now();
      auto& chr = this->_chromosomes[nchr];  // Alias for the chromosome that is being simulated
      // This random number generator is used to randomize sampling interval
      // (e.g. when --randomize-contact-sampling-interval is specified)
      std::bernoulli_distribution sample_contacts{1.0 / this->_sampling_interval};

      // Calculate the number of contacts after which the simulation for the current chromosome is
      // stopped. This is set to a very large number (2^64) when the simulation is set to run for a
      // fixed number of iterations
      const auto target_n_of_contacts =
          target_contact_density != 0.0
              ? target_contact_density *
                    static_cast<double>(chr.contacts.n_rows() * chr.contacts.n_cols())
              : std::numeric_limits<double>::max();

      // Variables to track the simulation progress for the current chromosome
      uint64_t local_extr_events_counter = 0;
      uint64_t ticks_local = 0;

      // This for loop is where the simulation actually takes place.
      // We basically keep iterating until one of the stopping condition is met
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
            lef->extrude(chr._rand_eng);
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

        // Add the number of extr. events to the global counter. This is used to calculate the n. of
        // extr. events per seconds, which is a good way to assess the simulation throughput
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
      if (++chromosomes_completed == this->get_n_chromosomes()) {
        {  // Notify the thread that is tracking the simulation progress that we are done simulating
          std::scoped_lock<std::mutex> lk(m);
          simulation_completed = true;
        }
        simulation_completed_cv.notify_all();
      }
      fmt::print(stderr, FMT_STRING("DONE simulating loop extrusion on '{}'! Simulation took {}\n"),
                 chr.name, absl::FormatDuration(absl::Now() - t0));
      /*
      if (write_contacts_to_file) {
        chr.write_contacts_to_tsv(output_dir, force_overwrite);
        chr.write_barriers_to_tsv(output_dir, force_overwrite);
        ++chromosomes_written_to_disk;
      }
       */
      // Free up resources
      // chr = Chromosome{"", 0, 0, 0};
    });
  }
  // This blocks until all the tasks posted to the thread_pool have been completed
  tpool.join();
  progress_tracker.join();
}

void Genome::simulate_extrusion() { this->simulate_extrusion(1, 0); }
void Genome::simulate_extrusion(uint32_t iterations) { this->simulate_extrusion(iterations, 0); }
void Genome::simulate_extrusion(double target_contact_density) {
  this->simulate_extrusion(0, target_contact_density);
}

template <typename I>
void Genome::randomly_generate_extrusion_barriers(I n_barriers, uint64_t seed) {
  static_assert(std::is_integral<I>::value, "I should be an integral numeric type.");

  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<std::size_t> chr_idx(weights.begin(), weights.end());
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  std::bernoulli_distribution strand_selector(0.5);
  std::seed_seq seeder{seed + n_barriers};
  std::mt19937 rand_eng{seeder};

  for (I i = 0UL; i < n_barriers; ++i) {
    // Randomly select a chromosome, barrier binding pos and direction
    auto& chr = this->_chromosomes[chr_idx(rand_eng)];
    std::uniform_int_distribution<uint64_t> uniform_rng(0, chr.simulated_length());
    const auto barrier_position = uniform_rng(rand_eng);
    const auto direction = strand_selector(rand_eng) ? dna::Direction::rev : dna::Direction::fwd;

    // Add the new extrusion barrier to the appropriate bin
    auto& bin = chr.dna.get_bin_from_pos(barrier_position);
    bin.add_extr_barrier(barrier_position, this->_probability_of_barrier_block, direction);
  }
  for (auto& chr : this->get_chromosomes()) {
    chr.sort_barriers_by_pos();
  }
}
}  // namespace modle
