#include "modle/genome.hpp"

#include <boost/asio/thread_pool.hpp>
#include <cassert>
#include <condition_variable>
#include <filesystem>
#include <functional>

#include "absl/container/flat_hash_map.h"
#include "fmt/printf.h"
#include "modle/bed.hpp"
#include "modle/chr_sizes.hpp"
#include "range/v3/view/chunk.hpp"

namespace modle {

Genome::Genome(const config& c)
    : _seed(c.seed),
      _path_to_chr_size_file(c.path_to_chr_sizes),
      _bin_size(c.bin_size),
      _avg_lef_processivity(c.average_lef_processivity),
      _probability_of_barrier_block(c.probability_of_extrusion_barrier_block),
      _probability_of_lef_rebind(c.probability_of_lef_rebind),
      _probability_of_extr_unit_bypass(c.probability_of_extrusion_unit_bypass),
      _lef_unloader_strength_coeff(c.lef_unloader_strength),
      _lefs(generate_lefs(c.number_of_lefs)),
      _chromosomes(init_chromosomes_from_file(c.diagonal_width)),
      _sampling_interval(c.contact_sampling_interval),
      _randomize_contact_sampling(c.randomize_contact_sampling_interval) {}

std::vector<uint32_t> Genome::get_chromosome_lengths() const {
  std::vector<uint32_t> lengths;
  lengths.reserve(this->get_n_chromosomes());
  for (const auto& chr : this->_chromosomes) {
    assert(chr.length() <= UINT32_MAX);  // NOLINT
    lengths.push_back(static_cast<uint32_t>(chr.length()));
  }
  return lengths;
}

std::vector<double> Genome::get_chromosome_lef_affinities() const {
  std::vector<double> lef_affinities(this->get_n_chromosomes());
  std::transform(this->_chromosomes.begin(), this->_chromosomes.end(), lef_affinities.begin(),
                 [](const auto& chr) { return chr.get_total_lef_affinity(); });
  return lef_affinities;
}

std::vector<std::string_view> Genome::get_chromosome_names() const {
  std::vector<std::string_view> names;
  names.reserve(this->get_n_chromosomes());
  for (const auto& chr : this->_chromosomes) {
    names.push_back(chr.name);
  }
  return names;
}

uint64_t Genome::get_n_lefs() const { return this->_lefs.size(); }

uint64_t Genome::get_n_of_free_lefs() const {
  return std::accumulate(
      this->_lefs.begin(), this->_lefs.end(), 0UL,
      // NOLINTNEXTLINE(readability-implicit-bool-conversion)
      [](uint32_t accumulator, const Lef& lef) { return accumulator + !lef.is_bound(); });
}

uint64_t Genome::get_n_of_busy_lefs() const {
  assert(this->get_n_lefs() >= this->get_n_of_free_lefs());  // NOLINT
  return this->get_n_lefs() - this->get_n_of_free_lefs();
}

void Genome::write_contacts_to_file(std::string_view output_dir, bool force_overwrite) const {
  auto nthreads = std::min(std::thread::hardware_concurrency(), this->get_n_chromosomes());
  boost::asio::thread_pool tpool(nthreads);
  fmt::print(stderr,
             "Writing contact matrices for {} chromosome(s) in folder '{} using {} threads'...\n",
             this->get_n_chromosomes(), output_dir, nthreads);
  auto t0 = absl::Now();
  std::filesystem::create_directories(output_dir);
  for (const auto& chr : this->_chromosomes) {
    boost::asio::post(tpool, [&]() { chr.write_contacts_to_tsv(output_dir, force_overwrite); });
  }
  tpool.join();
  fmt::fprintf(stderr, "DONE! Saved %lu contact matrices in %s\n", this->get_n_chromosomes(),
               absl::FormatDuration(absl::Now() - t0));
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
  std::vector<Chromosome> chromosomes;
  auto parser = modle::chr_sizes::Parser(this->_path_to_chr_size_file);
  for (const auto& chr : parser.parse()) {
    chromosomes.emplace_back(chr.name, chr.start, chr.end, this->_bin_size, diagonal_width);
  }
  return chromosomes;
}

std::vector<Lef> Genome::generate_lefs(uint32_t n) {
  std::vector<Lef> v;
  v.reserve(n);
  for (auto i = 0UL; i < n; ++i) {
    v.emplace_back(this->_bin_size, this->_avg_lef_processivity,
                   this->_probability_of_extr_unit_bypass, this->_lef_unloader_strength_coeff);
  }
  return v;
}

uint64_t Genome::n50() const {
  auto chr_lengths = this->get_chromosome_lengths();
  uint64_t n50_thresh = this->size() / 2;
  std::sort(chr_lengths.rbegin(), chr_lengths.rend(), std::greater<>());
  uint64_t tot = 0;
  for (const auto& len : chr_lengths) {
    if ((tot += len) >= n50_thresh) {
      return len;
    }
  }
  return chr_lengths.back();
}

void Genome::randomly_generate_extrusion_barriers(uint32_t n_barriers) {
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<std::size_t> chr_idx(weights.begin(), weights.end());
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  std::bernoulli_distribution strand_selector(0.5);
  std::seed_seq seeder{n_barriers};
  std::mt19937 rand_eng{seeder};

  for (auto i = 0UL; i < n_barriers; ++i) {
    // Randomly select a chromosome, barrier binding pos and direction
    auto& chr = this->_chromosomes[chr_idx(rand_eng)];
    std::uniform_int_distribution<uint64_t> uniform_rng(0, chr.length());
    const auto barrier_position = uniform_rng(rand_eng);
    const DNA::Direction direction =
        strand_selector(rand_eng) ? DNA::Direction::rev : DNA::Direction::fwd;

    // Add the new extrusion barrier to the appropriate bin
    auto& bin = chr.dna.get_bin_from_pos(barrier_position);
    bin.add_extr_barrier(barrier_position, this->_probability_of_barrier_block, direction);
  }
  for (auto& chr : this->get_chromosomes()) {
    chr.sort_barriers_by_pos();
  }
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
    for (const auto& bin : chr.dna) {
      if (bin.has_extr_barrier()) {
        for (auto& b : bin.get_all_extr_barriers()) {
          chr.barriers.push_back(&b);
        }
      }
    }
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
      [&](uint64_t accumulator, const Chromosome& chr) { return accumulator + chr.length(); });
}

uint64_t Genome::get_n_bins() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0ULL,
      [](uint64_t accumulator, const Chromosome& chr) { return accumulator + chr.get_n_bins(); });
}

uint64_t Genome::get_n_barriers() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0ULL,
                         [](uint64_t accumulator, const Chromosome& chr) {
                           return accumulator + chr.get_n_barriers();
                         });
}

// Make sure we are not taking a ptr from a reference
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
                     [](const Chromosome& chr) { return chr.get_n_barriers() == 0; }),
      this->_chromosomes.end());
  return n_chromosomes - this->get_n_chromosomes();
}

std::pair<double, double> Genome::run_burnin(double prob_of_rebinding,
                                             uint32_t target_n_of_unload_events,
                                             uint64_t min_extr_rounds) {
  auto nthreads = std::min(std::thread::hardware_concurrency(), this->get_n_chromosomes());
  boost::asio::thread_pool tpool(nthreads);
  std::vector<double> burnin_rounds(this->get_n_chromosomes());

  for (auto nchr = 0U; nchr < this->get_n_chromosomes(); ++nchr) {
    boost::asio::post(tpool, [&, nchr]() {
      const auto& chr = this->_chromosomes[nchr];
      std::seed_seq seed{this->_seed + std::hash<std::string>{}(chr.name) +
                         std::hash<uint64_t>{}(chr.length())};
      std::mt19937 rand_eng{seed};

      double avg_num_of_extr_events_per_bind =
          (static_cast<double>(this->_avg_lef_processivity) / chr.get_bin_size()) /
          2 /* N of active extr. unit */;
      const uint64_t n_of_lefs_to_bind_each_round = static_cast<uint64_t>(std::lround(
          std::max(static_cast<double>(chr.get_n_lefs()) / avg_num_of_extr_events_per_bind, 1.0)));

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
            unload_events[i] += !lef.is_bound();  // NOLINT(readability-implicit-bool-conversion)
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
          burnin_rounds[nchr] = static_cast<double>(rounds);
          return;
        }
      }
    });
  }
  tpool.join();

  if (burnin_rounds.size() == 1) {
    return std::make_pair(burnin_rounds[0], 0.0);
  }

  const auto avg_burnin_rounds = std::accumulate(burnin_rounds.begin(), burnin_rounds.end(), 0.0) /
                                 static_cast<double>(this->get_n_chromosomes());

  const auto burnin_rounds_stdev = std::sqrt(
      std::accumulate(burnin_rounds.begin(), burnin_rounds.end(), 0.0,
                      [&](auto accumulator, auto rounds) { return accumulator + rounds; }) /
      static_cast<double>(burnin_rounds.size() - 1));

  return std::make_pair(avg_burnin_rounds, burnin_rounds_stdev);
}

void Genome::simulate_extrusion(uint32_t iterations, double target_contact_density) {
  if (target_contact_density != 0.0) {
    iterations = UINT32_MAX;
  }
  auto nthreads = std::min(std::thread::hardware_concurrency(), this->get_n_chromosomes());
  std::atomic<uint64_t> ticks_done{0};
  std::atomic<uint64_t> extrusion_events{0};
  bool simulation_completed{false};
  std::condition_variable cv;
  std::mutex m;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-int-float-conversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
  std::thread progress_tracker([&]() {
    const uint64_t tot_ticks =
        target_contact_density != 0.0
            ? std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0U,
                              [&](auto accumulator, const auto& chr) {
                                return accumulator +
                                       std::llround(target_contact_density * chr.contacts.n_cols() *
                                                    chr.contacts.n_rows());
                              })
            : iterations * this->get_n_chromosomes();

    while (true) {
      {
        std::unique_lock<std::mutex> lk(m);
        cv.wait_for(lk, std::chrono::seconds(5));  // NOLINT
        if (simulation_completed) {
          lk.unlock();  // Probably unnecessary, as there's nothing waiting on this mutex
          return;
        }
      }
      fmt::print(stderr, FMT_STRING("Approx. {:.2f}% done ({:.2f}M extr. events/s)\n"),
                 100.0 * ticks_done / tot_ticks, extrusion_events / 5.0e6 /* 5s * 1M */);  // NOLINT
      extrusion_events = 0;
    }
  });
#pragma GCC diagnostic pop

  boost::asio::thread_pool tpool(nthreads);
  for (auto nchr = 0U; nchr < this->_chromosomes.size(); ++nchr) {
    boost::asio::post(tpool, [&, nchr]() {
      const auto t0 = absl::Now();
      auto& chr = this->_chromosomes[nchr];
      std::bernoulli_distribution sample_contacts{1.0 / this->_sampling_interval};

      const uint64_t target_n_of_contacts =
          target_contact_density != 0.0
              ? static_cast<uint64_t>(std::llround(
                    target_contact_density *
                    static_cast<double>(chr.contacts.n_rows() * chr.contacts.n_cols())))
              : UINT64_MAX;

      uint64_t extr_events_local = 0;
      uint64_t ticks_local = 0;
      for (auto i = 1UL; i <= iterations; ++i) {
        bool register_contacts = this->_randomize_contact_sampling
                                     ? sample_contacts(chr._rand_eng)
                                     : i % this->_sampling_interval == 0;
        for (auto& lef : chr.lefs) {
          if (lef->is_bound()) {  // Register contact and extrude if LEF is bound
            if (register_contacts) {
              lef->register_contact();
            }
            lef->extrude(chr._rand_eng);
            ++extr_events_local;
          }
        }
        for (auto& lef : chr.lefs) {
          // Check whether the last round of extrusions caused a collision with another LEF or an
          // extr. boundary, apply the stall and/or increase LEF lifetime where appropriate
          if (lef->is_bound()) {
            lef->check_constraints(chr._rand_eng);
          } else {
            // Randomly select a chromosome and try to bind one of the free LEFs to it
            lef->try_rebind(chr._rand_eng, this->_probability_of_lef_rebind, register_contacts);
          }
        }

        extrusion_events.fetch_add(extr_events_local, std::memory_order_relaxed);
        extr_events_local = 0;

        if (register_contacts) {
          if (target_contact_density != 0.0) {
            assert(chr.contacts.get_tot_contacts() >= ticks_local);
            ticks_done.fetch_add(chr.contacts.get_tot_contacts() - ticks_local,
                                 std::memory_order_relaxed);
            ticks_local = chr.contacts.get_tot_contacts();
          } else {
            assert(i >= ticks_local);
            ticks_done.fetch_add(i - ticks_local, std::memory_order_relaxed);
            ticks_local = i;
          }
          if (chr.contacts.get_tot_contacts() >= target_n_of_contacts) {
            break;
          }
        }
      }
      fmt::print(stderr, FMT_STRING("DONE simulating loop extrusion on '{}'! Simulation took {}\n"),
                 chr.name, absl::FormatDuration(absl::Now() - t0));
    });
  }
  tpool.join();

  {
    std::scoped_lock<std::mutex> lk(m);
    simulation_completed = true;
  }
  cv.notify_all();
  progress_tracker.join();
}
}  // namespace modle