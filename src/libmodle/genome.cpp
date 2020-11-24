#include "modle/genome.hpp"

#include <boost/asio/thread_pool.hpp>
#include <cassert>
#include <filesystem>
#include <functional>

#include "fmt/printf.h"
#include "modle/bed.hpp"
#include "modle/chr_sizes.hpp"

namespace modle {

Genome::Genome(const config& c)
    : _seed(c.seed),
      _rand_eng(std::mt19937(_seed)),
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
      _randomize_contact_sampling(c.randomize_contact_sampling_interval),
      _sample_contacts(1.0 / _sampling_interval) {}

std::vector<uint32_t> Genome::get_chromosome_lengths() const {
  std::vector<uint32_t> lengths;
  lengths.reserve(this->get_n_chromosomes());
  for (const auto& chr : this->_chromosomes) {
    assert(chr.length() <= UINT32_MAX);  // NOLINT
    lengths.push_back(static_cast<uint32_t>(chr.length()));
  }
  return lengths;
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
  auto nthreads = std::thread::hardware_concurrency();
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

  for (auto i = 0UL; i < n_barriers; ++i) {
    // Randomly select a chromosome, barrier binding pos and direction
    auto& chr = this->_chromosomes[chr_idx(this->_rand_eng)];
    std::uniform_int_distribution<uint64_t> uniform_rng(0, chr.length());
    const auto barrier_position = uniform_rng(this->_rand_eng);
    const DNA::Direction direction =
        strand_selector(this->_rand_eng) ? DNA::Direction::rev : DNA::Direction::fwd;

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

void Genome::randomly_bind_lefs() {
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<std::size_t> chr_idx(weights.begin(), weights.end());

  for (auto& lef : this->_lefs) {
    // Randomly select a chromosome, barrier binding pos and direction
    Chromosome& chr = this->_chromosomes[chr_idx(this->_rand_eng)];
    lef.randomly_bind_to_chr(&chr, this->_rand_eng);  // And bind to it
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

uint32_t Genome::run_burnin(double prob_of_rebinding, uint32_t target_n_of_unload_events,
                            uint64_t min_extr_rounds) {
  double avg_num_of_extr_events_per_bind =
      (static_cast<double>(this->_avg_lef_processivity) / this->_bin_size) /
      2 /* N of active extr. unit */;
  const uint64_t n_of_lefs_to_bind_each_round = static_cast<uint64_t>(std::lround(
      std::max(static_cast<double>(this->get_n_lefs()) / avg_num_of_extr_events_per_bind, 1.0)));

  std::vector<uint16_t> unload_events(this->get_n_lefs(), 0);
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<std::size_t> chr_idx(weights.begin(), weights.end());

  uint64_t start_idx = 0;
  uint64_t end_idx = n_of_lefs_to_bind_each_round;
  for (auto rounds = 0U;; ++rounds) {
    for (; start_idx < end_idx && start_idx < this->get_n_lefs(); ++start_idx) {
      Chromosome& chr = this->_chromosomes[chr_idx(this->_rand_eng)];
      this->_lefs[start_idx].randomly_bind_to_chr(&chr, this->_rand_eng);
    }
    for (auto i = 0U; i < start_idx; ++i) {
      auto& lef = this->_lefs[i];

      if (lef.is_bound()) {
        lef.extrude(this->_rand_eng);
        unload_events[i] += !lef.is_bound();  // NOLINT(readability-implicit-bool-conversion)
      }
    }

    for (auto i = 0U; i < start_idx; ++i) {
      try {
        if (auto& lef = this->_lefs[i]; lef.is_bound()) {
          lef.check_constraints(this->_rand_eng);
        } else {
          auto& chr = this->_chromosomes[chr_idx(this->_rand_eng)];
          lef.try_rebind(chr, this->_rand_eng, prob_of_rebinding, false);
        }
      } catch (const std::runtime_error& err) {
        throw std::runtime_error(fmt::format("Error occurred at lef #%lu: %s!", i, err.what()));
      }
    }

    start_idx = std::min(end_idx, this->get_n_lefs());
    end_idx = std::min(end_idx + n_of_lefs_to_bind_each_round, this->get_n_lefs());

    if (rounds >= min_extr_rounds &&
        std::all_of(unload_events.begin(), unload_events.end(),
                    [&](uint16_t n) { return n >= target_n_of_unload_events; })) {
      return rounds;
    }
  }
}

void Genome::simulate_extrusion(uint32_t iterations) {
  const auto step = iterations / 500;
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<std::size_t> chr_idx(weights.begin(), weights.end());

  auto t0 = absl::Now();
  for (auto i = 1UL; i <= iterations; ++i) {
    bool register_contacts = this->_randomize_contact_sampling
                                 ? this->_sample_contacts(this->_rand_eng)
                                 : i % this->_sampling_interval == 0;

    for (auto& lef : this->_lefs) {
      if (lef.is_bound()) {  // Register contact and extrude if LEF is bound
        if (register_contacts) {
          lef.register_contact();
        }
        lef.extrude(this->_rand_eng);
      }
    }

    for (auto& lef : this->_lefs) {
      // Check whether the last round of extrusions caused a collision with another LEF or an
      // extr. boundary, apply the stall and/or increase LEF lifetime where appropriate
      if (lef.is_bound()) {
        lef.check_constraints(this->_rand_eng);
      } else {
        // Randomly select a chromosome and try to bind one of the free LEFs to it
        auto& chr = this->_chromosomes[chr_idx(this->_rand_eng)];
        lef.try_rebind(chr, this->_rand_eng, this->_probability_of_lef_rebind, register_contacts);
      }
    }

    if (i % step == 0) {
      fmt::fprintf(stderr, "Running iteration %lu/%lu (%.2f iterations/s)\n", i, iterations,
                   step / absl::ToDoubleSeconds(absl::Now() - t0));
      t0 = absl::Now();
    }
  }
}
}  // namespace modle