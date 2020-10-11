#include "modle/genome.hpp"

#include <functional>

#include "absl/strings/str_format.h"
#include "modle/extr_barrier.hpp"
#include "modle/parsers.hpp"

namespace modle {

Genome::Genome(std::string_view path_to_bed, uint32_t bin_size, uint32_t n_lefs,
               uint32_t avg_lef_processivity, double probability_of_barrier_block,
               double probability_of_lef_rebind, double probability_of_extr_unit_bypass,
               uint64_t seed)
    : _path_to_bed(path_to_bed),
      _bin_size(bin_size),
      _avg_lef_processivity(avg_lef_processivity),
      _probability_of_barrier_block(probability_of_barrier_block),
      _probability_of_lef_rebind(probability_of_lef_rebind),
      _probability_of_extr_unit_bypass(probability_of_extr_unit_bypass),
      _lefs(generate_lefs(n_lefs)),
      _chromosomes(init_chromosomes_from_bed()),
      _seed(seed),
      _rand_gen(std::mt19937(_seed)) {}

std::vector<uint32_t> Genome::get_chromosome_lengths() const {
  std::vector<uint32_t> lengths;
  lengths.reserve(this->get_n_chromosomes());
  for (const auto& chr : this->_chromosomes) lengths.push_back(chr.length());
  return lengths;
}

uint32_t Genome::get_n_lefs() const { return this->_lefs.size(); }

uint32_t Genome::get_n_of_free_lefs() const {
  return std::accumulate(
      this->_lefs.begin(), this->_lefs.end(), 0UL,
      [](uint32_t accumulator, const Lef& lef) { return accumulator + !lef.is_bound(); });
}

uint32_t Genome::get_n_of_busy_lefs() const {
  return this->get_n_lefs() - this->get_n_of_free_lefs();
}

std::vector<Chromosome> Genome::init_chromosomes_from_bed() const {
  std::vector<Chromosome> chromosomes;
  SimpleBEDParser parser(this->_path_to_bed);
  SimpleBED record = parser.parse_next();

  do {
    chromosomes.emplace_back(record.chr, record.chr_end - record.chr_start, this->_bin_size,
                             this->_avg_lef_processivity);
    record = parser.parse_next();
  } while (!record.empty());

  return chromosomes;
}

std::vector<Lef> Genome::generate_lefs(uint32_t n) {
  std::vector<Lef> v;
  v.reserve(n);
  for (auto i = 0UL; i < n; ++i)
    v.emplace_back(this->_bin_size, this->_avg_lef_processivity,
                   this->_probability_of_extr_unit_bypass);
  return v;
}

void Genome::randomly_generate_barriers(uint32_t n_barriers) {
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());
  std::bernoulli_distribution strand_selector(0.5);

  for (auto i = 0UL; i < n_barriers; ++i) {
    // Randomly select a chromosome, barrier binding pos and direction
    auto& chr = this->_chromosomes[chr_idx(this->_rand_gen)];
    std::uniform_int_distribution<uint32_t> uniform_rng(0, chr.length());
    const auto barrier_position = uniform_rng(this->_rand_gen);
    const DNA::Direction direction =
        strand_selector(this->_rand_gen) ? DNA::Direction::rev : DNA::Direction::fwd;

    // Add the new extrusion barrier to the appropriate bin
    auto& bin = chr.dna.get_bin_from_pos(barrier_position);
    bin.add_extr_barrier(this->_probability_of_barrier_block, direction);
    chr.barriers.emplace_back(&bin.get_all_extr_barriers()->back());
  }
}

const std::vector<Chromosome>& Genome::get_chromosomes() const { return this->_chromosomes; }
std::vector<Chromosome>& Genome::get_chromosomes() { return this->_chromosomes; }

std::vector<Lef>& Genome::get_lefs() { return this->_lefs; }
const std::vector<Lef>& Genome::get_lefs() const { return this->_lefs; }

uint32_t Genome::get_n_chromosomes() const { return this->_chromosomes.size(); }

uint64_t Genome::size() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0ULL,
      [&](uint64_t accumulator, const Chromosome& chr) { return accumulator + chr.length(); });
}

uint32_t Genome::get_n_bins() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
      [](uint64_t accumulator, const Chromosome& chr) { return accumulator + chr.get_n_bins(); });
}

uint32_t Genome::get_n_barriers() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
                         [](uint32_t accumulator, const Chromosome& chr) {
                           return accumulator + chr.get_n_barriers();
                         });
}

void Genome::randomly_bind_lefs() {
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());

  for (auto& lef : this->_lefs) {
    auto& chr = this->_chromosomes[chr_idx(this->_rand_gen)];
    std::uniform_int_distribution<uint32_t> uniform_rng(0, chr.length());
    auto pos = uniform_rng(this->_rand_gen);

    lef.bind_at_pos(chr, pos, this->_rand_gen);
  }
}

void Genome::simulate_extrusion(uint32_t iterations) {
  const auto step = iterations / 500;
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());

  auto t0 = absl::Now();
  for (auto i = 1UL; i <= iterations; ++i) {
    for (auto& lef : this->_lefs) {
      if (lef.is_bound()) {  // Register contact and extrude if LEF is bound
        lef.register_contact();
        lef.extrude();
      }
    }

    for (auto& lef : this->_lefs) {
      // Check whether the last round of extrusions caused a collision with another LEF or an extr.
      // boundary, apply the appropriate stall and increase LEF lifetime when appropriate
      if (lef.is_bound()) {
        lef.check_constraints(this->_rand_gen);
      } else {
        // Randomly select a chromosome and try to bind one of the free LEFs to it
        auto& chr = this->_chromosomes[chr_idx(this->_rand_gen)];
        lef.try_rebind(chr, this->_rand_gen, this->_probability_of_lef_rebind);
      }
    }

    if (i % step == 0) {
      absl::FPrintF(stderr, "Running iteration %lu/%lu (%.2f iterations/s)\n", i, iterations,
                    step / absl::ToDoubleSeconds(absl::Now() - t0));
      t0 = absl::Now();
    }
  }
}
}  // namespace modle