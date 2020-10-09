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
      _rand_gen(std::mt19937(seed)) {
}

std::vector<uint32_t> Genome::get_chromosome_lengths() const {
  std::vector<uint32_t> lengths;
  lengths.reserve(this->get_n_chromosomes());
  for (const auto& [chr, dna, x, y] : this->_chromosomes) lengths.push_back(dna.length());
  return lengths;
}

uint32_t Genome::get_n_of_free_lefs() const {
  return std::accumulate(
      this->_lefs.begin(), this->_lefs.end(), 0UL,
      [](uint32_t accumulator, const Lef& lef) { return accumulator + !lef.is_bound(); });
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
    auto& chr = this->_chromosomes[chr_idx(this->_rand_gen)];
    std::uniform_int_distribution<uint32_t> uniform_rng(0, chr.length());
    const auto barrier_position = uniform_rng(this->_rand_gen);
    auto& bin = chr.dna.get_bin_from_pos(barrier_position);
    const DNA::Direction direction =
        strand_selector(this->_rand_gen) ? DNA::Direction::rev : DNA::Direction::fwd;
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

uint32_t Genome::n_bins() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
      [](uint64_t accumulator, const Chromosome& chr) { return accumulator + chr.n_bins(); });
}

uint32_t Genome::n_lefs() const { return this->_lefs.size(); }

uint32_t Genome::n_barriers() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
      [](uint32_t accumulator, const Chromosome& chr) { return accumulator + chr.n_barriers(); });
}

/*
uint32_t Genome::bind_free_lefs() {
  std::vector<Lef*> free_lefs;
  for (auto& lef : this->_lefs) {
    if (!lef.is_bound()) {
      free_lefs.push_back(&lef);
    }
  }
  auto lef = free_lefs.begin();
  for (auto& [chr, dna] : this->_chromosomes) {
  }
}
*/

void Genome::randomly_bind_lefs() {
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());
  for (auto& lef : this->_lefs) {
    auto& chr = this->_chromosomes[chr_idx(this->_rand_gen)];
    std::uniform_int_distribution<uint32_t> uniform_rng(0, chr.length());
    const auto p = uniform_rng(this->_rand_gen);
    absl::FPrintF(stderr, "Binding LEF at pos %lu.\n", p);
    lef.bind_at_pos(chr, p, this->_rand_gen);
  }
}

void Genome::simulate_extrusion(uint32_t iterations) {
  const auto step = iterations / 1000;
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());
  auto t0 = absl::Now();
  for (uint32_t i = 0; i < iterations; ++i) {
    for (auto& lef : this->_lefs) {
      if (lef.is_bound()) {
        lef.register_contact();
        lef.extrude();
        /*
        absl::FPrintF(stderr, "lpos=%lu; rpos=%lu; l_stall=%s; r_stall=%s; loop_size=%lu\n",
                      lef.get_left_pos(), lef.get_right_pos(),
                      lef.left_is_stalled() ? "True" : "False",
                      lef.right_is_stalled() ? "True" : "False", lef.get_loop_size());
        usleep(100'000);
        */
      }
    }
    for (auto& lef : this->_lefs) {
      if (lef.is_bound()) {
        //        lef.register_contact();
        lef.check_constraints(this->_rand_gen);
      } else {
        auto& chr = this->_chromosomes[chr_idx(this->_rand_gen)];
        lef.try_rebind(chr, this->_rand_gen, this->_probability_of_lef_rebind);
      }
    }

    //  absl::FPrintF(stderr, "Solved %lu lef collisions.\n", collisions);
    if ((i + 1) % step == 0) {
      absl::FPrintF(stderr, "Running iteration %lu (%.2f iterations/s)\n", i + 1,
                    step / absl::ToDoubleSeconds(absl::Now() - t0));
      t0 = absl::Now();
    }
  }
}
}  // namespace modle