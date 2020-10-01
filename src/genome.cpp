#include "modle/genome.hpp"

#include <absl/strings/str_format.h>

#include <functional>

#include "modle/extr_barrier.hpp"
#include "modle/parsers.hpp"

namespace modle {

Genome::Genome(std::string_view path_to_bed, uint32_t bin_size, uint32_t n_lefs,
               uint32_t avg_lef_processivity, double probability_of_barrier_block, uint64_t seed)
    : _path_to_bed(path_to_bed),
      _bin_size(bin_size),
      _chromosomes(init_chromosomes_from_bed()),
      _rndev(std::default_random_engine(seed)),
      _lefs(generate_lefs(n_lefs, avg_lef_processivity)),
      _avg_lef_processivity(avg_lef_processivity),
      _probability_of_barrier_block(probability_of_barrier_block),
      _seed(seed) {}

std::vector<uint32_t> Genome::get_chromosome_lengths() const {
  std::vector<uint32_t> lengths;
  lengths.reserve(this->get_n_chromosomes());
  for (const auto& [chr, dna, _] : this->_chromosomes) lengths.push_back(dna.length());
  return lengths;
}

uint32_t Genome::get_n_of_free_lefs() const {
  return std::accumulate(
      this->_lefs.begin(), this->_lefs.end(), 0UL,
      [](uint32_t accumulator, const Lef& lef) { return accumulator + !lef.is_bound(); });
}

std::vector<Genome::Chromosome> Genome::init_chromosomes_from_bed() const {
  std::vector<Genome::Chromosome> chromosomes;
  SimpleBEDParser parser(this->_path_to_bed);
  SimpleBED record = parser.parse_next();
  do {
    chromosomes.emplace_back(record.chr, DNA{record.chr_end - record.chr_start, this->_bin_size});
    record = parser.parse_next();
  } while (!record.empty());

  return chromosomes;
}

std::vector<Lef> Genome::generate_lefs(uint32_t n, uint32_t avg_processivity) {
  return std::vector<Lef>{n, Lef{avg_processivity}};
}

void Genome::randomly_generate_barriers(uint32_t n_barriers, double prob_of_block) {
  std::hash<std::string_view> str_hasher;
  absl::btree_multimap<std::string_view, std::shared_ptr<ExtrusionBarrier>> barriers;
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());
  std::bernoulli_distribution strand_selector(0.5);
  for (uint32_t i = 0; i < n_barriers; ++i) {
    auto& [chr, dna, barriers] = this->_chromosomes[chr_idx(this->_rndev)];
    std::uniform_int_distribution<uint32_t> uniform_rng(0, dna.length());
    const auto barrier_position = uniform_rng(this->_rndev);
    // We are using the hash of the chr name + the position as seed for the rng that controls
    // whether a boundary will block or not
    auto barrier = barriers.emplace_back(barrier_position, this->_bin_size, prob_of_block,
                                         str_hasher(chr) + barrier_position);
    if (strand_selector(this->_rndev)) {
      dna.add_fwd_barrier(barrier, barrier_position);
    } else {
      dna.add_rev_barrier(barrier, barrier_position);
    }
  }
}

const std::vector<Genome::Chromosome>& Genome::get_chromosomes() const {
  return this->_chromosomes;
}

std::vector<Genome::Chromosome>& Genome::get_chromosomes() { return this->_chromosomes; }

std::vector<Lef>& Genome::get_lefs() { return this->_lefs; }
const std::vector<Lef>& Genome::get_lefs() const { return this->_lefs; }

uint32_t Genome::get_n_chromosomes() const { return this->_chromosomes.size(); }

uint64_t Genome::size() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0ULL,
                         [&](uint64_t accumulator, const Genome::Chromosome& chr) {
                           return accumulator + chr.length();
                         });
}

uint32_t Genome::n_bins() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
                         [](uint64_t accumulator, const Genome::Chromosome& chr) {
                           return accumulator + chr.n_bins();
                         });
}

uint32_t Genome::n_lefs() const { return this->_lefs.size(); }

uint32_t Genome::n_barriers() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
                         [](uint32_t accumulator, const Genome::Chromosome& chr) {
                           return accumulator + chr.n_barriers();
                         });
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
  for (auto lef_it = this->_lefs.begin(); lef_it != this->_lefs.end(); ++lef_it) {
    auto& [chr, dna, _] = this->_chromosomes[chr_idx(this->_rndev)];
    std::uniform_int_distribution<uint32_t> uniform_rng(0, dna.length());
    if (!lef_it->bind_at_pos(chr, dna, uniform_rng(this->_rndev)))
      --lef_it;  // Go to next lef only if binding was successful
  }
}

void Genome::simulate_extrusion(uint32_t iterations) {
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());
  auto t0 = absl::Now();
  for (uint32_t i = 0; i < iterations; ++i) {
    for (auto& lef : this->_lefs) {
      if (lef.is_bound()) lef.extrude();
    }
    for (auto& lef : this->_lefs) {
      if (lef.is_bound()) {
        lef.check_constrains();
        lef.try_unload(this->_rndev);
      }
      if (!lef.is_bound()) {
        auto& [chr, dna, _] = this->_chromosomes[chr_idx(this->_rndev)];
        lef.try_rebind(chr, dna, this->_rndev, 0.25);
      }
    }

    //  absl::FPrintF(stderr, "Solved %lu lef collisions.\n", collisions);
    if ((i + 1) % 1000 == 0) {
      absl::FPrintF(stderr, "1000 rounds of extrusion took %s\n",
                    absl::FormatDuration(absl::Now() - t0));
      absl::FPrintF(stderr, "# of free lefs: %lu/%lu\n", this->get_n_of_free_lefs(),
                    this->n_lefs());
      t0 = absl::Now();
    }
  }
}

Genome::Chromosome::Chromosome(std::string name, DNA dna)
    : name(std::move(name)), dna(std::move(dna)) {}

uint32_t Genome::Chromosome::length() const { return this->dna.length(); }
uint32_t Genome::Chromosome::n_bins() const { return this->dna.n_bins(); }
uint32_t Genome::Chromosome::n_barriers() const { return this->barriers.size(); }
}  // namespace modle