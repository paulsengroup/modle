#pragma once

#include <cstdint>
#include <random>

namespace modle {

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
    std::uniform_int_distribution<uint64_t> uniform_rng(0, chr.length());
    const auto barrier_position = uniform_rng(rand_eng);
    const auto direction = strand_selector(rand_eng) ? DNA::Direction::rev : DNA::Direction::fwd;

    // Add the new extrusion barrier to the appropriate bin
    auto& bin = chr.dna.get_bin_from_pos(barrier_position);
    bin.add_extr_barrier(barrier_position, this->_probability_of_barrier_block, direction);
  }
  for (auto& chr : this->get_chromosomes()) {
    chr.sort_barriers_by_pos();
  }
}
}  // namespace modle
