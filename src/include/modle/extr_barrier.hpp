#pragma once

#include <cstdint>
#include <random>

#include "modle/dna.hpp"

namespace modle {
class ExtrusionBarrier {
 public:
  ExtrusionBarrier(double prob_of_block, DNA::Direction direction);
  [[nodiscard]] double get_prob_of_block() const;
  [[nodiscard]] DNA::Direction get_direction() const;
  [[nodiscard]] uint32_t generate_num_of_blocking_events(std::mt19937& rand_gen);

 private:
  DNA::Direction _direction{DNA::Direction::none};
  std::geometric_distribution<uint32_t> _block_for;
};
}  // namespace modle
