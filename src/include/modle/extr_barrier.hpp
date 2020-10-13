#pragma once

#include <cstdint>
#include <random>

#include "modle/dna.hpp"

namespace modle {
/** \brief Simple class to model a loop extrusion barrier
 *
 * An instance of ExtrusionBarrier can have a direction (DNA::Direction) and a probability of
 * blocking. The probability of blocking is used to initialize a random number generator that
 * outputs numbers under the geometric distribution.
 * The numbers returned by this generator will be used as a base to determine for how long a given
 * ExtrusionUnit that reached this Extrusionbarrier instance should be stalled.
 */
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
