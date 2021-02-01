#pragma once

#ifdef USE_XOSHIRO
#include <XoshiroCpp.hpp>
#endif
#include <cstdint>  // for uint64_t
#include <random>   // for geometric_distribution

#include "modle/common.hpp"

namespace modle {

#ifdef USE_XOSHIRO
using PRNG = XoshiroCpp::Xoshiro256PlusPlus;
using seeder = XoshiroCpp::SplitMix64;
#else
using PRNG = std::mt19937_64;
using seeder = std::seed_seq;
#endif

/** \brief Simple class to model a loop extrusion barrier
 *
 * An instance of ExtrusionBarrier can have a direction (dna::Direction) and a probability of
 * blocking. The probability of blocking is used to initialize a random number generator that
 * outputs numbers under the geometric distribution.
 * The numbers returned by this generator will be used as a base to determine for how long a given
 * ExtrusionUnit that reached this ExtrusionBarrier instance should be stalled.
 */
class ExtrusionBarrier {
 public:
  inline ExtrusionBarrier(std::size_t pos, double prob_of_block, dna::Direction direction);
  [[nodiscard]] inline std::size_t get_pos() const;
  [[nodiscard]] inline double get_prob_of_block() const;
  [[nodiscard]] inline dna::Direction get_direction_of_block() const;
  [[nodiscard]] inline dna::Direction get_motif_direction() const;
  [[nodiscard]] inline uint32_t generate_num_stalls(modle::PRNG& rand_eng);
  [[nodiscard]] inline bool operator<(const ExtrusionBarrier& other) const;

 private:
  std::size_t _pos;
  dna::Direction _direction{dna::Direction::none};
  std::geometric_distribution<uint32_t> _n_stall_generator;
};
}  // namespace modle

#include "../../extr_barrier_impl.hpp"
