#pragma once

#include <random>

#include "modle/common.hpp"

namespace modle {

class ExtrusionBarrier {
 public:
  inline ExtrusionBarrier() = default;
  inline ExtrusionBarrier(bp_t pos, double prob_of_block, dna::Direction motif_direction);
  inline ExtrusionBarrier(bp_t pos, double prob_of_block, char motif_direction);

  [[nodiscard]] inline bp_t pos() const;
  [[nodiscard]] inline double pblock() const;
  [[nodiscard]] inline dna::Direction blocking_direction() const;

 protected:
  bp_t _pos;
  double _prob_of_block;
  dna::Direction _blocking_direction;
};

namespace CTCF {
enum State : uint8_t { NOT_OCCUPIED = 0, OCCUPIED = 1 };
[[nodiscard]] inline State next_state(State current_state, double occupied_self_transition_prob,
                                      double not_occupied_self_transition_prob, PRNG& rand_eng);
using state_gen_t = std::uniform_real_distribution<double>;
}  // namespace CTCF
}  // namespace modle

#include "../../extrusion_barriers_impl.hpp"  // IWYU pragma: keep
