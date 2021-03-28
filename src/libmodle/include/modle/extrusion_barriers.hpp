#pragma once

#include <random>

#include "modle/common.hpp"

namespace modle {

class ExtrusionUnit;

class ExtrusionBarrier {
 public:
  inline ExtrusionBarrier() = default;
  inline ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                          double transition_prob_non_blocking_to_non_blocking,
                          dna::Direction motif_direction);
  inline ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                          double transition_prob_non_blocking_to_non_blocking,
                          char motif_direction);

  [[nodiscard]] inline bp_t pos() const;
  [[nodiscard]] inline double prob_block_to_block() const;
  [[nodiscard]] inline double prob_block_to_no_block() const;
  [[nodiscard]] inline double prob_no_block_to_no_block() const;
  [[nodiscard]] inline double prob_no_block_to_block() const;
  [[nodiscard]] inline dna::Direction blocking_direction() const;
  [[nodiscard]] inline bool operator<(const ExtrusionBarrier& other) const;

 protected:
  bp_t _pos;
  double _blocking_to_blocking_transition_prob;
  double _non_blocking_to_non_blocking_transition_prob;
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
