#pragma once

// IWYU pragma: private, include "modle/extrusion_barriers.hpp"

#include <cassert>  // for assert

#include "modle/common.hpp"  // for Bp, Direction

namespace modle {
ExtrusionBarrier::ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                                   double transition_prob_non_blocking_to_non_blocking,
                                   dna::Direction motif_direction)
    : _pos(pos),
      _blocking_to_blocking_transition_prob(transition_prob_blocking_to_blocking),
      _non_blocking_to_non_blocking_transition_prob(transition_prob_non_blocking_to_non_blocking),
      _blocking_direction(motif_direction == dna::Direction::fwd ? dna::Direction::rev
                                                                 : dna::Direction::fwd) {
  // NOLINTNEXTLINE
  assert(motif_direction == dna::Direction::fwd || motif_direction == dna::Direction::rev);
}

ExtrusionBarrier::ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                                   double transition_prob_non_blocking_to_non_blocking,
                                   char motif_direction)
    : _pos(pos),
      _blocking_to_blocking_transition_prob(transition_prob_blocking_to_blocking),
      _non_blocking_to_non_blocking_transition_prob(transition_prob_non_blocking_to_non_blocking),
      _blocking_direction(motif_direction == '+' ? dna::Direction::rev : dna::Direction::fwd) {
  // NOLINTNEXTLINE
  assert(motif_direction == '+' || motif_direction == '-');
}

bp_t ExtrusionBarrier::pos() const { return this->_pos; }
double ExtrusionBarrier::prob_block_to_block() const {
  return this->_blocking_to_blocking_transition_prob;
}

double ExtrusionBarrier::prob_block_to_no_block() const { return 1.0 - prob_block_to_block(); }

double ExtrusionBarrier::prob_no_block_to_no_block() const {
  return this->_non_blocking_to_non_blocking_transition_prob;
}

double ExtrusionBarrier::prob_no_block_to_block() const {
  return 1.0 - prob_no_block_to_no_block();
}

dna::Direction ExtrusionBarrier::blocking_direction() const { return this->_blocking_direction; }

CTCF::State CTCF::next_state(CTCF::State current_state, double occupied_self_transition_prob,
                             double not_occupied_self_transition_prob, PRNG& rand_eng) {
  assert(occupied_self_transition_prob >= 0 && occupied_self_transition_prob <= 1);
  assert(not_occupied_self_transition_prob >= 0 && not_occupied_self_transition_prob <= 1);

  const auto p = CTCF::state_gen_t{0.0, 1.0}(rand_eng);
  if (current_state == NOT_OCCUPIED && p > not_occupied_self_transition_prob) {
    return OCCUPIED;
  } else if (p > occupied_self_transition_prob) {
    assert(current_state == OCCUPIED);
    return NOT_OCCUPIED;
  }

  return current_state;
}

}  // namespace modle
