#pragma once

#include <cassert>

#include "modle/common.hpp"

namespace modle {
ExtrusionBarrier::ExtrusionBarrier(Bp pos, double prob_of_block, dna::Direction motif_direction)
    : _pos(pos),
      _prob_of_block(prob_of_block),
      _blocking_direction(motif_direction == dna::Direction::fwd ? dna::Direction::rev
                                                                 : dna::Direction::fwd) {
  // NOLINTNEXTLINE
  assert(motif_direction == dna::Direction::fwd || motif_direction == dna::Direction::rev);
}

ExtrusionBarrier::ExtrusionBarrier(Bp pos, double prob_of_block, char motif_direction)
    : _pos(pos),
      _prob_of_block(prob_of_block),
      _blocking_direction(motif_direction == '+' ? dna::Direction::rev : dna::Direction::fwd) {
  // NOLINTNEXTLINE
  assert(motif_direction == '+' || motif_direction == '-');
}

Bp ExtrusionBarrier::pos() const { return this->_pos; }
double ExtrusionBarrier::pblock() const { return this->_prob_of_block; }
dna::Direction ExtrusionBarrier::blocking_direction() const { return this->_blocking_direction; }

}  // namespace modle
