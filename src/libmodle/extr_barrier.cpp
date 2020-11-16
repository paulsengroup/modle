#include "modle/extr_barrier.hpp"

#include <cassert>
#include <cstdint>
#include <random>

namespace modle {
ExtrusionBarrier::ExtrusionBarrier(uint64_t pos, double prob_of_block, DNA::Direction direction)
    : _pos(pos), _direction(direction), _n_stall_generator(1 - prob_of_block) {
  assert(_direction == DNA::Direction::fwd || _direction == DNA::Direction::rev);  // NOLINT
}

uint64_t ExtrusionBarrier::get_pos() const { return this->_pos; }

double ExtrusionBarrier::get_prob_of_block() const { return this->_n_stall_generator.p(); }

DNA::Direction ExtrusionBarrier::get_motif_direction() const {
  assert(this->_direction == DNA::Direction::fwd ||  // NOLINT
         this->_direction == DNA::Direction::rev);   // NOLINT
  return this->_direction;
}

DNA::Direction ExtrusionBarrier::get_direction_of_block() const {
  assert(this->_direction == DNA::Direction::fwd ||  // NOLINT
         this->_direction == DNA::Direction::rev);   // NOLINT
  return this->_direction == DNA::Direction::fwd ? DNA::Direction::rev : DNA::Direction::fwd;
}

uint32_t ExtrusionBarrier::generate_num_stalls(std::mt19937 &rand_eng) {
  return this->_n_stall_generator(rand_eng);
}

bool ExtrusionBarrier::operator<(const ExtrusionBarrier &other) const {
  return this->_pos < other._pos;
}

}  // namespace modle
