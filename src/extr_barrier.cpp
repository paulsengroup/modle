#include "modle/extr_barrier.hpp"

#include "absl/strings/str_format.h"

namespace modle {
ExtrusionBarrier::ExtrusionBarrier(double prob_of_block, DNA::Direction direction)
    : _direction(direction), _block_for(1 - prob_of_block) {}

double ExtrusionBarrier::get_prob_of_block() const { return this->_block_for.p(); }
DNA::Direction ExtrusionBarrier::get_direction() const {
  assert(this->_direction != DNA::Direction::both && this->_direction != DNA::Direction::none);
  return this->_direction;
}

uint32_t ExtrusionBarrier::generate_num_of_blocking_events(std::mt19937 &rand_dev) {
  return this->_block_for(rand_dev);
}

}  // namespace modle
