#include "modle/extr_barrier.hpp"

#include "absl/strings/str_format.h"

namespace modle {
ExtrusionBarrier::ExtrusionBarrier(uint32_t abs_position, uint32_t bin_size, double prob_of_block,
                                   uint64_t seed)
    : _abs_position(abs_position),
      _rel_position(abs_to_rel_position(abs_position, bin_size)),
      _rng(seed),
      _block(std::bernoulli_distribution(prob_of_block)),
      _seed(seed) {}

double ExtrusionBarrier::get_prob_of_blocking() const { return this->_block.p(); }

uint32_t ExtrusionBarrier::abs_to_rel_position(uint32_t abs_pos, uint32_t bin_size) {
  const uint32_t bin_start = (abs_pos / bin_size) * bin_size;
  return abs_pos - bin_start;
}

bool ExtrusionBarrier::is_blocking() { return this->_block(this->_rng); }

void ExtrusionBarrier::relocate(uint32_t new_abs_pos, uint32_t bin_size) {
  this->_abs_position = new_abs_pos;
  this->_rel_position = abs_to_rel_position(new_abs_pos, bin_size);
}

uint32_t ExtrusionBarrier::get_abs_pos() const { return this->_abs_position; }
uint32_t ExtrusionBarrier::get_rel_pos() const { return this->_rel_position; }

}  // namespace modle
