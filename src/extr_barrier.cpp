#include "modle/extr_barrier.hpp"

namespace modle {
ExtrusionBarrier::ExtrusionBarrier(uint32_t abs_position, uint32_t bin_size, double prob_of_block,
                                   uint64_t seed)
    : _abs_position(abs_position),
      _rel_position(abs_to_rel_position(abs_position, bin_size)),
      _prob_of_block(prob_of_block),
      _rng(seed),
      _block(prob_of_block),
      _seed(seed) {}

double ExtrusionBarrier::get_prob_of_blocking() const { return this->_block.p(); }

uint32_t ExtrusionBarrier::abs_to_rel_position(uint32_t abs_pos, uint32_t bin_size) {
  const uint32_t bin_start = (abs_pos / bin_size) * bin_size;
  return abs_pos - bin_start;
}

void ExtrusionBarrier::activate() {
  this->_block = std::bernoulli_distribution(this->_prob_of_block);
}
void ExtrusionBarrier::inactivate() { this->_block = std::bernoulli_distribution(0); }

bool ExtrusionBarrier::is_active() const { return this->_block.p() > 0; }

bool ExtrusionBarrier::is_blocking() { return this->_block(this->_rng); }

void ExtrusionBarrier::relocate(uint32_t new_abs_pos, uint32_t bin_size) {
  this->_abs_position = new_abs_pos;
  this->_rel_position = abs_to_rel_position(new_abs_pos, bin_size);
}

void ExtrusionBarrier::relocate_and_activate(uint32_t new_abs_pos, uint32_t bin_size) {
  this->relocate(new_abs_pos, bin_size);
  this->activate();
}

uint32_t ExtrusionBarrier::get_abs_pos() const { return this->_abs_position; }
uint32_t ExtrusionBarrier::get_rel_pos() const { return this->_rel_position; }

}  // namespace modle
