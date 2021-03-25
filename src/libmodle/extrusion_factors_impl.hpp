#pragma once

// IWYU pragma: private, include "modle/extrusion_factors.hpp"

#include <algorithm>  // for min, clamp
#include <cassert>    // for assert
#include <cstddef>    // IWYU pragma: keep for size_t
#include <limits>     // for numeric_limits

#include "modle/common.hpp"  // for Bp

namespace modle {

ExtrusionUnit::ExtrusionUnit(bp_t pos) : _pos(pos) {}

void ExtrusionUnit::bind_at_pos(bp_t pos) { this->_pos = pos; }

bp_t ExtrusionUnit::pos() const { return this->_pos; }

bool ExtrusionUnit::operator<(const ExtrusionUnit& other) { return this->_pos < other._pos; }
bool ExtrusionUnit::operator<(std::size_t other_pos) { return this->_pos < other_pos; }

void ExtrusionUnit::release() { this->_pos = std::numeric_limits<bp_t>::max(); }

bool Lef::is_bound() const {
  assert((this->rev_unit._pos == std::numeric_limits<bp_t>::max()) ==
         (this->fwd_unit._pos == std::numeric_limits<bp_t>::max()));
  return this->active;
}

void Lef::bind_at_pos(bp_t pos) {
  this->rev_unit.bind_at_pos(pos);
  this->fwd_unit.bind_at_pos(pos);

  this->active = true;
}

void Lef::release() {
  this->rev_unit.release();
  this->fwd_unit.release();
  this->active = false;
}

void Lef::reset() { this->release(); }

}  // namespace modle
