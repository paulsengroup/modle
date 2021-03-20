#pragma once

// IWYU pragma: private, include "modle/extrusion_factors.hpp"

#include <algorithm>  // for min, clamp
#include <cassert>    // for assert
#include <cstddef>    // IWYU pragma: keep for size_t
#include <limits>     // for numeric_limits

#include "modle/common.hpp"  // for Bp

namespace modle {

ExtrusionUnit::ExtrusionUnit(bp_t pos, bp_t nstalls_lef_lef, bp_t nstalls_lef_bar)
    : _pos(pos), _nstalls_lef_lef(nstalls_lef_lef), _nstalls_lef_bar(nstalls_lef_bar) {}

void ExtrusionUnit::bind_at_pos(bp_t pos) {
  this->_pos = pos;
  this->_nstalls_lef_lef = 0;
  this->_nstalls_lef_bar = 0;
}

bp_t ExtrusionUnit::pos() const { return this->_pos; }

bool ExtrusionUnit::stalled() const { return this->_nstalls_lef_lef + this->_nstalls_lef_bar > 0; }

void ExtrusionUnit::increment_stalls(bp_t nstalls) {
  assert(this->stalled());
  if (this->_nstalls_lef_lef > 0) {
    this->_nstalls_lef_lef +=
        std::min(nstalls, std::numeric_limits<bp_t>::max() - this->_nstalls_lef_lef);
  } else {
    this->_nstalls_lef_bar +=
        std::min(nstalls, std::numeric_limits<bp_t>::max() - this->_nstalls_lef_bar);
  }
}

void ExtrusionUnit::decrement_stalls(bp_t nstalls) {
  assert(this->stalled());
  if (this->_nstalls_lef_lef > 0) {
    this->_nstalls_lef_lef -= std::min(nstalls, this->_nstalls_lef_lef);
  } else {
    this->_nstalls_lef_bar -= std::min(nstalls, this->_nstalls_lef_bar);
  }
}

bp_t ExtrusionUnit::lef_lef_stalls() const { return this->_nstalls_lef_lef; }

bp_t ExtrusionUnit::lef_bar_stalls() const { return this->_nstalls_lef_bar; }

void ExtrusionUnit::operator-(bp_t n) { this->_pos = std::clamp(this->_pos - n, 0UL, this->_pos); }

void ExtrusionUnit::operator+(bp_t n) {
  this->_pos = std::clamp(this->_pos + n, this->_pos, std::numeric_limits<bp_t>::max());
}

bool ExtrusionUnit::operator<(const ExtrusionUnit& other) { return this->_pos < other._pos; }
bool ExtrusionUnit::operator<(std::size_t other_pos) { return this->_pos < other_pos; }

void ExtrusionUnit::release() {
  this->_pos = std::numeric_limits<bp_t>::max();
  this->_nstalls_lef_lef = 0;
  this->_nstalls_lef_bar = 0;
}

bool Lef::is_bound() const {
  assert((this->rev_unit._pos == std::numeric_limits<bp_t>::max()) ==
         (this->fwd_unit._pos == std::numeric_limits<bp_t>::max()));
  return this->lifetime != 0;
}

void Lef::bind_at_pos(bp_t pos, bp_t lifetime_) {
  assert(this->lifetime == 0);

  this->rev_unit.bind_at_pos(pos);
  this->fwd_unit.bind_at_pos(pos);

  this->lifetime = lifetime_;
}

void Lef::release() {
  this->rev_unit.release();
  this->fwd_unit.release();
  this->lifetime = 0;
}

void Lef::reset() { this->release(); }

}  // namespace modle
