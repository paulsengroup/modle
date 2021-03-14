#pragma once

#include <cstdint>

#include "modle/common.hpp"

namespace modle {

ExtrusionUnit::ExtrusionUnit(Bp pos, Bp nstalls_lef_lef, Bp nstalls_lef_bar)
    : _pos(pos), _nstalls_lef_lef(nstalls_lef_lef), _nstalls_lef_bar(nstalls_lef_bar) {}

void ExtrusionUnit::bind_at_pos(Bp pos) {
  this->_pos = pos;
  this->_nstalls_lef_lef = 0;
  this->_nstalls_lef_bar = 0;
}

Bp ExtrusionUnit::pos() const { return this->_pos; }

bool ExtrusionUnit::stalled() const { return this->_nstalls_lef_lef + this->_nstalls_lef_bar > 0; }

void ExtrusionUnit::increment_stalls(Bp nstalls) {
  assert(this->stalled());
  if (this->_nstalls_lef_lef > 0) {
    this->_nstalls_lef_lef +=
        std::min(nstalls, std::numeric_limits<Bp>::max() - this->_nstalls_lef_lef);
  } else {
    this->_nstalls_lef_bar +=
        std::min(nstalls, std::numeric_limits<Bp>::max() - this->_nstalls_lef_bar);
  }
}

void ExtrusionUnit::decrement_stalls(Bp nstalls) {
  assert(this->stalled());
  if (this->_nstalls_lef_lef > 0) {
    this->_nstalls_lef_lef -= std::min(nstalls, this->_nstalls_lef_lef);
  } else {
    this->_nstalls_lef_bar -= std::min(nstalls, this->_nstalls_lef_bar);
  }
}

Bp ExtrusionUnit::lef_lef_stalls() const { return this->_nstalls_lef_lef; }

Bp ExtrusionUnit::lef_bar_stalls() const { return this->_nstalls_lef_bar; }

void ExtrusionUnit::operator-(Bp n) { this->_pos = std::clamp(this->_pos - n, 0UL, this->_pos); }

void ExtrusionUnit::operator+(Bp n) {
  this->_pos = std::clamp(this->_pos + n, this->_pos, std::numeric_limits<Bp>::max());
}

bool ExtrusionUnit::operator<(const ExtrusionUnit& other) { return this->_pos < other._pos; }
bool ExtrusionUnit::operator<(std::size_t other_pos) { return this->_pos < other_pos; }

void ExtrusionUnit::release() {
  this->_pos = std::numeric_limits<Bp>::max();
  this->_nstalls_lef_lef = 0;
  this->_nstalls_lef_bar = 0;
}

bool Lef::is_bound() const {
  assert((this->rev_unit._pos == std::numeric_limits<Bp>::max()) ==
         (this->fwd_unit._pos == std::numeric_limits<Bp>::max()));
  return this->lifetime != 0;
}

void Lef::bind_at_pos(Bp pos, Bp lifetime_) {
  assert(this->lifetime == 0);

  this->rev_unit._pos = pos;
  this->fwd_unit._pos = pos;

  this->rev_unit._nstalls_lef_lef = 0;
  this->fwd_unit._nstalls_lef_lef = 0;

  this->rev_unit._nstalls_lef_bar = 0;
  this->fwd_unit._nstalls_lef_bar = 0;

  this->lifetime = lifetime_;
}

void Lef::release() {
  this->rev_unit.release();
  this->fwd_unit.release();
  this->lifetime = 0;
}

}  // namespace modle
