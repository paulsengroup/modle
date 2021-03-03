#pragma once

#include <cstdint>

#include "modle/common.hpp"

namespace modle {

ExtrusionUnit::ExtrusionUnit(Bp pos, Bp lifetime, std::size_t rank, Bp nstalls_lef_lef,
                             Bp nstalls_lef_bar)
    : _pos(pos),
      _rank(rank),
      _lifetime(lifetime),
      _nstalls_lef_lef(nstalls_lef_bar),
      _nstalls_lef_bar(nstalls_lef_bar) {}

void ExtrusionUnit::bind_at_pos(Bp pos, Bp lifetime) {
  this->_pos = pos;
  this->_lifetime = lifetime;
  this->_nstalls_lef_lef = 0;
  this->_nstalls_lef_bar = 0;
}

Bp ExtrusionUnit::pos() const { return this->_pos; }
std::size_t ExtrusionUnit::rank() const { return this->_rank; }

Bp ExtrusionUnit::lifetime() const { return this->_lifetime; }

void ExtrusionUnit::operator-(Bp n) { this->_pos = std::clamp(this->_pos - n, 0UL, this->_pos); }

void ExtrusionUnit::operator+(Bp n) {
  this->_pos = std::clamp(this->_pos + n, this->_pos, std::numeric_limits<Bp>::max());
}

bool ExtrusionUnit::operator<(const ExtrusionUnit& other) { return this->_pos < other._pos; }
bool ExtrusionUnit::operator<(std::size_t other_pos) { return this->_pos < other_pos; }

double ExtrusionUnit::compute_probability_of_release(std::size_t avg_lef_lifetime,
                                                     std::size_t nactive_units) {
  return static_cast<double>(nactive_units) / (static_cast<double>(avg_lef_lifetime));
}

}  // namespace modle
