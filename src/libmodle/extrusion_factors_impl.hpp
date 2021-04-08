#pragma once

// IWYU pragma: private, include "modle/extrusion_factors.hpp"

#include <algorithm>  // for min, clamp
#include <cassert>    // for assert
#include <cstddef>    // IWYU pragma: keep for size_t
#include <limits>     // for numeric_limits

#include "modle/common.hpp"              // for bp_t
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier
#include "modle/utils.hpp"               // for ndebug_defined

namespace modle {

ExtrusionUnit::ExtrusionUnit(bp_t pos) noexcept(utils::ndebug_defined()) : _pos(pos) {}

void ExtrusionUnit::bind_at_pos(bp_t pos) noexcept(utils::ndebug_defined()) { this->_pos = pos; }

bp_t ExtrusionUnit::pos() const noexcept(utils::ndebug_defined()) { return this->_pos; }

bool ExtrusionUnit::operator<(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined()) {
  return this->pos() < other.pos();
}
bool ExtrusionUnit::operator<(const ExtrusionBarrier& barrier) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() < barrier.pos();
}
template <typename I>
bool ExtrusionUnit::operator<(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() < other_pos;
}

bool ExtrusionUnit::operator>(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined()) {
  return this->pos() > other.pos();
}
bool ExtrusionUnit::operator>(const ExtrusionBarrier& barrier) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() > barrier.pos();
}
template <typename I>
bool ExtrusionUnit::operator>(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() > other_pos;
}

bool ExtrusionUnit::operator<=(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined()) {
  return this->pos() <= other.pos();
}
bool ExtrusionUnit::operator<=(const ExtrusionBarrier& barrier) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() <= barrier.pos();
}
template <typename I>
bool ExtrusionUnit::operator<=(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() <= other_pos;
}

bool ExtrusionUnit::operator>=(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined()) {
  return this->pos() >= other.pos();
}
bool ExtrusionUnit::operator>=(const ExtrusionBarrier& barrier) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() >= barrier.pos();
}
template <typename I>
bool ExtrusionUnit::operator>=(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() >= other_pos;
}

bool ExtrusionUnit::operator==(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined()) {
  return this->pos() == other.pos();
}
bool ExtrusionUnit::operator==(const ExtrusionBarrier& barrier) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() == barrier.pos();
}
template <typename I>
bool ExtrusionUnit::operator==(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() == other_pos;
}

int64_t ExtrusionUnit::operator-(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return static_cast<int64_t>(this->pos()) - static_cast<int64_t>(other.pos());
}

int64_t ExtrusionUnit::operator-(const ExtrusionBarrier& other_barrier) const
    noexcept(utils::ndebug_defined()) {
  return static_cast<int64_t>(this->pos()) - static_cast<int64_t>(other_barrier.pos());
}

template <typename I>
int64_t ExtrusionUnit::operator-(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return static_cast<int64_t>(this->pos()) - static_cast<int64_t>(other_pos);
}

template <typename I>
int64_t ExtrusionUnit::operator+(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return static_cast<int64_t>(this->pos() + other_pos);
}

void ExtrusionUnit::release() noexcept(utils::ndebug_defined()) {
  this->_pos = std::numeric_limits<bp_t>::max();
}

bool Lef::is_bound() const noexcept(utils::ndebug_defined()) {
  assert((this->rev_unit._pos == std::numeric_limits<bp_t>::max()) ==
         (this->fwd_unit._pos == std::numeric_limits<bp_t>::max()));
  return this->active;
}

void Lef::bind_at_pos(bp_t pos) noexcept(utils::ndebug_defined()) {
  this->rev_unit.bind_at_pos(pos);
  this->fwd_unit.bind_at_pos(pos);

  this->active = true;
}

void Lef::release() noexcept(utils::ndebug_defined()) {
  this->rev_unit.release();
  this->fwd_unit.release();
  this->active = false;
}

void Lef::reset() noexcept(utils::ndebug_defined()) { this->release(); }

}  // namespace modle
