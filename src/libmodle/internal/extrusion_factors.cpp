// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/extrusion_factors.hpp"

#include <cassert>  // for assert
#include <cstddef>  // IWYU pragma: keep for size_t
#include <limits>   // for numeric_limits

#include "modle/common/common.hpp"       // for bp_t
#include "modle/common/utils.hpp"        // for ndebug_defined
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier

namespace modle {

ExtrusionUnit::ExtrusionUnit(bp_t pos) noexcept : _pos(pos) {}

void ExtrusionUnit::bind_at_pos(bp_t pos) noexcept(utils::ndebug_defined()) { this->_pos = pos; }

bp_t ExtrusionUnit::pos() const noexcept(utils::ndebug_defined()) { return this->_pos; }

bool ExtrusionUnit::operator<(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined()) {
  return this->pos() < other.pos();
}
bool ExtrusionUnit::operator<(const ExtrusionBarrier& barrier) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() < barrier.pos();
}

bool ExtrusionUnit::operator>(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined()) {
  return this->pos() > other.pos();
}
bool ExtrusionUnit::operator>(const ExtrusionBarrier& barrier) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() > barrier.pos();
}

bool ExtrusionUnit::operator<=(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined()) {
  return this->pos() <= other.pos();
}
bool ExtrusionUnit::operator<=(const ExtrusionBarrier& barrier) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() <= barrier.pos();
}

bool ExtrusionUnit::operator>=(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined()) {
  return this->pos() >= other.pos();
}
bool ExtrusionUnit::operator>=(const ExtrusionBarrier& barrier) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() >= barrier.pos();
}

bool ExtrusionUnit::operator==(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined()) {
  return this->pos() == other.pos();
}
bool ExtrusionUnit::operator==(const ExtrusionBarrier& barrier) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() == barrier.pos();
}

int64_t ExtrusionUnit::operator-(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return static_cast<int64_t>(this->pos()) - static_cast<int64_t>(other.pos());
}

int64_t ExtrusionUnit::operator-(const ExtrusionBarrier& other_barrier) const
    noexcept(utils::ndebug_defined()) {
  return static_cast<int64_t>(this->pos()) - static_cast<int64_t>(other_barrier.pos());
}

void ExtrusionUnit::release() noexcept(utils::ndebug_defined()) {
  this->_pos = (std::numeric_limits<bp_t>::max)();
}

bool Lef::is_bound() const noexcept(utils::ndebug_defined()) {
  // NOLINTNEXTLINE
  assert((this->rev_unit._pos == (std::numeric_limits<bp_t>::max)()) ==
         (this->fwd_unit._pos == (std::numeric_limits<bp_t>::max)()));
  return this->binding_epoch != (std::numeric_limits<size_t>::max)();
}

void Lef::bind_at_pos(size_t current_epoch, bp_t pos) noexcept(utils::ndebug_defined()) {
  this->rev_unit.bind_at_pos(pos);
  this->fwd_unit.bind_at_pos(pos);

  this->binding_epoch = current_epoch;
}

void Lef::release() noexcept(utils::ndebug_defined()) {
  this->rev_unit.release();
  this->fwd_unit.release();
  this->binding_epoch = (std::numeric_limits<size_t>::max)();
}

void Lef::reset() noexcept(utils::ndebug_defined()) { this->release(); }

bp_t Lef::loop_size() const noexcept {
  assert(this->rev_unit.pos() <= this->fwd_unit.pos());  // NOLINT
  return this->fwd_unit.pos() - this->rev_unit.pos();
}
}  // namespace modle
