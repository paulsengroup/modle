// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "modle/extrusion_factors.hpp"

#include "modle/common/common.hpp"
#include "modle/common/utils.hpp"

namespace modle {

template <typename I>
constexpr bool ExtrusionUnit::operator<(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() < other_pos;
}

template <typename I>
constexpr bool ExtrusionUnit::operator>(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() > other_pos;
}

template <typename I>
constexpr bool ExtrusionUnit::operator<=(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() <= other_pos;
}

template <typename I>
constexpr bool ExtrusionUnit::operator>=(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() >= other_pos;
}

template <typename I>
constexpr bool ExtrusionUnit::operator==(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() == other_pos;
}

template <typename I>
constexpr i64 ExtrusionUnit::operator-(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return static_cast<i64>(this->pos()) - static_cast<i64>(other_pos);
}

template <typename I>
constexpr i64 ExtrusionUnit::operator+(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return static_cast<i64>(this->pos() + other_pos);
}

constexpr ExtrusionUnit::ExtrusionUnit(bp_t pos) noexcept : _pos(pos) {}

constexpr void ExtrusionUnit::bind_at_pos(bp_t pos) noexcept(utils::ndebug_defined()) {
  this->_pos = pos;
}

constexpr bp_t ExtrusionUnit::pos() const noexcept(utils::ndebug_defined()) { return this->_pos; }

constexpr bool ExtrusionUnit::operator<(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() < other.pos();
}

constexpr bool ExtrusionUnit::operator>(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() > other.pos();
}

constexpr bool ExtrusionUnit::operator<=(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() <= other.pos();
}

constexpr bool ExtrusionUnit::operator>=(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() >= other.pos();
}

constexpr bool ExtrusionUnit::operator==(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return this->pos() == other.pos();
}

constexpr i64 ExtrusionUnit::operator-(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return static_cast<i64>(this->pos()) - static_cast<i64>(other.pos());
}

constexpr void ExtrusionUnit::release() noexcept(utils::ndebug_defined()) {
  this->_pos = (std::numeric_limits<bp_t>::max)();
}

inline Lef::Lef(usize binding_epoch_, ExtrusionUnit rev_unit_, ExtrusionUnit fwd_unit_) noexcept
    : binding_epoch(binding_epoch_),
      rev_unit(std::move(rev_unit_)),
      fwd_unit(std::move(fwd_unit_)) {}

constexpr bool Lef::is_bound() const noexcept(utils::ndebug_defined()) {
  assert((this->rev_unit._pos == (std::numeric_limits<bp_t>::max)()) ==
         (this->fwd_unit._pos == (std::numeric_limits<bp_t>::max)()));
  return this->binding_epoch != (std::numeric_limits<usize>::max)();
}

constexpr void Lef::bind_at_pos(usize current_epoch, bp_t pos) noexcept(utils::ndebug_defined()) {
  this->rev_unit.bind_at_pos(pos);
  this->fwd_unit.bind_at_pos(pos);

  this->binding_epoch = current_epoch;
}

constexpr void Lef::release() noexcept(utils::ndebug_defined()) {
  this->rev_unit.release();
  this->fwd_unit.release();
  this->binding_epoch = (std::numeric_limits<usize>::max)();
}

constexpr void Lef::reset() noexcept(utils::ndebug_defined()) { this->release(); }

constexpr bp_t Lef::loop_size() const noexcept {
  assert(this->rev_unit.pos() <= this->fwd_unit.pos());
  return this->fwd_unit.pos() - this->rev_unit.pos();
}

}  // namespace modle
