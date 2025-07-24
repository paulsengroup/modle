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
  return pos() < other_pos;
}

template <typename I>
constexpr bool ExtrusionUnit::operator>(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return pos() > other_pos;
}

template <typename I>
constexpr bool ExtrusionUnit::operator<=(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return pos() <= other_pos;
}

template <typename I>
constexpr bool ExtrusionUnit::operator>=(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return pos() >= other_pos;
}

template <typename I>
constexpr bool ExtrusionUnit::operator==(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return pos() == other_pos;
}

template <typename I>
constexpr std::int64_t ExtrusionUnit::operator-(I other_pos) const
    noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return static_cast<std::int64_t>(pos()) - static_cast<std::int64_t>(other_pos);
}

template <typename I>
constexpr std::int64_t ExtrusionUnit::operator+(I other_pos) const
    noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return static_cast<std::int64_t>(pos() + other_pos);
}

constexpr ExtrusionUnit::ExtrusionUnit(bp_t pos) noexcept : _pos(pos) {}

constexpr void ExtrusionUnit::bind_at_pos(bp_t pos) noexcept(utils::ndebug_defined()) {
  _pos = pos;
}

constexpr bp_t ExtrusionUnit::pos() const noexcept(utils::ndebug_defined()) { return _pos; }

constexpr bool ExtrusionUnit::operator<(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return pos() < other.pos();
}

constexpr bool ExtrusionUnit::operator>(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return pos() > other.pos();
}

constexpr bool ExtrusionUnit::operator<=(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return pos() <= other.pos();
}

constexpr bool ExtrusionUnit::operator>=(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return pos() >= other.pos();
}

constexpr bool ExtrusionUnit::operator==(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return pos() == other.pos();
}

constexpr std::int64_t ExtrusionUnit::operator-(const ExtrusionUnit& other) const
    noexcept(utils::ndebug_defined()) {
  return static_cast<std::int64_t>(pos()) - static_cast<std::int64_t>(other.pos());
}

constexpr void ExtrusionUnit::release() noexcept(utils::ndebug_defined()) {
  _pos = (std::numeric_limits<bp_t>::max)();
}

inline Lef::Lef(std::size_t binding_epoch_, ExtrusionUnit rev_unit_,
                ExtrusionUnit fwd_unit_) noexcept
    : binding_epoch(binding_epoch_),
      rev_unit(std::move(rev_unit_)),
      fwd_unit(std::move(fwd_unit_)) {}

constexpr bool Lef::is_bound() const noexcept(utils::ndebug_defined()) {
  assert((rev_unit._pos == (std::numeric_limits<bp_t>::max)()) ==
         (fwd_unit._pos == (std::numeric_limits<bp_t>::max)()));
  return binding_epoch != (std::numeric_limits<std::size_t>::max)();
}

constexpr void Lef::bind_at_pos(std::size_t current_epoch,
                                bp_t pos) noexcept(utils::ndebug_defined()) {
  rev_unit.bind_at_pos(pos);
  fwd_unit.bind_at_pos(pos);

  binding_epoch = current_epoch;
}

constexpr void Lef::release() noexcept(utils::ndebug_defined()) {
  rev_unit.release();
  fwd_unit.release();
  binding_epoch = (std::numeric_limits<std::size_t>::max)();
}

constexpr void Lef::reset() noexcept(utils::ndebug_defined()) { release(); }

constexpr bp_t Lef::loop_size() const noexcept {
  assert(rev_unit.pos() <= fwd_unit.pos());
  return fwd_unit.pos() - rev_unit.pos();
}

}  // namespace modle
