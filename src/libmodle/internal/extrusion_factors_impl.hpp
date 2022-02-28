// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "modle/extrusion_factors.hpp"

#include "modle/common/common.hpp"  // for bp_t, i64
#include "modle/common/utils.hpp"   // for ndebug_defined

namespace modle {

template <typename I>
bool ExtrusionUnit::operator<(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() < other_pos;
}

template <typename I>
bool ExtrusionUnit::operator>(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() > other_pos;
}

template <typename I>
bool ExtrusionUnit::operator<=(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() <= other_pos;
}

template <typename I>
bool ExtrusionUnit::operator>=(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() >= other_pos;
}

template <typename I>
bool ExtrusionUnit::operator==(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return this->pos() == other_pos;
}

template <typename I>
i64 ExtrusionUnit::operator-(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return static_cast<i64>(this->pos()) - static_cast<i64>(other_pos);
}

template <typename I>
i64 ExtrusionUnit::operator+(I other_pos) const noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "I should be an integral number");
  return static_cast<i64>(this->pos() + other_pos);
}
}  // namespace modle
