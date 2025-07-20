// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <limits>  // for numeric_limits

#include "modle/common/common.hpp"  // for bp_t, i64
#include "modle/common/utils.hpp"   // for ndebug_defined

namespace modle {

class ExtrusionUnit {
  friend struct Lef;
  friend class Simulation;

 public:
  constexpr ExtrusionUnit() noexcept = default;
  constexpr explicit ExtrusionUnit(bp_t pos) noexcept;
  constexpr void bind_at_pos(bp_t pos) noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr bp_t pos() const noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr bool operator<(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] constexpr bool operator<(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr bool operator>(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] constexpr bool operator>(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr bool operator<=(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] constexpr bool operator<=(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr bool operator>=(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] constexpr bool operator>=(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr bool operator==(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] constexpr bool operator==(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr i64 operator-(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] constexpr i64 operator-(I other_pos) const noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] constexpr i64 operator+(I other_pos) const noexcept(utils::ndebug_defined());
  constexpr void release() noexcept(utils::ndebug_defined());

 private:
  bp_t _pos{(std::numeric_limits<bp_t>::max)()};
};

struct Lef {
  constexpr Lef() = default;
  inline Lef(usize binding_epoch_, ExtrusionUnit rev_unit_, ExtrusionUnit fwd_unit_) noexcept;
  [[nodiscard]] constexpr bool is_bound() const noexcept(utils::ndebug_defined());
  constexpr void bind_at_pos(usize current_epoch, bp_t pos) noexcept(utils::ndebug_defined());
  constexpr void release() noexcept(utils::ndebug_defined());
  constexpr void reset() noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr bp_t loop_size() const noexcept;

  usize binding_epoch{(std::numeric_limits<usize>::max)()};
  ExtrusionUnit rev_unit{};
  ExtrusionUnit fwd_unit{};
};

}  // namespace modle

#include "../../extrusion_factors_impl.hpp"  // IWYU pragma: export
