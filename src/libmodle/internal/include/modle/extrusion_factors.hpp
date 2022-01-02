// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <limits>  // for numeric_limits

#include "modle/common/common.hpp"  // for bp_t, i64
#include "modle/common/utils.hpp"   // for ndebug_defined

namespace modle {

class ExtrusionBarrier;

class ExtrusionUnit {
  friend struct Lef;
  friend class Simulation;

 public:
  ExtrusionUnit() noexcept = default;
  explicit ExtrusionUnit(bp_t pos) noexcept;
  void bind_at_pos(bp_t pos) noexcept(utils::ndebug_defined());
  [[nodiscard]] bp_t pos() const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator<(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator<(const ExtrusionBarrier& barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline bool operator<(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator>(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator>(const ExtrusionBarrier& barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline bool operator>(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator<=(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator<=(const ExtrusionBarrier& barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline bool operator<=(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator>=(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator>=(const ExtrusionBarrier& barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline bool operator>=(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator==(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator==(const ExtrusionBarrier& barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline bool operator==(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] i64 operator-(const ExtrusionUnit& other) const noexcept(utils::ndebug_defined());
  [[nodiscard]] i64 operator-(const ExtrusionBarrier& other_barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline i64 operator-(I other_pos) const noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline i64 operator+(I other_pos) const noexcept(utils::ndebug_defined());
  void release() noexcept(utils::ndebug_defined());

 private:
  bp_t _pos{(std::numeric_limits<bp_t>::max)()};
};

struct Lef {
  Lef() = default;
  [[nodiscard]] bool is_bound() const noexcept(utils::ndebug_defined());
  void bind_at_pos(usize current_epoch, bp_t pos) noexcept(utils::ndebug_defined());
  void release() noexcept(utils::ndebug_defined());
  void reset() noexcept(utils::ndebug_defined());
  [[nodiscard]] bp_t loop_size() const noexcept;

  usize binding_epoch{(std::numeric_limits<usize>::max)()};
  ExtrusionUnit rev_unit{};
  ExtrusionUnit fwd_unit{};
};

}  // namespace modle

#include "../../extrusion_factors_impl.hpp"  // IWYU pragma: export
