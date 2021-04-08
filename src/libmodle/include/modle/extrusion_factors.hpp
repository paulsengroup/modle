#pragma once

#include <cstddef>  // IWYU pragma: keep for size_t
#include <limits>   // for numeric_limits

#include "modle/common.hpp"  // for bp_t_t
#include "modle/utils.hpp"   // for ndebug_defined

namespace modle {

class ExtrusionBarrier;

class ExtrusionUnit {
  friend struct Lef;
  friend class Simulation;

 public:
  inline ExtrusionUnit() noexcept(utils::ndebug_defined()) = default;
  inline explicit ExtrusionUnit(bp_t pos) noexcept(utils::ndebug_defined());
  inline void bind_at_pos(bp_t pos) noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bp_t pos() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator<(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator<(const ExtrusionBarrier& barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline bool operator<(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator>(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator>(const ExtrusionBarrier& barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline bool operator>(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator<=(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator<=(const ExtrusionBarrier& barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline bool operator<=(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator>=(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator>=(const ExtrusionBarrier& barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline bool operator>=(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator==(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator==(const ExtrusionBarrier& barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline bool operator==(I other_pos) const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline int64_t operator-(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline int64_t operator-(const ExtrusionBarrier& other_barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline int64_t operator-(I other_pos) const noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline int64_t operator+(I other_pos) const noexcept(utils::ndebug_defined());
  inline void release() noexcept(utils::ndebug_defined());

 private:
  bp_t _pos{std::numeric_limits<bp_t>::max()};
};

struct Lef {
  [[nodiscard]] inline bool is_bound() const noexcept(utils::ndebug_defined());
  inline void bind_at_pos(bp_t pos) noexcept(utils::ndebug_defined());
  inline void release() noexcept(utils::ndebug_defined());
  inline void reset() noexcept(utils::ndebug_defined());
  ExtrusionUnit rev_unit{};
  ExtrusionUnit fwd_unit{};
  bool active{false};
};

}  // namespace modle

#include "../../extrusion_factors_impl.hpp"  // IWYU pragma: keep
