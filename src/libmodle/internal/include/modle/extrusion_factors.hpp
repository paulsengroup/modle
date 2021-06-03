#pragma once

#include <cstddef>  // IWYU pragma: keep for size_t
#include <cstdint>  // for int64_t
#include <limits>   // for numeric_limits

#include "modle/common.hpp"  // for bp_t_t
#include "modle/utils.hpp"   // for ndebug_defined

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
  [[nodiscard]] int64_t operator-(const ExtrusionUnit& other) const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] int64_t operator-(const ExtrusionBarrier& other_barrier) const
      noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline int64_t operator-(I other_pos) const noexcept(utils::ndebug_defined());
  template <typename I>
  [[nodiscard]] inline int64_t operator+(I other_pos) const noexcept(utils::ndebug_defined());
  void release() noexcept(utils::ndebug_defined());

 private:
  bp_t _pos{std::numeric_limits<bp_t>::max()};
};

struct Lef {
  [[nodiscard]] bool is_bound() const noexcept(utils::ndebug_defined());
  void bind_at_pos(bp_t pos) noexcept(utils::ndebug_defined());
  void release() noexcept(utils::ndebug_defined());
  void reset() noexcept(utils::ndebug_defined());
  ExtrusionUnit rev_unit{};
  ExtrusionUnit fwd_unit{};
  bool active{false};
};

}  // namespace modle

#include "../../extrusion_factors_impl.hpp"  // IWYU pragma: keep
