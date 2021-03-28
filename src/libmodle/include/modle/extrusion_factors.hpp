#pragma once

#include <cstddef>  // IWYU pragma: keep for size_t
#include <limits>   // for numeric_limits

#include "modle/common.hpp"  // for Bp

namespace modle {

class ExtrusionBarrier;

class ExtrusionUnit {
  friend struct Lef;
  friend class Simulation;

 public:
  inline ExtrusionUnit() = default;
  inline ExtrusionUnit(bp_t pos);
  inline void bind_at_pos(bp_t pos);
  [[nodiscard]] inline bp_t pos() const;
  [[nodiscard]] inline bool operator<(const ExtrusionUnit& other) const;
  [[nodiscard]] inline bool operator<(const ExtrusionBarrier& barrier) const;
  template <typename I>
  [[nodiscard]] inline bool operator<(I other_pos) const;
  [[nodiscard]] inline bool operator>(const ExtrusionUnit& other) const;
  [[nodiscard]] inline bool operator>(const ExtrusionBarrier& barrier) const;
  template <typename I>
  [[nodiscard]] inline bool operator>(I other_pos) const;
  [[nodiscard]] inline bool operator<=(const ExtrusionUnit& other) const;
  [[nodiscard]] inline bool operator<=(const ExtrusionBarrier& barrier) const;
  template <typename I>
  [[nodiscard]] inline bool operator<=(I other_pos) const;
  [[nodiscard]] inline bool operator>=(const ExtrusionUnit& other) const;
  [[nodiscard]] inline bool operator>=(const ExtrusionBarrier& barrier) const;
  template <typename I>
  [[nodiscard]] inline bool operator>=(I other_pos) const;
  [[nodiscard]] inline bool operator==(const ExtrusionUnit& other) const;
  [[nodiscard]] inline bool operator==(const ExtrusionBarrier& barrier) const;
  template <typename I>
  [[nodiscard]] inline bool operator==(I other_pos) const;
  [[nodiscard]] inline int64_t operator-(const ExtrusionUnit& other) const;
  [[nodiscard]] inline int64_t operator-(const ExtrusionBarrier& barrier) const;
  template <typename I>
  [[nodiscard]] inline int64_t operator-(I other_pos) const;
  template <typename I>
  [[nodiscard]] inline int64_t operator+(I other_pos) const;
  inline void release();

 private:
  bp_t _pos{std::numeric_limits<bp_t>::max()};
};

struct Lef {
  [[nodiscard]] inline bool is_bound() const;
  inline void bind_at_pos(bp_t pos);
  inline void release();
  inline void reset();
  ExtrusionUnit rev_unit{};
  ExtrusionUnit fwd_unit{};
  bool active{false};
};

}  // namespace modle

#include "../../extrusion_factors_impl.hpp"  // IWYU pragma: keep
