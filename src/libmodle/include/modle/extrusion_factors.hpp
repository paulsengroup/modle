#pragma once

#include <cstddef>  // IWYU pragma: keep for size_t
#include <limits>   // for numeric_limits

#include "modle/common.hpp"  // for Bp

namespace modle {

class ExtrusionUnit {
  friend struct Lef;
  friend class Simulation;

 public:
  inline ExtrusionUnit() = default;
  inline ExtrusionUnit(bp_t pos, bp_t nstalls_lef_lef = 0, bp_t nstalls_lef_bar = 0);
  // inline ExtrusionUnit(Bp pos, Iter lifetime);
  inline void bind_at_pos(bp_t pos);
  [[nodiscard]] inline bp_t pos() const;
  [[nodiscard]] inline bool stalled() const;
  inline void increment_stalls(bp_t nstalls);
  inline void decrement_stalls(bp_t nstalls);
  [[nodiscard]] inline bp_t lef_lef_stalls() const;
  [[nodiscard]] inline bp_t lef_bar_stalls() const;
  inline void operator-(bp_t n);
  inline void operator+(bp_t n);
  [[nodiscard]] inline bool operator<(const ExtrusionUnit& other);
  [[nodiscard]] inline bool operator<(std::size_t other_pos);
  inline void release();

 private:
  bp_t _pos{std::numeric_limits<bp_t>::max()};
  bp_t _nstalls_lef_lef{0};
  bp_t _nstalls_lef_bar{0};
};

struct Lef {
  [[nodiscard]] inline bool is_bound() const;
  inline void bind_at_pos(bp_t pos, bp_t lifetime_);
  inline void release();
  inline void reset();
  ExtrusionUnit rev_unit{};
  ExtrusionUnit fwd_unit{};
  bp_t lifetime{0};
};

}  // namespace modle

#include "../../extrusion_factors_impl.hpp"  // IWYU pragma: keep
