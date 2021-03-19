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
  inline ExtrusionUnit(Bp pos, Bp nstalls_lef_lef = 0, Bp nstalls_lef_bar = 0);
  // inline ExtrusionUnit(Bp pos, Iter lifetime);
  inline void bind_at_pos(Bp pos);
  [[nodiscard]] inline Bp pos() const;
  [[nodiscard]] inline bool stalled() const;
  inline void increment_stalls(Bp nstalls);
  inline void decrement_stalls(Bp nstalls);
  [[nodiscard]] inline Bp lef_lef_stalls() const;
  [[nodiscard]] inline Bp lef_bar_stalls() const;
  inline void operator-(Bp n);
  inline void operator+(Bp n);
  [[nodiscard]] inline bool operator<(const ExtrusionUnit& other);
  [[nodiscard]] inline bool operator<(std::size_t other_pos);
  inline void release();

 private:
  Bp _pos{std::numeric_limits<Bp>::max()};
  Bp _nstalls_lef_lef{0};
  Bp _nstalls_lef_bar{0};
};

struct Lef {
  [[nodiscard]] inline bool is_bound() const;
  inline void bind_at_pos(Bp pos, Bp lifetime_);
  inline void release();
  inline void reset();
  ExtrusionUnit rev_unit{};
  ExtrusionUnit fwd_unit{};
  Bp lifetime{0};
};

}  // namespace modle

#include "../../extrusion_factors_impl.hpp"  // IWYU pragma: keep
