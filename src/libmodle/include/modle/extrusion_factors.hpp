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
  inline ExtrusionUnit(bp_t pos);
  inline void bind_at_pos(bp_t pos);
  [[nodiscard]] inline bp_t pos() const;
  [[nodiscard]] inline bool operator<(const ExtrusionUnit& other);
  [[nodiscard]] inline bool operator<(std::size_t other_pos);
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
