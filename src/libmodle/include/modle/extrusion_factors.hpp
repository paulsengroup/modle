#pragma once

#include <limits>

#include "modle/common.hpp"

namespace modle {

class ExtrusionUnit {
  friend struct Lef;
  friend struct Genome;

 public:
  inline ExtrusionUnit() = default;
  inline ExtrusionUnit(Bp pos, Bp lifetime, std::size_t rank = 0, Bp nstalls_lef_lef = 0,
                       Bp nstalls_lef_bar = 0);
  // inline ExtrusionUnit(Bp pos, Iter lifetime);
  inline void bind_at_pos(Bp pos, Bp lifetime);
  [[nodiscard]] inline Bp pos() const;
  [[nodiscard]] inline std::size_t rank() const;
  [[nodiscard]] inline Bp lifetime() const;
  inline void operator-(Bp n);
  inline void operator+(Bp n);
  [[nodiscard]] inline bool operator<(const ExtrusionUnit& other);
  [[nodiscard]] inline bool operator<(std::size_t other_pos);

  [[nodiscard]] inline static double compute_probability_of_release(std::size_t avg_unit_lifetime,
                                                                    std::size_t nactive_units);

 private:
  Bp _pos{std::numeric_limits<Bp>::max()};
  std::size_t _rank{std::numeric_limits<std::size_t>::max()};
  Bp _lifetime{0};
  Bp _nstalls_lef_lef{0};
  Bp _nstalls_lef_bar{0};

#ifdef ENABLE_TESTING
 public:
  template <typename I>
  inline void set_rank(I n) {
    this->_rank = n;
  };
#endif
};

struct Lef {
  ExtrusionUnit rev_unit{};
  ExtrusionUnit fwd_unit{};
};

}  // namespace modle

#include "../../extrusion_factors_impl.hpp"
