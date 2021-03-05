#pragma once

#include "modle/common.hpp"

namespace modle {
class ExtrusionBarrier {
 public:
  inline ExtrusionBarrier() = default;
  inline ExtrusionBarrier(Bp pos, double prob_of_block, dna::Direction motif_direction);
  inline ExtrusionBarrier(Bp pos, double prob_of_block, char motif_direction);

  [[nodiscard]] inline Bp pos() const;
  [[nodiscard]] inline double pblock() const;
  [[nodiscard]] inline dna::Direction blocking_direction() const;

 private:
  Bp _pos;
  double _prob_of_block;
  dna::Direction _blocking_direction;
};
}  // namespace modle

#include "../../extrusion_barriers_impl.hpp"
