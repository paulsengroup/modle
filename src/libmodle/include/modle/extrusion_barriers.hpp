#pragma once

#include "modle/common.hpp"

namespace modle {
class ExtrusionBarrier {
 public:
  inline ExtrusionBarrier() = default;
  inline ExtrusionBarrier(bp_t pos, double prob_of_block, dna::Direction motif_direction);
  inline ExtrusionBarrier(bp_t pos, double prob_of_block, char motif_direction);

  [[nodiscard]] inline bp_t pos() const;
  [[nodiscard]] inline double pblock() const;
  [[nodiscard]] inline dna::Direction blocking_direction() const;

 private:
  bp_t _pos;
  double _prob_of_block;
  dna::Direction _blocking_direction;
};
}  // namespace modle

#include "../../extrusion_barriers_impl.hpp"  // IWYU pragma: keep
