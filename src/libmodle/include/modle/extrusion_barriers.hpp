#pragma once

#include <absl/types/span.h>  // for Span

#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cstdint>                                  // for uint_fast8_t
#include <random>                                   // for uniform_real_distribution

#include "modle/common.hpp"  // for bp_t_t, Direction, PRNG_t
#include "modle/utils.hpp"   // for ndebug_defined

namespace modle {

class ExtrusionUnit;

class ExtrusionBarrier {
 public:
  inline ExtrusionBarrier() = default;
  inline ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                          double transition_prob_non_blocking_to_non_blocking,
                          dna::Direction motif_direction);
  inline ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                          double transition_prob_non_blocking_to_non_blocking,
                          char motif_direction);

  [[nodiscard]] inline bp_t pos() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline double prob_occupied_to_occupied() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline double prob_occupied_to_not_occupied() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline double prob_not_occupied_to_not_occupied() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline double prob_not_occupied_to_occupied() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline dna::Direction blocking_direction_major() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline dna::Direction blocking_direction_minor() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator<(const ExtrusionBarrier& other) const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline static double
  compute_blocking_to_blocking_transition_probabilities_from_pblock(
      double probability_of_barrier_block,
      double non_blocking_to_non_blocking_transition_prob) noexcept(utils::ndebug_defined());

 protected:
  bp_t _pos;                                             // NOLINT
  double _occupied_to_occupied_transition_prob;          // NOLINT
  double _non_occupied_to_not_occupied_transition_prob;  // NOLINT
  dna::Direction _blocking_direction;                    // NOLINT
};

namespace CTCF {
enum State : uint_fast8_t { NOT_OCCUPIED = 0, OCCUPIED = 1 };
[[nodiscard]] inline State next_state(State current_state, double occupied_self_transition_prob,
                                      double not_occupied_self_transition_prob, PRNG_t& rand_eng);
using state_gen_t = std::uniform_real_distribution<double>;

//! Update CTCF states for the current iteration based on the states from the previous
//! iteration.

//! \param extr_barriers
//! \param mask bitset used to store CTCF states. States will be updated inplace. Bits set to 0
//!        represent CTCFs in NOT_OCCUPIED state, while bits set to 1 represent CTCFs in
//!        OCCUPIED state
//! \param rand_eng
inline void update_states(absl::Span<const ExtrusionBarrier> extr_barriers,
                          boost::dynamic_bitset<>& mask,
                          modle::PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

}  // namespace CTCF
}  // namespace modle

#include "../../extrusion_barriers_impl.hpp"  // IWYU pragma: keep
