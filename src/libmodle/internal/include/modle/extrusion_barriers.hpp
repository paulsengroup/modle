#pragma once

#include <absl/types/span.h>  // for Span

#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cstdint>                                  // for uint_fast8_t

#include "modle/common/common.hpp"  // for bp_t, Direction
#include "modle/common/random.hpp"  // for random::PRNG_t
#include "modle/common/utils.hpp"   // for ndebug_defined

namespace modle {

class ExtrusionUnit;

class ExtrusionBarrier {
 public:
  ExtrusionBarrier() = default;
  ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                   double transition_prob_non_blocking_to_non_blocking,
                   dna::Direction motif_direction);
  ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                   double transition_prob_non_blocking_to_non_blocking, char motif_direction);

  [[nodiscard]] bp_t pos() const noexcept(utils::ndebug_defined());
  [[nodiscard]] bp_t& pos() noexcept(utils::ndebug_defined());
  [[nodiscard]] double prob_occupied_to_occupied() const noexcept(utils::ndebug_defined());
  [[nodiscard]] double prob_occupied_to_not_occupied() const noexcept(utils::ndebug_defined());
  [[nodiscard]] double prob_not_occupied_to_not_occupied() const noexcept(utils::ndebug_defined());
  [[nodiscard]] double prob_not_occupied_to_occupied() const noexcept(utils::ndebug_defined());
  [[nodiscard]] dna::Direction blocking_direction_major() const noexcept(utils::ndebug_defined());
  [[nodiscard]] dna::Direction blocking_direction_minor() const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator==(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] bool operator!=(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] bool operator<(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] bool operator>(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] bool operator>=(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] bool operator<=(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] static double compute_blocking_to_blocking_transition_probabilities_from_pblock(
      double probability_of_barrier_block,
      double non_blocking_to_non_blocking_transition_prob) noexcept(utils::ndebug_defined());
  [[nodiscard]] double occupancy() const noexcept(utils::ndebug_defined());

 protected:
  bp_t _pos;                                             // NOLINT
  double _occupied_to_occupied_transition_prob;          // NOLINT
  double _non_occupied_to_not_occupied_transition_prob;  // NOLINT
  dna::Direction _blocking_direction{dna::none};         // NOLINT
};

namespace CTCF {
enum State : uint_fast8_t { NOT_OCCUPIED = 0, OCCUPIED = 1 };
[[nodiscard]] State next_state(State current_state, double occupied_self_transition_prob,
                               double not_occupied_self_transition_prob, random::PRNG_t& rand_eng);
//! Update CTCF states for the current iteration based on the states from the previous
//! iteration.

//! \param extr_barriers
//! \param mask bitset used to store CTCF states. States will be updated inplace. Bits set to 0
//!        represent CTCFs in NOT_OCCUPIED state, while bits set to 1 represent CTCFs in
//!        OCCUPIED state
//! \param rand_eng
void update_states(absl::Span<const ExtrusionBarrier> extr_barriers, boost::dynamic_bitset<>& mask,
                   random::PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

}  // namespace CTCF
}  // namespace modle
