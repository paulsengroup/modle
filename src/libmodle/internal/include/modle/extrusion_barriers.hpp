// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span
#include <fmt/format.h>       // for format_parse_context, format_error

#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cassert>                                  // for assert

#include "modle/common/common.hpp"  // for bp_t, Direction, fwd, none
#include "modle/common/random.hpp"  // for PRNG_t
#include "modle/common/utils.hpp"   // for ndebug_defined

namespace modle {

class ExtrusionUnit;

class ExtrusionBarrier {
 public:
  constexpr ExtrusionBarrier() = default;
  constexpr ExtrusionBarrier(bp_t pos, double transition_prob_occupied_to_occupied,
                             double transition_prob_non_occupied_to_not_occupied,
                             dna::Direction motif_direction);
  constexpr ExtrusionBarrier(bp_t pos, double transition_prob_occupied_to_occupied,
                             double transition_prob_not_occupied_to_not_occupied,
                             char motif_direction);

  [[nodiscard]] constexpr bp_t pos() const noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr bp_t& pos() noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr double prob_occupied_to_occupied() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr double prob_occupied_to_not_occupied() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr double prob_not_occupied_to_not_occupied() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr double prob_not_occupied_to_occupied() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr dna::Direction blocking_direction_major() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr dna::Direction blocking_direction_minor() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr bool operator==(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr bool operator!=(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr bool operator<(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr bool operator>(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr bool operator>=(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr bool operator<=(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr static double
  compute_blocking_to_blocking_transition_probabilities_from_pblock(
      double barrier_occupancy,
      double non_blocking_to_non_blocking_transition_prob) noexcept(utils::ndebug_defined());
  [[nodiscard]] constexpr double occupancy() const noexcept(utils::ndebug_defined());

 protected:
  // clang-format off
  bp_t _pos{(std::numeric_limits<bp_t>::max)()};           // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)
  double _occupied_to_occupied_transition_prob{};          // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)
  double _non_occupied_to_not_occupied_transition_prob{};  // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)
  dna::Direction _blocking_direction{dna::none};           // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)
  // clang-format on
};

namespace CTCF {
enum State : u8f { NOT_OCCUPIED = 0, OCCUPIED = 1 };
[[nodiscard]] inline State next_state(State current_state, double occupied_self_transition_prob,
                                      double not_occupied_self_transition_prob,
                                      random::PRNG_t& rand_eng);
//! Update CTCF states for the current iteration based on the states from the previous
//! iteration.

//! \param extr_barriers
//! \param mask bitset used to store CTCF states. States will be updated inplace. Bits set to 0
//!        represent CTCFs in NOT_OCCUPIED state, while bits set to 1 represent CTCFs in
//!        OCCUPIED state
//! \param rand_eng
inline void update_states(absl::Span<const ExtrusionBarrier> extr_barriers,
                          boost::dynamic_bitset<>& mask,
                          random::PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

}  // namespace CTCF
}  // namespace modle

template <>
struct fmt::formatter<modle::ExtrusionBarrier> {
  char presentation = 's';  // s == short, f == full
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  template <typename FormatContext>
  inline auto format(const modle::ExtrusionBarrier& b, FormatContext& ctx) -> decltype(ctx.out());
};

#include "../../extrusion_barriers_impl.hpp"
