// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <limits>

#include "modle/common/common.hpp"
#include "modle/common/dna.hpp"
#include "modle/common/random.hpp"

namespace modle {

template <class BarrierIt, class StateIt>
ExtrusionBarriers::ExtrusionBarriers(BarrierIt first_barrier, BarrierIt last_barrier,
                                     StateIt first_state, bool sort_barriers)
    : ExtrusionBarriers(static_cast<usize>(std::distance(first_barrier, last_barrier))) {
  clear();

  auto it1 = first_barrier;
  auto it2 = first_state;

  while (it1 != last_barrier) {
    push_back(*it1++, *it2++);
  }

  if (sort_barriers) {
    sort();
  }
}

template <class BarrierIt>
ExtrusionBarriers::ExtrusionBarriers(BarrierIt first_barrier, BarrierIt last_barrier, State state,
                                     bool sort_barriers)
    : ExtrusionBarriers(static_cast<usize>(std::distance(first_barrier, last_barrier))) {
  clear();

  std::for_each(first_barrier, last_barrier,
                [&](const auto& barrier) { push_back(barrier, state); });

  if (sort_barriers) {
    sort();
  }
}

constexpr ExtrusionBarrier::ExtrusionBarrier(bp_t pos_, TP transition_prob_active_to_active,
                                             TP transition_prob_inactive_to_inactive,
                                             dna::Direction motif_direction)
    : pos(pos_),
      stp_active(transition_prob_active_to_active),
      stp_inactive(transition_prob_inactive_to_inactive),
      blocking_direction(motif_direction == dna::FWD ? dna::REV : dna::FWD) {
  assert(motif_direction == dna::FWD || motif_direction == dna::REV);
  assert(transition_prob_active_to_active() >= 0.0 && transition_prob_active_to_active() <= 1.0);
  assert(transition_prob_inactive_to_inactive() >= 0.0 &&
         transition_prob_inactive_to_inactive() <= 1.0);
}

constexpr ExtrusionBarrier::ExtrusionBarrier(bp_t pos_, TP transition_prob_active_to_active,
                                             TP transition_prob_inactive_to_inactive,
                                             char motif_direction)
    : pos(pos_),
      stp_active(transition_prob_active_to_active),
      stp_inactive(transition_prob_inactive_to_inactive),
      blocking_direction(motif_direction == '+' ? dna::REV : dna::FWD) {
  assert(motif_direction == '+' || motif_direction == '-');
  assert(transition_prob_active_to_active() >= 0.0 && transition_prob_active_to_active() <= 1.0);
  assert(transition_prob_inactive_to_inactive() >= 0.0 &&
         transition_prob_inactive_to_inactive() <= 1.0);
}

constexpr bool ExtrusionBarrier::operator==(const ExtrusionBarrier& other) const noexcept {
  return pos == other.pos;
}

constexpr bool ExtrusionBarrier::operator!=(const ExtrusionBarrier& other) const noexcept {
  return !(*this == other);
}

constexpr bool ExtrusionBarrier::operator<(const ExtrusionBarrier& other) const noexcept {
  return pos < other.pos;
}

constexpr bool ExtrusionBarrier::operator>(const ExtrusionBarrier& other) const noexcept {
  return pos > other.pos;
}

constexpr bool ExtrusionBarrier::operator<=(const ExtrusionBarrier& other) const noexcept {
  return pos <= other.pos;
}

constexpr bool ExtrusionBarrier::operator>=(const ExtrusionBarrier& other) const noexcept {
  return pos >= other.pos;
}

// NOLINTNEXTLINE(hicpp-explicit-conversions)
constexpr internal::TransitionProbability::TransitionProbability(double p) noexcept : _p(p) {
  assert(_p >= 0.0);
  assert(_p <= 1.0);
}

constexpr double internal::TransitionProbability::operator()() const noexcept { return _p; }

constexpr double ExtrusionBarrier::compute_stp_active_from_occupancy(TP stp_inactive,
                                                                     double occupancy) noexcept {
  if (MODLE_UNLIKELY(occupancy == 0)) {
    return 0.0;
  }

  const auto tp_inactive_to_active = 1.0 - stp_inactive();
  const auto tp_active_to_inactive =
      (tp_inactive_to_active - (occupancy * tp_inactive_to_active)) / occupancy;
  return std::clamp(1.0 - tp_active_to_inactive, 0.0, 1.0);
}

constexpr double ExtrusionBarrier::compute_occupancy_from_stp(TP stp_active,
                                                              TP stp_inactive) noexcept {
  if (MODLE_UNLIKELY(stp_active() + stp_inactive() == 0)) {
    return 0.0;
  }

  const auto tp_inactive_to_active = 1.0 - stp_inactive();
  const auto tp_active_to_inactive = 1.0 - stp_active();
  const auto occupancy = tp_inactive_to_active / (tp_inactive_to_active + tp_active_to_inactive);
  return std::clamp(occupancy, 0.0, 1.0);
}

}  // namespace modle

constexpr auto fmt::formatter<modle::ExtrusionBarrier>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  const auto* it = ctx.begin();
  const auto* end = ctx.end();
  if (it != end && (*it == 's' || *it == 'f')) {
    presentation = *it++;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  }

  // Check if reached the end of the range:
  if (it != end && *it != '}') {
    throw fmt::format_error("invalid format");
  }

  // Return an iterator past the end of the parsed range:
  return it;
}

// Formats the point p using the parsed format specification (presentation)
// stored in this formatter.
template <typename FormatContext>
inline auto fmt::formatter<modle::ExtrusionBarrier>::format(const modle::ExtrusionBarrier& b,
                                                            FormatContext& ctx) const
    -> decltype(ctx.out()) {
  // ctx.out() is an output iterator to write to.
  if (presentation == 's') {
    return fmt::format_to(ctx.out(), "ExtrusionBarrier{{pos={}; motif_direction={}}}", b.pos,
                          b.blocking_direction == modle::dna::FWD ? "rev" : "fwd");
  }
  assert(presentation == 'f');
  return fmt::format_to(
      ctx.out(), "ExtrusionBarrier{{pos={}; motif_direction={}; Pbb={:.4f}; Puu={:.4f}}}", b.pos,
      b.blocking_direction == modle::dna::FWD ? "rev" : "fwd", b.stp_active, b.stp_inactive);
}
