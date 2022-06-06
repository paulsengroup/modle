// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span

#include <algorithm>                                // for clamp
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset, dynamic_bitset<>::ref...
#include <cassert>                                  // for assert
#include <limits>                                   // for numeric_limits, numeric_limits<>::digits

#include "modle/common/common.hpp"  // for Direction, fwd, bp_t, rev, none
#include "modle/common/random.hpp"  // for generate_canonical, PRNG_t
#include "modle/common/utils.hpp"   // for ndebug_defined

namespace modle {
constexpr ExtrusionBarrier::ExtrusionBarrier(bp_t pos, double transition_prob_occupied_to_occupied,
                                             double transition_prob_non_occupied_to_not_occupied,
                                             dna::Direction motif_direction)
    : _pos(pos),
      _occupied_to_occupied_transition_prob(transition_prob_occupied_to_occupied),
      _non_occupied_to_not_occupied_transition_prob(transition_prob_non_occupied_to_not_occupied),
      _blocking_direction(motif_direction == dna::Direction::fwd ? dna::Direction::rev
                                                                 : dna::Direction::fwd) {
  assert(motif_direction == dna::Direction::fwd || motif_direction == dna::Direction::rev);
  assert(transition_prob_occupied_to_occupied >= 0.0 &&
         transition_prob_occupied_to_occupied <= 1.0);
  assert(transition_prob_non_occupied_to_not_occupied >= 0.0 &&
         transition_prob_non_occupied_to_not_occupied <= 1.0);
}

constexpr ExtrusionBarrier::ExtrusionBarrier(bp_t pos, double transition_prob_occupied_to_occupied,
                                             double transition_prob_not_occupied_to_not_occupied,
                                             char motif_direction)
    : _pos(pos),
      _occupied_to_occupied_transition_prob(transition_prob_occupied_to_occupied),
      _non_occupied_to_not_occupied_transition_prob(transition_prob_not_occupied_to_not_occupied),
      _blocking_direction(motif_direction == '+' ? dna::Direction::rev : dna::Direction::fwd) {
  assert(motif_direction == '+' || motif_direction == '-');
  assert(transition_prob_occupied_to_occupied >= 0.0 &&
         transition_prob_occupied_to_occupied <= 1.0);
  assert(transition_prob_not_occupied_to_not_occupied >= 0.0 &&
         transition_prob_not_occupied_to_not_occupied <= 1.0);
}

constexpr bp_t ExtrusionBarrier::pos() const noexcept(utils::ndebug_defined()) {
  return this->_pos;
}
constexpr bp_t& ExtrusionBarrier::pos() noexcept(utils::ndebug_defined()) { return this->_pos; }
constexpr double ExtrusionBarrier::prob_occupied_to_occupied() const
    noexcept(utils::ndebug_defined()) {
  return this->_occupied_to_occupied_transition_prob;
}

constexpr double ExtrusionBarrier::prob_occupied_to_not_occupied() const
    noexcept(utils::ndebug_defined()) {
  return 1.0 - prob_occupied_to_occupied();
}

constexpr double ExtrusionBarrier::prob_not_occupied_to_not_occupied() const
    noexcept(utils::ndebug_defined()) {
  return this->_non_occupied_to_not_occupied_transition_prob;
}

constexpr double ExtrusionBarrier::prob_not_occupied_to_occupied() const
    noexcept(utils::ndebug_defined()) {
  return 1.0 - prob_not_occupied_to_not_occupied();
}

constexpr dna::Direction ExtrusionBarrier::blocking_direction_major() const
    noexcept(utils::ndebug_defined()) {
  return this->_blocking_direction;
}

constexpr dna::Direction ExtrusionBarrier::blocking_direction_minor() const
    noexcept(utils::ndebug_defined()) {
  return this->_blocking_direction == dna::fwd ? dna::rev : dna::fwd;
}

constexpr bool ExtrusionBarrier::operator==(const ExtrusionBarrier& other) const noexcept {
  return this->pos() == other.pos();
}

constexpr bool ExtrusionBarrier::operator!=(const ExtrusionBarrier& other) const noexcept {
  return !(*this == other);
}

constexpr bool ExtrusionBarrier::operator<(const ExtrusionBarrier& other) const noexcept {
  return this->pos() < other.pos();
}

constexpr bool ExtrusionBarrier::operator>(const ExtrusionBarrier& other) const noexcept {
  return this->pos() > other.pos();
}

constexpr bool ExtrusionBarrier::operator<=(const ExtrusionBarrier& other) const noexcept {
  return this->pos() <= other.pos();
}

constexpr bool ExtrusionBarrier::operator>=(const ExtrusionBarrier& other) const noexcept {
  return this->pos() >= other.pos();
}

constexpr double
ExtrusionBarrier::compute_blocking_to_blocking_transition_probabilities_from_pblock(
    double barrier_occupancy,
    double non_blocking_to_non_blocking_transition_prob) noexcept(utils::ndebug_defined()) {
  // pno = Transition prob. from non-occupied to occupied
  // pon = Transition prob. from occupied to non-occupied
  // occ = Occupancy
  const auto pno = 1.0 - non_blocking_to_non_blocking_transition_prob;
  const auto occ = barrier_occupancy;
  const auto pon = (pno - (occ * pno)) / occ;
  return std::clamp(1.0 - pon, 0.0, 1.0);
}
constexpr double ExtrusionBarrier::occupancy() const noexcept(utils::ndebug_defined()) {
  // Make sure this was not default constructed
  assert(this->_blocking_direction != dna::none);

  // pno = Transition prob. from non-occupied to occupied
  // pon = Transition prob. from occupied to non-occupied
  // occ = Occupancy
  const auto pno = this->prob_not_occupied_to_occupied();
  const auto pon = this->prob_occupied_to_not_occupied();
  const auto occ = pno / (pon + pno);
  return std::clamp(occ, 0.0, 1.0);
}

CTCF::State CTCF::next_state(CTCF::State current_state, double occupied_self_transition_prob,
                             double not_occupied_self_transition_prob, random::PRNG_t& rand_eng) {
  assert(occupied_self_transition_prob >= 0 && occupied_self_transition_prob <= 1);
  assert(not_occupied_self_transition_prob >= 0 && not_occupied_self_transition_prob <= 1);

  const auto p = random::generate_canonical<double, std::numeric_limits<double>::digits>(rand_eng);
  if (current_state == NOT_OCCUPIED && p > not_occupied_self_transition_prob) {
    return OCCUPIED;
  }
  if (current_state == OCCUPIED && p > occupied_self_transition_prob) {
    return NOT_OCCUPIED;
  }

  return current_state;
}

void CTCF::update_states(absl::Span<const ExtrusionBarrier> extr_barriers,
                         boost::dynamic_bitset<>& mask,
                         random::PRNG_t& rand_eng) noexcept(utils::ndebug_defined()) {
  assert(extr_barriers.size() == mask.size());
  for (usize i = 0; i < extr_barriers.size(); ++i) {
    // NOLINTNEXTLINE(readability-implicit-bool-conversion)
    mask[i] = CTCF::next_state(mask[i] ? CTCF::OCCUPIED : CTCF::NOT_OCCUPIED,
                               extr_barriers[i].prob_occupied_to_occupied(),
                               extr_barriers[i].prob_not_occupied_to_not_occupied(), rand_eng);
  }
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
                                                            FormatContext& ctx)
    -> decltype(ctx.out()) {
  // ctx.out() is an output iterator to write to.
  if (presentation == 's') {
    return fmt::format_to(ctx.out(), "ExtrusionBarrier{{pos={}; motif_direction={}}}", b.pos(),
                          b.blocking_direction_major() == modle::dna::fwd ? "rev" : "fwd");
  }
  assert(presentation == 'f');
  return fmt::format_to(ctx.out(),
                        "ExtrusionBarrier{{pos={}; motif_direction={}; Pbb={:.4f}; Puu={:.4f}}}",
                        b.pos(), b.blocking_direction_major() == modle::dna::fwd ? "rev" : "fwd",
                        b.prob_occupied_to_occupied(), b.prob_not_occupied_to_not_occupied());
}