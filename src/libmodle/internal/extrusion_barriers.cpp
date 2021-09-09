#include "modle/extrusion_barriers.hpp"

#include <absl/types/span.h>  // for Span

#include <algorithm>                                // for clamp
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cassert>                                  // for assert

#include "modle/common/common.hpp"  // for bp_t, Direction
#include "modle/common/random.hpp"  // for random::generate_canonical, random::PRNG_t
#include "modle/common/utils.hpp"   // for ndebug_defined

namespace modle {
ExtrusionBarrier::ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                                   double transition_prob_non_blocking_to_non_blocking,
                                   dna::Direction motif_direction)
    : _pos(pos),
      _occupied_to_occupied_transition_prob(transition_prob_blocking_to_blocking),
      _non_occupied_to_not_occupied_transition_prob(transition_prob_non_blocking_to_non_blocking),
      _blocking_direction(motif_direction == dna::Direction::fwd ? dna::Direction::rev
                                                                 : dna::Direction::fwd) {
  // NOLINTNEXTLINE
  assert(motif_direction == dna::Direction::fwd || motif_direction == dna::Direction::rev);
  assert(transition_prob_blocking_to_blocking >= 0.0 &&  // NOLINT
         transition_prob_blocking_to_blocking <= 1.0);
  assert(transition_prob_non_blocking_to_non_blocking >= 0.0 &&  // NOLINT
         transition_prob_non_blocking_to_non_blocking <= 1.0);
}

ExtrusionBarrier::ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                                   double transition_prob_non_blocking_to_non_blocking,
                                   char motif_direction)
    : _pos(pos),
      _occupied_to_occupied_transition_prob(transition_prob_blocking_to_blocking),
      _non_occupied_to_not_occupied_transition_prob(transition_prob_non_blocking_to_non_blocking),
      _blocking_direction(motif_direction == '+' ? dna::Direction::rev : dna::Direction::fwd) {
  // NOLINTNEXTLINE
  assert(motif_direction == '+' || motif_direction == '-');
  assert(transition_prob_blocking_to_blocking >= 0.0 &&  // NOLINT
         transition_prob_blocking_to_blocking <= 1.0);
  assert(transition_prob_non_blocking_to_non_blocking >= 0.0 &&  // NOLINT
         transition_prob_non_blocking_to_non_blocking <= 1.0);
}

bp_t ExtrusionBarrier::pos() const noexcept(utils::ndebug_defined()) { return this->_pos; }
bp_t& ExtrusionBarrier::pos() noexcept(utils::ndebug_defined()) { return this->_pos; }
double ExtrusionBarrier::prob_occupied_to_occupied() const noexcept(utils::ndebug_defined()) {
  return this->_occupied_to_occupied_transition_prob;
}

double ExtrusionBarrier::prob_occupied_to_not_occupied() const noexcept(utils::ndebug_defined()) {
  return 1.0 - prob_occupied_to_occupied();
}

double ExtrusionBarrier::prob_not_occupied_to_not_occupied() const
    noexcept(utils::ndebug_defined()) {
  return this->_non_occupied_to_not_occupied_transition_prob;
}

double ExtrusionBarrier::prob_not_occupied_to_occupied() const noexcept(utils::ndebug_defined()) {
  return 1.0 - prob_not_occupied_to_not_occupied();
}

dna::Direction ExtrusionBarrier::blocking_direction_major() const
    noexcept(utils::ndebug_defined()) {
  return this->_blocking_direction;
}

dna::Direction ExtrusionBarrier::blocking_direction_minor() const
    noexcept(utils::ndebug_defined()) {
  return this->_blocking_direction == dna::fwd ? dna::rev : dna::fwd;
}

bool ExtrusionBarrier::operator==(const ExtrusionBarrier& other) const noexcept {
  return this->pos() == other.pos();
}

bool ExtrusionBarrier::operator!=(const ExtrusionBarrier& other) const noexcept {
  return !(*this == other);
}

bool ExtrusionBarrier::operator<(const ExtrusionBarrier& other) const noexcept {
  return this->pos() < other.pos();
}

bool ExtrusionBarrier::operator>(const ExtrusionBarrier& other) const noexcept {
  return this->pos() > other.pos();
}

bool ExtrusionBarrier::operator<=(const ExtrusionBarrier& other) const noexcept {
  return this->pos() <= other.pos();
}

bool ExtrusionBarrier::operator>=(const ExtrusionBarrier& other) const noexcept {
  return this->pos() >= other.pos();
}

double ExtrusionBarrier::compute_blocking_to_blocking_transition_probabilities_from_pblock(
    double probability_of_barrier_block,
    double non_blocking_to_non_blocking_transition_prob) noexcept(utils::ndebug_defined()) {
  // pno = Transition prob. from non-occupied to occupied
  // pon = Transition prob. from occupied to non-occupied
  // occ = Occupancy
  const auto pno = 1.0 - non_blocking_to_non_blocking_transition_prob;
  const auto occ = probability_of_barrier_block;
  const auto pon = (pno - (occ * pno)) / occ;
  return std::clamp(1.0 - pon, 0.0, 1.0);
}
double ExtrusionBarrier::occupancy() const noexcept(utils::ndebug_defined()) {
  // Make sure this was not default constructed
  assert(this->_blocking_direction != dna::none);  // NOLINT

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
  assert(occupied_self_transition_prob >= 0 && occupied_self_transition_prob <= 1);  // NOLINT
  assert(not_occupied_self_transition_prob >= 0 &&
         not_occupied_self_transition_prob <= 1);  // NOLINT

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
  assert(extr_barriers.size() == mask.size());  // NOLINT
  for (auto i = 0UL; i < extr_barriers.size(); ++i) {
    mask[i] = CTCF::next_state(mask[i] ? CTCF::OCCUPIED : CTCF::NOT_OCCUPIED,
                               extr_barriers[i].prob_occupied_to_occupied(),
                               extr_barriers[i].prob_not_occupied_to_not_occupied(), rand_eng);
  }
}

}  // namespace modle
