#pragma once

#include <absl/types/span.h>  // for Span

#include <boost/range/adaptor/reversed.hpp>  // for reversed_range, reverse
#include <cassert>                           // for assert
#include <cstddef>                           // for size_t

#include "modle/common.hpp"             // for bp_t, MODLE_UNLIKELY
#include "modle/extrusion_factors.hpp"  // for ExtrusionUnit, Lef
#include "modle/utils.hpp"              // for ndebug_defined

namespace modle {

void Simulation::correct_moves_for_lef_bar_collisions(
    const absl::Span<const Lef> lefs, const absl::Span<const ExtrusionBarrier> barriers,
    const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
    const absl::Span<const collision_t> rev_collisions,
    const absl::Span<const collision_t> fwd_collisions) noexcept(utils::ndebug_defined()) {
  // LEF-BAR collisions are encoded with a number between 0 and nbarriers
  // This number corresponds to the index of the barrier that is causing the collision.
  const auto upper_bound = barriers.size();

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (rev_collisions[i] < upper_bound) MODLE_UNLIKELY {  // Process rev collisions
        const auto& barrier_idx = rev_collisions[i];
        const auto& barrier = barriers[barrier_idx];
        assert(lefs[i].rev_unit.pos() > barrier.pos());  // NOLINT

        // Compute the distance and update the rev_move such that after extruding once, the rev unit
        // of the current LEF will be located 1bp downstream of the extr. barrier.
        const auto distance = lefs[i].rev_unit.pos() - barrier.pos();
        rev_moves[i] = distance > 1UL ? distance - 1 : 0UL;
      }

    if (fwd_collisions[i] < upper_bound) MODLE_UNLIKELY {  // Process fwd collisions
        const auto& barrier_idx = fwd_collisions[i];
        const auto& barrier = barriers[barrier_idx];
        assert(lefs[i].fwd_unit.pos() < barrier.pos());  // NOLINT

        // Same as above. In this case the unit will be located 1bp upstream of the extr. barrier.
        const auto distance = barrier.pos() - lefs[i].fwd_unit.pos();
        fwd_moves[i] = distance > 1UL ? distance - 1 : 0UL;
      }
  }
}

inline void Simulation::correct_moves_for_primary_lef_lef_collisions(
    const absl::Span<const Lef> lefs, const absl::Span<const ExtrusionBarrier> barriers,
    const absl::Span<const size_t> rev_ranks, const absl::Span<const size_t> fwd_ranks,
    const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
    const absl::Span<const collision_t> rev_collisions,
    const absl::Span<const collision_t> fwd_collisions) noexcept(utils::ndebug_defined()) {
  // Primary LEF-LEF collisions are encoded with a number between nbarriers and nbarriers + nlefs.
  // Given a pair of extr. units that are moving in opposite directions, the index i corresponding
  // to the extr. unit that is causing the collision is encoded as nbarriers + i.
  const auto lower_bound = barriers.size();
  const auto upper_bound = lower_bound + lefs.size();

  auto is_lef_lef_primary_collision = [&](const auto i) constexpr {
    return i >= lower_bound && i < upper_bound;
  };

  auto is_lef_bar_collision = [&](const auto i) constexpr { return i < lower_bound; };

  for (auto rev_idx : rev_ranks) {  // Loop over rev units in 5'-3' order
    if (auto rev_collision = rev_collisions[rev_idx]; is_lef_lef_primary_collision(rev_collision))
      MODLE_UNLIKELY {
        // Decode the index of the fwd unit involved in the current collision
        rev_collision -= lower_bound;
        const auto& fwd_idx = rev_collision;

        if (auto fwd_collision = fwd_collisions[fwd_idx];
            is_lef_lef_primary_collision(fwd_collision)) {
          // This branch handles the typical case, where a pair of extr. units moving in opposite
          // direction bumped into each other, causing a primary LEF-LEF collision.
          const auto& rev_unit = lefs[rev_idx].rev_unit;
          const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

          auto& rev_move = rev_moves[rev_idx];
          auto& fwd_move = fwd_moves[fwd_idx];

          const auto [p1, p2] =
              compute_lef_lef_collision_pos(rev_unit, fwd_unit, rev_move, fwd_move);
          assert(rev_unit.pos() >= p1);  // NOLINT
          assert(fwd_unit.pos() <= p2);  // NOLINT

          // Update the moves of the involved units, so that after the next call to extrude these
          // two units will be located nex to each others at the collision site.
          rev_move = rev_unit.pos() - p1;
          fwd_move = p2 - fwd_unit.pos();

        } else if (is_lef_bar_collision(fwd_collision)) {
          // This branch handles the special case where the fwd unit involved in the collision is
          // blocked by an extrusion barrier, and thus cannot be moved.
          // In this case we update the move such that after extrusion, the rev unit will be located
          // next to the fwd unit.
          const auto& rev_unit = lefs[rev_idx].rev_unit;
          const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

          auto& rev_move = rev_moves[rev_idx];
          const auto& fwd_move = fwd_moves[fwd_idx];

          assert(rev_unit.pos() >= fwd_unit.pos() + fwd_move);  // NOLINT
          rev_move = rev_unit.pos() - (fwd_unit.pos() + fwd_move);
        }
      }
  }

  // This second loop is required in order to properly handle the scenario where a fwd unit collides
  // with a rev unit that is stalled due to a LEF-BAR collision.
  // The logic is the same as what was described in the previous section.
  // There may be a way to handle this case directly in the first pass, but for the time being, this
  // will have to do.
  for (auto fwd_idx : fwd_ranks) {
    if (auto fwd_collision = fwd_collisions[fwd_idx]; is_lef_lef_primary_collision(fwd_collision))
      MODLE_UNLIKELY {
        fwd_collision -= lower_bound;
        const auto& rev_idx = fwd_collision;

        if (const auto& rev_collision = rev_collisions[rev_idx];
            is_lef_bar_collision(rev_collision)) {
          const auto& rev_unit = lefs[rev_idx].rev_unit;
          const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

          const auto& rev_move = rev_moves[rev_idx];
          auto& fwd_move = fwd_moves[fwd_idx];

          assert(rev_unit.pos() >= fwd_unit.pos() + rev_move);  // NOLINT
          fwd_move = (rev_unit.pos() - rev_move) - fwd_unit.pos();
        }
      }
  }
}

inline void Simulation::correct_moves_for_secondary_lef_lef_collisions(
    const absl::Span<const Lef> lefs, const size_t nbarriers, absl::Span<const size_t> rev_ranks,
    absl::Span<const size_t> fwd_ranks, const absl::Span<bp_t> rev_moves,
    const absl::Span<bp_t> fwd_moves, const absl::Span<const collision_t> rev_collisions,
    const absl::Span<const collision_t> fwd_collisions, const size_t nrev_units_at_5prime,
    const size_t nfwd_units_at_3prime) noexcept(utils::ndebug_defined()) {
  // Secondary LEF-LEF collisions are encoded with a number between nbarriers + lefs and nbarriers +
  // 2*nlefs.
  // Given an extrusion unit involved in a secondary LEF-LEF collision, the index of the extr. unit
  // that is causing the collision is encoded as nbarriers + nlefs + i.
  const auto lower_bound = nbarriers + lefs.size();
  const auto upper_bound = lower_bound + lefs.size();

  auto is_secondary_lef_lef_collision = [&](const auto i) constexpr noexcept {
    return i >= lower_bound && i < upper_bound;
  };

  // Skip over rev units that are already at the 5'-end. Furthermore, the left-most rev unit cannot
  // be involved in a secondary LEF-LEF collision, so we skip it.
  rev_ranks.remove_prefix(nrev_units_at_5prime + 1);
  for (const auto rev_idx2 : rev_ranks) {
    if (auto rev_idx1 = rev_collisions[rev_idx2]; is_secondary_lef_lef_collision(rev_idx1))
      MODLE_UNLIKELY {
        rev_idx1 -= lower_bound;  // Decode the index
        const auto& rev_unit1 = lefs[rev_idx1].rev_unit;
        const auto& rev_unit2 = lefs[rev_idx2].rev_unit;

        assert(rev_unit1.pos() >= rev_moves[rev_idx1]);  // NOLINT
        assert(rev_unit2.pos() >= rev_moves[rev_idx2]);  // NOLINT
        assert(rev_unit1.pos() - rev_moves[rev_idx1] >=
               rev_unit2.pos() - rev_moves[rev_idx2]);  // NOLINT

        // Update the move of the second unit such that after extruding once, this unit will be
        // located one bp downstream of the first rev unit
        const auto move = rev_unit2.pos() - (rev_unit1.pos() - rev_moves[rev_idx1]);
        rev_moves[rev_idx2] = move > 0UL ? move - 1 : 0UL;
      }
  }

  // Look above for detailed comments.
  fwd_ranks.remove_suffix(nfwd_units_at_3prime + 1);
  for (const auto fwd_idx1 : boost::adaptors::reverse(fwd_ranks)) {
    if (auto fwd_idx2 = fwd_collisions[fwd_idx1]; is_secondary_lef_lef_collision(fwd_idx2))
      MODLE_UNLIKELY {
        fwd_idx2 -= lower_bound;
        const auto& fwd_unit1 = lefs[fwd_idx1].fwd_unit;
        const auto& fwd_unit2 = lefs[fwd_idx2].fwd_unit;

        assert(fwd_unit1.pos() <= fwd_unit2.pos());  // NOLINT
        assert(fwd_unit1.pos() + fwd_moves[fwd_idx1] >=
               fwd_unit2.pos() + fwd_moves[fwd_idx2]);  // NOLINT

        const auto delta = (fwd_unit2.pos() + fwd_moves[fwd_idx2]) - fwd_unit1.pos();
        fwd_moves[fwd_idx1] = delta > 0UL ? delta - 1 : 0UL;
      }
  }
}

}  // namespace modle