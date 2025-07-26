// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "modle/simulation.hpp"
// clang-format on

#include <cassert>
#include <span>

#include "modle/common/common.hpp"
#include "modle/common/utils.hpp"
#include "modle/extrusion_barriers.hpp"
#include "modle/extrusion_factors.hpp"

namespace modle {

void Simulation::correct_moves_for_lef_bar_collisions(
    const std::span<const Lef> lefs, const ExtrusionBarriers& barriers,
    const std::span<bp_t> rev_moves, const std::span<bp_t> fwd_moves,
    const std::span<const CollisionT> rev_collisions,
    const std::span<const CollisionT> fwd_collisions) noexcept(utils::ndebug_defined()) {
  for (std::size_t i = 0; i < lefs.size(); ++i) {
    if (rev_collisions[i].collision_occurred(CollisionT::LEF_BAR))
        [[unlikely]] {  // Process rev collisions
      const auto barrier_idx = rev_collisions[i].decode_index();
      const auto barrier_pos = barriers.pos(barrier_idx);
      assert(lefs[i].rev_unit.pos() > barrier_pos);

      // Compute the distance and update the rev_move such that after extruding once, the rev unit
      // of the current LEF will be located 1bp downstream of the extr. barrier.
      const auto distance = lefs[i].rev_unit.pos() - barrier_pos;
      assert(distance != 0);
      rev_moves[i] = distance - 1;
    }

    if (fwd_collisions[i].collision_occurred(CollisionT::LEF_BAR))
        [[unlikely]] {  // Process fwd collisions
      const auto barrier_idx = fwd_collisions[i].decode_index();
      const auto barrier_pos = barriers.pos(barrier_idx);
      assert(lefs[i].fwd_unit.pos() < barrier_pos);

      // Same as above. In this case the unit will be located 1bp upstream of the extr. barrier.
      const auto distance = barrier_pos - lefs[i].fwd_unit.pos();
      assert(distance != 0);
      fwd_moves[i] = distance - 1;
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Simulation::correct_moves_for_primary_lef_lef_collisions(
    const std::span<const Lef> lefs, const std::span<const std::size_t> rev_ranks,
    const std::span<const std::size_t> fwd_ranks, const std::span<bp_t> rev_moves,
    const std::span<bp_t> fwd_moves, const std::span<const CollisionT> rev_collisions,
    const std::span<const CollisionT> fwd_collisions) noexcept(utils::ndebug_defined()) {
  // Primary LEF-LEF collisions are encoded with a number between nbarriers and nbarriers + nlefs.
  // Given a pair of extr. units that are moving in opposite directions, the index i corresponding
  // to the extr. unit that is causing the collision is encoded as nbarriers + i.
  for (auto rev_idx : rev_ranks) {  // Loop over rev units in 5'-3' order
    if (rev_collisions[rev_idx].collision_occurred(CollisionT::LEF_LEF_PRIMARY)) [[unlikely]] {
      const auto fwd_idx = rev_collisions[rev_idx].decode_index();

      if (fwd_collisions[fwd_idx].collision_occurred(CollisionT::LEF_LEF_PRIMARY)) {
        // This branch handles the typical case, where a pair of extr. units moving in opposite
        // direction bumped into each other, causing a primary LEF-LEF collision.
        const auto& rev_unit = lefs[rev_idx].rev_unit;
        const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

        auto& rev_move = rev_moves[rev_idx];
        auto& fwd_move = fwd_moves[fwd_idx];

        const auto [p1, p2] = compute_lef_lef_collision_pos(rev_unit, fwd_unit, rev_move, fwd_move);
        assert(rev_unit.pos() >= p1);
        assert(fwd_unit.pos() <= p2);

        // Update the moves of the involved units, so that after the next call to extrude these
        // two units will be located nex to each others at the collision site.
        rev_move = rev_unit.pos() - p1;
        fwd_move = p2 - fwd_unit.pos();

      } else if (fwd_collisions[fwd_idx].collision_occurred(CollisionT::LEF_BAR)) {
        // This branch handles the special case where the fwd unit involved in the collision is
        // blocked by an extrusion barrier, and thus cannot be moved.
        // In this case we update the move such that after extrusion, the rev unit will be located
        // next to the fwd unit.
        const auto& rev_unit = lefs[rev_idx].rev_unit;
        const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

        auto& rev_move = rev_moves[rev_idx];
        const auto& fwd_move = fwd_moves[fwd_idx];

        assert(rev_unit.pos() >= fwd_unit.pos() + fwd_move);
        rev_move = rev_unit.pos() - (fwd_unit.pos() + fwd_move) - 1;
      }
    }
  }

  // This second loop is required in order to properly handle the scenario where a fwd unit collides
  // with a rev unit that is stalled due to a LEF-BAR collision.
  // The logic is the same as what was described in the previous section.
  // There may be a way to handle this case directly in the first pass, but for the time being, this
  // will have to do.
  for (auto fwd_idx : fwd_ranks) {
    if (fwd_collisions[fwd_idx].collision_occurred(CollisionT::LEF_LEF_PRIMARY)) [[unlikely]] {
      const auto rev_idx = fwd_collisions[fwd_idx].decode_index();

      if (rev_collisions[rev_idx].collision_occurred(CollisionT::LEF_BAR)) {
        const auto& rev_unit = lefs[rev_idx].rev_unit;
        const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

        const auto& rev_move = rev_moves[rev_idx];
        auto& fwd_move = fwd_moves[fwd_idx];

        assert(rev_unit.pos() >= fwd_unit.pos() + rev_move);
        fwd_move = (rev_unit.pos() - rev_move) - fwd_unit.pos() - 1;
      }
    }
  }
}
}  // namespace modle
