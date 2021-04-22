#pragma once

#include <absl/types/span.h>  // for Span

#include <boost/range/adaptor/reversed.hpp>  // for reversed_range, reverse
#include <cassert>                           // for assert
#include <cstddef>                           // for size_t

#include "modle/common.hpp"             // for bp_t, MODLE_UNLIKELY
#include "modle/extrusion_factors.hpp"  // for ExtrusionUnit, Lef
#include "modle/utils.hpp"              // for ndebug_defined

namespace modle {

template <typename I>
void Simulation::adjust_moves_for_lef_bar_collisions(
    absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
    absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
    absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");

  const auto upper_bound = barriers.size();

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (rev_collisions[i] < upper_bound) MODLE_UNLIKELY {
        const auto& barrier_idx = rev_collisions[i];
        const auto& barrier = barriers[barrier_idx];
        assert(lefs[i].rev_unit.pos() > barrier.pos());  // NOLINT

        const auto delta = lefs[i].rev_unit.pos() - barrier.pos();
        rev_moves[i] = delta > 1UL ? delta - 1 : 0UL;
      }

    if (fwd_collisions[i] < upper_bound) MODLE_UNLIKELY {
        const auto& barrier_idx = fwd_collisions[i];
        const auto& barrier = barriers[barrier_idx];
        assert(lefs[i].fwd_unit.pos() < barrier.pos());  // NOLINT

        const auto delta = barrier.pos() - lefs[i].fwd_unit.pos();
        fwd_moves[i] = delta > 1UL ? delta - 1 : 0UL;
      }
  }
}

template <typename I>
inline void Simulation::adjust_moves_for_primary_lef_lef_collisions(
    absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
    absl::Span<const std::size_t> rev_ranks, absl::Span<const std::size_t> fwd_ranks,
    absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
    absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");

  (void)fwd_ranks;  // TODO Remove me
  const auto lower_bound = barriers.size();
  const auto upper_bound = lower_bound + lefs.size();

  auto is_lef_lef_primary_collision = [&](const auto i) constexpr {
    return i >= lower_bound && i < upper_bound;
  };

  auto is_lef_bar_collision = [&](const auto i) constexpr { return i < lower_bound; };

  for (auto rev_idx : rev_ranks) {
    if (auto rev_collision = rev_collisions[rev_idx]; is_lef_lef_primary_collision(rev_collision))
      MODLE_UNLIKELY {
        rev_collision -= lower_bound;
        const auto& fwd_idx = rev_collision;
        if (auto fwd_collision = fwd_collisions[fwd_idx];
            is_lef_lef_primary_collision(fwd_collision)) {
          const auto& rev_unit = lefs[rev_idx].rev_unit;
          const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

          auto& rev_move = rev_moves[rev_idx];
          auto& fwd_move = fwd_moves[fwd_idx];

          const auto [p1, p2] =
              compute_lef_lef_collision_pos(rev_unit, fwd_unit, rev_move, fwd_move);
          assert(rev_unit.pos() >= p1);  // NOLINT
          assert(fwd_unit.pos() <= p2);  // NOLINT

          rev_move = rev_unit.pos() - p1;
          fwd_move = p2 - fwd_unit.pos();
        } else if (is_lef_bar_collision(fwd_collision)) {
          const auto& rev_unit = lefs[rev_idx].rev_unit;
          const auto& fwd_unit = lefs[fwd_idx].fwd_unit;

          auto& rev_move = rev_moves[rev_idx];
          const auto& fwd_move = fwd_moves[fwd_idx];

          assert(rev_unit.pos() >= fwd_unit.pos() + fwd_move);  // NOLINT
          rev_move = rev_unit.pos() - (fwd_unit.pos() + fwd_move);
        }
      }
  }

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

template <typename I>
inline void Simulation::adjust_moves_for_secondary_lef_lef_collisions(
    absl::Span<const Lef> lefs, std::size_t nbarriers, absl::Span<const std::size_t> rev_ranks,
    absl::Span<const std::size_t> fwd_ranks, absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
    absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined()) {
  static_assert(std::is_integral_v<I>, "Collision buffers should be of integral type.");

  const auto lower_bound = nbarriers + lefs.size();
  const auto upper_bound = lower_bound + lefs.size();

  for (const auto rev_idx2 : rev_ranks) {  // TODO skip the first rank
    if (auto rev_idx1 = rev_collisions[rev_idx2]; rev_idx1 >= lower_bound && rev_idx1 < upper_bound)
      MODLE_UNLIKELY {
        rev_idx1 -= lower_bound;
        const auto& rev_unit1 = lefs[rev_idx1].rev_unit;
        const auto& rev_unit2 = lefs[rev_idx2].rev_unit;

        assert(rev_unit1.pos() >= rev_moves[rev_idx1]);  // NOLINT
        assert(rev_unit2.pos() >= rev_moves[rev_idx2]);  // NOLINT
        assert(rev_unit1.pos() - rev_moves[rev_idx1] >=
               rev_unit2.pos() - rev_moves[rev_idx2]);  // NOLINT

        const auto move = rev_unit2.pos() - (rev_unit1.pos() - rev_moves[rev_idx1]);
        rev_moves[rev_idx2] = move > 0UL ? move - 1 : 0UL;
      }
  }

  for (const auto fwd_idx1 : boost::adaptors::reverse(fwd_ranks)) {  // TODO skip the last rank
    if (auto fwd_idx2 = fwd_collisions[fwd_idx1]; fwd_idx2 >= lower_bound && fwd_idx2 < upper_bound)
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

template <typename I>
void Simulation::adjust_moves_based_on_collisions(
    absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
    absl::Span<const std::size_t> rev_ranks, absl::Span<const std::size_t> fwd_ranks,
    absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
    absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined()) {
  Simulation::adjust_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                  rev_collisions, fwd_collisions);
  Simulation::adjust_moves_for_primary_lef_lef_collisions(
      lefs, barriers, rev_ranks, fwd_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);
  Simulation::adjust_moves_for_secondary_lef_lef_collisions(lefs, barriers.size(), rev_ranks,
                                                            fwd_ranks, rev_moves, fwd_moves,
                                                            rev_collisions, fwd_collisions);
}
}  // namespace modle
