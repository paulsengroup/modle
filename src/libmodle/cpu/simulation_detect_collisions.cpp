// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "modle/simulation.hpp"
// clang-format on

#include <absl/types/span.h>

#include <algorithm>
#include <cassert>
#include <limits>
#include <utility>

#include "modle/collision_encoding.hpp"
#include "modle/common/common.hpp"
#include "modle/common/dna.hpp"
#include "modle/extrusion_barriers.hpp"
#include "modle/extrusion_factors.hpp"
#include "modle/genome.hpp"

namespace modle {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
std::pair<std::size_t, std::size_t> Simulation::detect_units_at_interval_boundaries(
    const GenomicInterval& interval, absl::Span<const Lef> lefs,
    absl::Span<const std::size_t> rev_lef_ranks, absl::Span<const std::size_t> fwd_lef_ranks,
    absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
    absl::Span<CollisionT> rev_collisions, absl::Span<CollisionT> fwd_collisions) {
  {
    assert(lefs.size() == fwd_lef_ranks.size());
    assert(lefs.size() == rev_lef_ranks.size());
    assert(lefs.size() == fwd_moves.size());
    assert(lefs.size() == rev_moves.size());
    assert(lefs.size() == fwd_collisions.size());
    assert(lefs.size() == rev_collisions.size());
    assert(std::is_sorted(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                          }));
    assert(std::is_sorted(rev_lef_ranks.begin(), rev_lef_ranks.end(),
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                          }));
    assert(std::all_of(rev_collisions.begin(), rev_collisions.end(),
                       [&](const CollisionT c) { return c() == 0; }));
    assert(std::all_of(fwd_collisions.begin(), fwd_collisions.end(),
                       [&](const CollisionT c) { return c() == 0; }));
  }
  // Detect if the first rev unit or last fwd unit are about to fall off interval. boundaries
  // Also detect extr. units that are already at interval boundaries

  std::size_t num_rev_units_at_5prime = 0;
  std::size_t num_fwd_units_at_3prime = 0;

  const auto& first_active_fwd_unit = lefs[fwd_lef_ranks[0]].fwd_unit;
  const auto& last_active_rev_unit =
      lefs[*std::find_if(rev_lef_ranks.rbegin(), rev_lef_ranks.rend(), [&](const auto i) {
        return lefs[i].is_bound();
      })].rev_unit;

  // Detect and count the number of rev units located at the 5'-end
  for (std::size_t i = 0; i < lefs.size(); ++i) {
    const auto& rev_idx = rev_lef_ranks[i];
    const auto& rev_unit = lefs[rev_idx].rev_unit;
    const auto& rev_move = rev_moves[rev_idx];

    assert(lefs[rev_idx].is_bound());
    assert(interval.start() + rev_move <= rev_unit.pos());

    const std::size_t _5prime_idx = 5;
    if (rev_unit.pos() == interval.start()) {  // Unit already at the 5'-end
      assert(rev_moves[rev_idx] == 0);
      ++num_rev_units_at_5prime;
      rev_collisions[rev_idx].set(_5prime_idx, CollisionT::COLLISION | CollisionT::CHROM_BOUNDARY);

    } else if (rev_unit.pos() > first_active_fwd_unit.pos()) {
      // As the rev_unit is located downstream of a fwd unit, a LEF-LEF collision might take place
      // before the rev unit reaches the 5'-end
      break;

    } else if (rev_unit.pos() - rev_move == interval.start()) {
      // Unit will reach the 5'-end by the end of the current epoch
      rev_collisions[rev_idx].set(_5prime_idx, CollisionT::COLLISION | CollisionT::CHROM_BOUNDARY);
      ++num_rev_units_at_5prime;
      break;
    }
  }

  // Detect and count the number of fwd units located at the 3'-end. See above for detailed comments
  for (auto i = lefs.size() - 1; i > 0; --i) {
    const auto& fwd_idx = fwd_lef_ranks[i];
    const auto& fwd_unit = lefs[fwd_idx].fwd_unit;
    const auto& fwd_move = fwd_moves[fwd_idx];

    if (!lefs[fwd_idx].is_bound()) {
      ++num_fwd_units_at_3prime;  // Inactive units are technically not at the 3'-end, but we count
      continue;                   // them anyway so that we can shrink the spans on LEF-related
                                  // buffers to avoid doing some work in later steps
    }

    const std::size_t _3prime_idx = 3;
    assert(fwd_unit.pos() + fwd_move < interval.end());
    if (fwd_unit.pos() == interval.end() - 1) {
      assert(fwd_moves[fwd_idx] == 0);
      ++num_fwd_units_at_3prime;
      fwd_collisions[fwd_idx].set(_3prime_idx, CollisionT::COLLISION | CollisionT::CHROM_BOUNDARY);

    } else if (fwd_unit.pos() < last_active_rev_unit.pos()) {
      break;

    } else if (fwd_unit.pos() + fwd_move == interval.end() - 1) {
      fwd_collisions[fwd_idx].set(_3prime_idx, CollisionT::COLLISION | CollisionT::CHROM_BOUNDARY);
      ++num_fwd_units_at_3prime;
      break;
    }
  }

  return std::make_pair(num_rev_units_at_5prime, num_fwd_units_at_3prime);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Simulation::detect_lef_bar_collisions(
    const absl::Span<const Lef> lefs, const absl::Span<const std::size_t> rev_lef_ranks,
    const absl::Span<const std::size_t> fwd_lef_ranks, const absl::Span<const bp_t> rev_moves,
    const absl::Span<const bp_t> fwd_moves, const ExtrusionBarriers& barriers,
    const absl::Span<CollisionT> rev_collisions, const absl::Span<CollisionT> fwd_collisions,
    random::PRNG_t& rand_eng, std::size_t num_rev_units_at_5prime,
    std::size_t num_fwd_units_at_3prime) const noexcept(utils::ndebug_defined()) {
  {
    assert(lefs.size() == fwd_lef_ranks.size());
    assert(lefs.size() == rev_lef_ranks.size());
    assert(lefs.size() == fwd_moves.size());
    assert(lefs.size() == rev_moves.size());
    assert(lefs.size() == fwd_collisions.size());
    assert(lefs.size() == rev_collisions.size());
    assert(lefs.size() >= num_rev_units_at_5prime);
    assert(lefs.size() >= num_fwd_units_at_3prime);
    assert(std::is_sorted(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                          }));
    assert(std::is_sorted(rev_lef_ranks.begin(), rev_lef_ranks.end(),
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                          }));
  }
  // Loop over LEFs, using a procedure similar to merge in mergesort.
  // The idea here is that if we have a way to visit extr. units in 5'-3' order, then detecting
  // LEF-BAR collisions boils down to the following:
  //  - Loop over all extrusion barriers:
  //    - Find the first rev. unit located downstream of the current extr. barrier
  //    - Find the last fwd. unit located upstream of the current extr. barrier
  //    - If the distance between the extr. barrier and rev/fwd extr. units is less than the
  //      distance that will be covered by the rev/fwd extr. unit in the current iteration, then
  //      there may be a collision.
  //      To determine whether or not there will be a collision, run a bernoulli trial with a
  //      probability of success equal to the probability of block of the current barrier
  //
  // When a collision is detected the appropriate entry in rev/fwd_collision buffer is set to the
  // index of the extrusion barrier that caused the collision.
  // Furthermore the appropriate entry in the move vector is updated such that after calling
  // Simulation::extrude, the extr. unit that is being blocked will be located 1pb up/downstream
  // of the extr. barrier that caused the collision.

  // Init indices and position with the rev and fwd extr. units that are the closest to the 5'-end
  auto j = std::min(num_rev_units_at_5prime, num_rev_units_at_5prime - 1);

  auto unit_idx = rev_lef_ranks[j];
  auto unit_pos = lefs[unit_idx].rev_unit.pos();

  // Loop over extr. barriers and find the first, possibly colliding extr. unit
  for (std::size_t i = 0; i < barriers.size(); ++i) {
    if (barriers.is_not_active(i)) {  // Inactive barriers are trnsparent to extr. units
      continue;
    }

    assert(j < lefs.size());
    // Probability of block is set based on the extr. barrier blocking direction
    const auto& pblock = barriers.direction(i) == dna::REV ? c().lef_bar_major_collision_pblock
                                                           : c().lef_bar_minor_collision_pblock;

    // Look for the first rev extr. unit that comes after the current barrier
    while (unit_pos <= barriers.pos(i)) {
      if (++j == lefs.size()) [[unlikely]] {  // All rev units have been processed
        goto process_fwd_unit;                // Move to the next section
      }

      // Update idx and position with those corresponding to the rev extr. unit that comes next
      unit_idx = rev_lef_ranks[j];  // in 5'-3' order
      unit_pos = lefs[unit_idx].rev_unit.pos();
    }

    if (lefs[unit_idx].is_bound()) [[unlikely]] {
      // We have a LEF-BAR collision event if the distance between the rev. unit and the extr.
      // barrier is less or equal than the distance that the rev extr. unit is set to move in
      // the current iteration. If pblock != 1, then we also require a successful bernoulli
      // trial before calling a collision
      assert(unit_pos >= barriers.pos(i));
      const auto delta = unit_pos - barriers.pos(i);
      if (delta > 0 && delta <= rev_moves[unit_idx] &&
          Simulation::run_lef_bar_collision_trial(pblock, rand_eng)) {
        // Collision detected
        rev_collisions[unit_idx].set(i, CollisionT::COLLISION | CollisionT::LEF_BAR);
      }
    }
  }

// Look in the previous section for detailed comments
process_fwd_unit:
  // Init indices and position with the rev and fwd extr. units that are the closest to the 5'-end
  j = lefs.size() - std::min(num_fwd_units_at_3prime, num_fwd_units_at_3prime - 1);

  assert(j != 0);
  unit_idx = fwd_lef_ranks[--j];
  unit_pos = lefs[unit_idx].fwd_unit.pos();

  // Loop over extr. barriers and find the first, possibly colliding extr. unit
  assert(!barriers.empty());
  const auto sentinel_idx = (std::numeric_limits<std::size_t>::max)();
  for (auto i = barriers.size() - 1; i != sentinel_idx; --i) {
    if (barriers.is_not_active(i)) {  // Inactive barriers are trnsparent to extr. units
      continue;
    }

    // Probability of block is set based on the extr. barrier blocking direction
    const auto& pblock = barriers.direction(i) == dna::FWD ? c().lef_bar_major_collision_pblock
                                                           : c().lef_bar_minor_collision_pblock;
    // Look for the next fwd unit that comes strictly before the current extr. barrier
    while (unit_pos >= barriers.pos(i)) {
      if (--j == sentinel_idx) [[unlikely]] {
        return;
      }

      unit_idx = fwd_lef_ranks[j];
      unit_pos = lefs[unit_idx].fwd_unit.pos();
    }

    if (lefs[unit_idx].is_bound()) [[unlikely]] {
      const auto delta = barriers.pos(i) - unit_pos;
      if (delta > 0 && delta <= fwd_moves[unit_idx] &&
          Simulation::run_lef_bar_collision_trial(pblock, rand_eng)) {
        fwd_collisions[unit_idx].set(i, CollisionT::COLLISION | CollisionT::LEF_BAR);
      }
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Simulation::detect_primary_lef_lef_collisions(
    const absl::Span<const Lef> lefs, const ExtrusionBarriers& barriers,
    const absl::Span<const std::size_t> rev_lef_ranks,
    const absl::Span<const std::size_t> fwd_lef_ranks, const absl::Span<const bp_t> rev_moves,
    const absl::Span<const bp_t> fwd_moves, const absl::Span<CollisionT> rev_collisions,
    const absl::Span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng,
    std::size_t num_rev_units_at_5prime, std::size_t num_fwd_units_at_3prime) const
    noexcept(utils::ndebug_defined()) {
  {
    assert(!lefs.empty());
    assert(lefs.size() == fwd_lef_ranks.size());
    assert(lefs.size() == rev_lef_ranks.size());
    assert(lefs.size() == fwd_moves.size());
    assert(lefs.size() == rev_moves.size());
    assert(lefs.size() == fwd_collisions.size());
    assert(lefs.size() == rev_collisions.size());
    assert(lefs.size() >= num_rev_units_at_5prime);
    assert(lefs.size() >= num_fwd_units_at_3prime);
    assert(std::is_sorted(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                          }));
    assert(std::is_sorted(rev_lef_ranks.begin(), rev_lef_ranks.end(),
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                          }));
  }
  // Loop over LEFs, using a procedure similar to merge in mergesort
  // The idea here is that if we have a way to visit extr. units in 5'-3' order, then detecting
  // LEF-LEF collisions boils down to:
  //  - Starting from the first fwd extr. unit in 5'-3' order:
  //    - Look for the first rev extr. unit that is located downstream of the fwd unit that is being
  //    processed
  //    - If the distance between the pair of extr. units identified in the previous step is below a
  //      certain threshold, then we may have a LEF-LEF collision.
  //      If the probability of unit bypass is 0, then we always predict a LEF-LEF collision. If
  //      this probability is larger than 0, then we predict a LEF-LEF collision based on the
  //      outcome of a Bernoulli trial with probability of success equal to 1 - the prob. of bypass
  //    - If a LEF-LEF collision caused by two extr. unit moving in opposite direction is detected,
  //      encode the index of the extr. unit that caused the collision in the appropriate collision
  //      mask. This kind of collisions are encoded as offset + i, where offset = nbarriers and i =
  //      the index of the unit that is colliding.

  if (num_rev_units_at_5prime == lefs.size() || num_fwd_units_at_3prime == lefs.size())
      [[unlikely]] {
    return;
  }

  //    Initialize indexes so that we skip over rev units at the 5' and fwd units at the 3' (if any)
  std::size_t i1 = 0;
  auto j1 = num_rev_units_at_5prime;
  const auto i2 = lefs.size() - std::min(num_fwd_units_at_3prime, num_fwd_units_at_3prime - 1);
  const auto j2 = lefs.size();

  while (true) {
    auto rev_idx = rev_lef_ranks[j1];             // index of the jth rev unit in 5'-3' order
    auto rev_pos = lefs[rev_idx].rev_unit.pos();  // pos of the jth unit
    auto fwd_idx = fwd_lef_ranks[i1];             // index of the ith fwd unit in 5'-3' order
    auto fwd_pos = lefs[fwd_idx].fwd_unit.pos();  // pos of the ith unit

    // Find the first rev unit that comes right after the ith fwd unit
    while (rev_pos <= fwd_pos) {
      if (++j1 == j2) [[unlikely]] {
        return;  // all rev units have been processed
      }
      rev_idx = rev_lef_ranks[j1];
      rev_pos = lefs[rev_idx].rev_unit.pos();
    }

    // Find the last fwd unit that comes right before the jth rev unit
    // This is necessary in order to handle ties in fwd units ranking, as well as the case where
    // there are many fwd units between the pos of the jth-1 and jth rev extr. units
    while (fwd_pos < rev_pos) {
      if (++i1 == i2) [[unlikely]] {
        return;  // all fwd units have been processed
      }
      fwd_idx = fwd_lef_ranks[i1];             // index of the ith fwd unit in 5'-3' order
      fwd_pos = lefs[fwd_idx].fwd_unit.pos();  // pos of the ith fwd unit
    }

    // The previous while loop finds the first fwd unit that comes at or after the rev unit that
    // is being processed, but we are actually interested in the fwd unit before this one.
    // This is why we take the unit at index i1-1
    fwd_idx = fwd_lef_ranks[std::min(i1, i1 - 1)];
    fwd_pos = lefs[fwd_idx].fwd_unit.pos();

    // We have a LEF-LEF collision event if the distance between the two extr. units is less than
    // the sum of their moves.
    // If probability_of_extrusion_unit_bypass != 0, then we also require a successful bernoulli
    // trial before calling a collision
    if (const auto delta = rev_pos - fwd_pos; delta > 0 &&
                                              delta < rev_moves[rev_idx] + fwd_moves[fwd_idx] &&
                                              run_lef_lef_collision_trial(rand_eng)) {
      // Declare few aliases to reduce code verbosity later on
      const auto& rev_move = rev_moves[rev_idx];
      const auto& fwd_move = fwd_moves[fwd_idx];
      auto& rev_collision = rev_collisions[rev_idx];
      auto& fwd_collision = fwd_collisions[fwd_idx];

      // Note that the rev collision pos will always be 1bp upstream of the fwd collision pos
      auto [collision_pos_rev, collision_pos_fwd] = compute_lef_lef_collision_pos(
          lefs[rev_idx].rev_unit, lefs[fwd_idx].fwd_unit, rev_move, fwd_move);

      if (!rev_collision.collision_occurred() && !fwd_collision.collision_occurred()) {
        // In the simplest case both units are free to move.
        rev_collision.set(fwd_idx, CollisionT::COLLISION | CollisionT::LEF_LEF_PRIMARY);
        fwd_collision.set(rev_idx, CollisionT::COLLISION | CollisionT::LEF_LEF_PRIMARY);

      } else if (rev_collision.collision_occurred() && !fwd_collision.collision_occurred()) {
        // In this case only the fwd unit is free to move.
        // This case is a bit more complicated than the first one, because we have to handle the
        // edge case where we have mistakenly predicted a LEF-BAR collision. This usually happens
        // when two extr. units are moving in opposite directions, and the unit extruding in rev
        // direction is located downstream of the unit extruding in fwd direction. If located
        // between the two units there's an extr. barrier, and one of the two units manages to
        // bypass the extr. barrier, then there's a chance that the LEF-LEF collision will take
        // place before the LEF-BAR collision (whether this happens or not depends on the distance
        // between the extr. units and the extr. barriers, as well as their respective extr. speed)

        assert(rev_collision.collision_occurred(CollisionT::LEF_BAR));
        assert(collision_pos_rev != 0 && collision_pos_fwd != 0);
        const auto barrier_pos = barriers.pos(rev_collision.decode_index());
        if (collision_pos_fwd > barrier_pos) [[unlikely]] {
          // Detected the mis-prediction mentioned above: make the LEF-BAR collision a LEF-LEF
          // collision
          rev_collision.set(fwd_idx, CollisionT::COLLISION | CollisionT::LEF_LEF_PRIMARY);
          fwd_collision.set(rev_idx, CollisionT::COLLISION | CollisionT::LEF_LEF_PRIMARY);
        } else {
          assert(rev_collision.collision_occurred(CollisionT::LEF_BAR));

          // fwd extr unit is being blocked by a rev unit that is itself being blocked by an extr.
          // barrier
          fwd_collision.set(rev_idx, CollisionT::COLLISION | CollisionT::LEF_LEF_PRIMARY);
        }
        // This branch follows the same logic as the previous one. In this case the rev unit is free
        // to move, while the fwd unit has been predicted to be stalled by an extr. barrier
      } else if (!rev_collision.collision_occurred() && fwd_collision.collision_occurred()) {
        assert(fwd_collision.collision_occurred(CollisionT::LEF_BAR));
        assert(collision_pos_rev != 0 && collision_pos_fwd != 0);
        const auto barrier_pos = barriers.pos(fwd_collision.decode_index());
        rev_collision.set(fwd_idx, CollisionT::COLLISION | CollisionT::LEF_LEF_PRIMARY);
        if (collision_pos_rev < barrier_pos) [[unlikely]] {
          fwd_collision.set(rev_idx, CollisionT::COLLISION | CollisionT::LEF_LEF_PRIMARY);
        }
      }
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Simulation::process_secondary_lef_lef_collisions(
    [[maybe_unused]] const GenomicInterval& interval, const absl::Span<const Lef> lefs,
    const absl::Span<const std::size_t> rev_lef_ranks,
    const absl::Span<const std::size_t> fwd_lef_ranks, const absl::Span<bp_t> rev_moves,
    const absl::Span<bp_t> fwd_moves, const absl::Span<CollisionT> rev_collisions,
    const absl::Span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng,
    std::size_t num_rev_units_at_5prime, std::size_t num_fwd_units_at_3prime) const
    noexcept(utils::ndebug_defined()) {
  {
    assert(lefs.size() == fwd_lef_ranks.size());
    assert(lefs.size() == rev_lef_ranks.size());
    assert(lefs.size() == fwd_moves.size());
    assert(lefs.size() == rev_moves.size());
    assert(lefs.size() == fwd_collisions.size());
    assert(lefs.size() == rev_collisions.size());
    assert(lefs.size() >= num_rev_units_at_5prime);
    assert(lefs.size() >= num_fwd_units_at_3prime);
  }

  // Loop over pairs of consecutive fwd units.
  // Throughout the comments for this function we will use the following notation:
  //  - U1 = unit at index i - 1 (i.e. the unit closer to the 5'-end);
  //  - U2 = unit at index i (i.e. the unit closer to the 3'-end);
  //
  // 0. Keep advancing the index i until we find the first U2 unit that is stalled. Then:
  //     1. Check if calling Simulation::extrude would cause U1 to end up at the same position or
  //        downstream of U2.
  //        If this is the case, then we may have a LEF-LEF collision.
  //          2. If the probability of bypass is 0, then we certainly have a collision, otherwise
  //             whether a collision will occur or not is decided by a bernoulli trial with
  //             probability of success equal to 1 - prob. of bypass.
  //          3. If the previous step determines that a LEF-LEF collision is about to take place,
  //             update the appropriate collision mask with the encoded index of the extr. unit that
  //             is causing the collision. Index i is encoded as nbarriers + nlefs + i.
  //          4. Keep backtracking repeating steps 1-4 until step 1 fails
  // 5. Continue iterating from where we left at step 0

  // Loop over rev units in 5'-3' direction, and skip over units that are located at the 5'-end
  for (auto i = std::max(std::size_t(1), num_rev_units_at_5prime); i < lefs.size(); ++i) {
    // index of the ith-1 rev unit in 5'-3' order (i.e. U1)
    const auto& rev_idx1 = rev_lef_ranks[i - 1];
    const auto& rev_pos1 = lefs[rev_idx1].rev_unit.pos();
    if (!rev_collisions[rev_idx1].collision_occurred()) [[likely]] {
      // If U1 is not stalled, then it is not possible to have a secondary LEF-LEF collision
      continue;
    }

    // index of the ith rev unit in 5'-3' order (i.e. U2)
    const auto& rev_idx2 = rev_lef_ranks[i];
    const auto& rev_pos2 = lefs[rev_idx2].rev_unit.pos();
    if (rev_collisions[rev_idx2].collision_occurred()) {
      // LEF-BAR and primary LEF-LEF collisions have precedence over secondary LEF-LEF collisions
      continue;
    }

    auto& move1 = rev_moves[rev_idx1];
    auto& move2 = rev_moves[rev_idx2];

    assert(rev_pos2 >= rev_pos1);
    assert(rev_pos1 >= move1);
    assert(rev_pos2 >= move2);
    assert(rev_pos1 - move1 >= interval.start());
    assert(rev_pos2 - move2 >= interval.start());

    if (rev_pos2 - move2 <= rev_pos1 - move1) {
      if (run_lef_lef_collision_trial(rand_eng)) {
        rev_collisions[rev_idx2].set(rev_idx1,
                                     CollisionT::COLLISION | CollisionT::LEF_LEF_SECONDARY);
        const auto move = rev_pos2 - (rev_pos1 - move1);
        move2 = std::min(move, move - 1);
      } else {
        assert(rev_collisions[rev_idx1].collision_occurred());
        rev_collisions[rev_idx2].set(rev_idx1, CollisionT::LEF_LEF_SECONDARY);
      }
    }
  }

  // Loop over fwd units in 3'-5' direction, and skip over units that are located at the 3'-end.
  // Look in the section above for detailed comments
  auto i = lefs.size() - std::min(num_fwd_units_at_3prime, num_fwd_units_at_3prime - 1) - 1;
  assert(i < lefs.size());
  for (; i > 0; --i) {
    // index of the ith fwd unit in 3'-5' order (i.e. U2)
    const auto& fwd_idx2 = fwd_lef_ranks[i];
    const auto& fwd_pos2 = lefs[fwd_idx2].fwd_unit.pos();
    if (!fwd_collisions[fwd_idx2].collision_occurred()) [[likely]] {
      // If U2 is not stalled, then it is not possible to have a secondary LEF-LEF collision
      continue;
    }

    const auto& fwd_idx1 = fwd_lef_ranks[i - 1];
    const auto& fwd_pos1 = lefs[fwd_idx1].fwd_unit.pos();
    if (fwd_collisions[fwd_idx1].collision_occurred()) {
      continue;
    }

    auto& move1 = fwd_moves[fwd_idx1];
    auto& move2 = fwd_moves[fwd_idx2];

    assert(fwd_pos2 >= fwd_pos1);
    assert(fwd_pos1 + move1 < interval.end());
    assert(fwd_pos2 + move2 < interval.end());

    if (fwd_pos1 + move1 >= fwd_pos2 + move2) {
      if (run_lef_lef_collision_trial(rand_eng)) {
        fwd_collisions[fwd_idx1].set(fwd_idx2,
                                     CollisionT::COLLISION | CollisionT::LEF_LEF_SECONDARY);
        const auto move = (fwd_pos2 + move2) - fwd_pos1;
        move1 = std::min(move, move - 1);
      } else {
        assert(fwd_collisions[fwd_idx2].collision_occurred());
        fwd_collisions[fwd_idx1].set(fwd_idx2, CollisionT::LEF_LEF_SECONDARY);
      }
    }
  }
}

void Simulation::fix_secondary_lef_lef_collisions(
    [[maybe_unused]] const GenomicInterval& interval, const absl::Span<Lef> lefs,
    const absl::Span<std::size_t> rev_lef_ranks, const absl::Span<std::size_t> fwd_lef_ranks,
    const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
    const absl::Span<CollisionT> rev_collisions, const absl::Span<CollisionT> fwd_collisions,
    const std::size_t num_rev_units_at_5prime,
    const std::size_t num_fwd_units_at_3prime) noexcept(utils::ndebug_defined()) {
  // Let us consider the following scenario:
  // 1>   2>X
  // Where 1> and 2> are two consecutive EUs moving in fwd direction and X is an obstacle (e.g. an
  // extrusion barrier) which is currently blocking EU2.
  // Let us also assume that the probability of LEF bypass is != 0 and that if unobstructed, EU1 is
  // set to move past EU2 and X.
  //
  // In this scenario EU2 is blocked by obstacle X, and so its move is set to 0.
  // Let us also assume that the Bernoulli trial to determine whether EU1 is involved in a
  // LEF-LEF secondary collision caused by EU2 fails (i.e. collision is avoided).
  //
  // If EU1's move is left unmodified, EU1 will bypass EU2, X and anything else long its path.
  //
  // This function is meant to detect and address this scenario by correcting EU1's move so that
  // after move EU1 is located downstream of EU2.
  //
  // The case shown above (where EU2 is 1bp upstream of X) is a bit more complicated for the
  // following reason:
  // MoDLE currently does not allow extrusion units to move backwards (i.e. against the direction of
  // extrusion). This means that it is not possible to modify EU1 and EU2 moves in such a way that
  // after extrusion EU1 is located downstream of EU2 but upstream of X.
  // To properly deal with this special case we swap the positions of EU1 and EU2 and we update the
  // moves so that after extrusion EU2 will be located 1bp upstream of EU1 which will itself ne
  // located 1bp upstream of X

  const auto num_active_fwd_units =
      lefs.size() - std::min(num_fwd_units_at_3prime, num_fwd_units_at_3prime - 1);

  // Loop over rev units in 5'-3' direction, and skip over units that are located at the 5'-end
  for (auto i = std::max(std::size_t(1), num_rev_units_at_5prime); i < lefs.size(); ++i) {
    // The case addressed by If EU2 (the unit corresponding to rev_idx2) avoided a secondary LEF-LEF
    // collision
    const auto& rev_idx2 = rev_lef_ranks[i];
    const auto& rev_collision2 = rev_collisions[rev_idx2];
    if (rev_collision2.collision_avoided(CollisionT::LEF_LEF_SECONDARY)) [[unlikely]] {
      const auto& rev_idx1 = rev_lef_ranks[i - 1];
      [[maybe_unused]] const auto& rev_collision1 = rev_collisions[rev_idx1];
      assert(rev_idx1 == rev_collision2.decode_index());  // idx of the EU blocking EU1
      assert(rev_collision1.collision_occurred());

      auto& lef1 = lefs[rev_idx1];
      auto& lef2 = lefs[rev_idx2];
      // Position after move of EU1
      const auto pos1 = lef1.rev_unit.pos() - rev_moves[rev_idx1];
      // If possible, adjust EU2's move so that after extrusion EU2 will be located 1bp upstream of
      // EU1
      if (lef2.rev_unit.pos() > pos1 + 1) [[likely]] {
        rev_moves[rev_idx2] = lef2.rev_unit.pos() - (pos1 + 1);
      } else {
        rev_moves[rev_idx2] = 0;
      }

      // Mark EU2 as involved in a secondary LEF-LEF collision
      rev_collisions[rev_idx2].set(rev_idx1, CollisionT::COLLISION | CollisionT::LEF_LEF_SECONDARY);

      // Swap EUs position while ensuring that the rev unit does not cross its corresponding fwd
      // unit
      const auto p1 = lef1.rev_unit._pos;
      const auto p2 = lef2.rev_unit._pos;
      lef1.rev_unit._pos = std::min(lef1.fwd_unit._pos, p2);
      lef2.rev_unit._pos = std::min(lef2.fwd_unit._pos, p1);

      // Swap collisions, moves and ranks
      std::swap(rev_collisions[rev_idx1], rev_collisions[rev_idx2]);
      std::swap(rev_moves[rev_idx1], rev_moves[rev_idx2]);
      std::swap(rev_lef_ranks[i - 1], rev_lef_ranks[i]);

      // Sometime when swapping moves it is possible that the resulting move will cause an extrusion
      // unit to cross chromosomal boundaries. Right now I can't think of a way around this other
      // than clamping the moves once again
      rev_moves[rev_lef_ranks[i - 1]] =
          std::min(lefs[rev_lef_ranks[i - 1]].rev_unit.pos() - interval.start(),
                   rev_moves[rev_lef_ranks[i - 1]]);
      rev_moves[rev_lef_ranks[i]] = std::min(
          lefs[rev_lef_ranks[i]].rev_unit.pos() - interval.start(), rev_moves[rev_lef_ranks[i]]);
      assert(lefs[rev_lef_ranks[i - 1]].rev_unit.pos() >= rev_moves[rev_lef_ranks[i - 1]]);
      assert(lefs[rev_lef_ranks[i]].rev_unit.pos() >= rev_moves[rev_lef_ranks[i]]);
    }
  }

  // Look above for detailed comments
  for (std::size_t i = 0; i < num_active_fwd_units - 1; ++i) {
    const auto& fwd_idx1 = fwd_lef_ranks[i];
    const auto& fwd_collision1 = fwd_collisions[fwd_idx1];
    if (fwd_collision1.collision_avoided(CollisionT::LEF_LEF_SECONDARY)) {
      const auto& fwd_idx2 = fwd_lef_ranks[i + 1];  // idx of the EU blocking EU1
      [[maybe_unused]] const auto& fwd_collision2 = fwd_collisions[fwd_idx2];
      assert(fwd_idx2 == fwd_collision1.decode_index());
      assert(fwd_collision2.collision_occurred());

      auto& lef1 = lefs[fwd_idx1];
      auto& lef2 = lefs[fwd_idx2];

      const auto pos2 = lef2.fwd_unit.pos() + fwd_moves[fwd_idx2];
      if (pos2 > lef1.fwd_unit.pos() + 1) {
        fwd_moves[fwd_idx1] = pos2 - (lef1.fwd_unit.pos() + 1);
      } else {
        fwd_moves[fwd_idx1] = 0;
      }
      fwd_collisions[fwd_idx1].set(fwd_idx2, CollisionT::COLLISION | CollisionT::LEF_LEF_SECONDARY);

      const auto p1 = lef1.fwd_unit._pos;
      const auto p2 = lef2.fwd_unit._pos;
      lef1.fwd_unit._pos = std::max(lef1.rev_unit._pos, p2);
      lef2.fwd_unit._pos = std::max(lef2.rev_unit._pos, p1);
      std::swap(fwd_collisions[fwd_idx1], fwd_collisions[fwd_idx2]);
      std::swap(fwd_moves[fwd_idx1], fwd_moves[fwd_idx2]);
      std::swap(fwd_lef_ranks[i], fwd_lef_ranks[i + 1]);

      fwd_moves[fwd_lef_ranks[i]] = std::min(
          interval.end() - 1 - lefs[fwd_lef_ranks[i]].fwd_unit.pos(), fwd_moves[fwd_lef_ranks[i]]);
      fwd_moves[fwd_lef_ranks[i + 1]] =
          std::min(interval.end() - 1 - lefs[fwd_lef_ranks[i + 1]].fwd_unit.pos(),
                   fwd_moves[fwd_lef_ranks[i + 1]]);

      assert(lefs[fwd_lef_ranks[i]].fwd_unit.pos() + fwd_moves[fwd_lef_ranks[i]] < interval.end());
      assert(lefs[fwd_lef_ranks[i + 1]].fwd_unit.pos() + fwd_moves[fwd_lef_ranks[i + 1]] <
             interval.end());
    }
  }
}
}  // namespace modle
