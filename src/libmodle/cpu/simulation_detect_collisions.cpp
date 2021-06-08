#include <absl/types/span.h>  // for Span

#include <algorithm>                                // for is_sorted, all_of, max, find_if
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/random/bernoulli_distribution.hpp>  // for bernoulli_distribution
#include <cassert>                                  // for assert
#include <cstddef>                                  // for size_t
#include <limits>                                   // for numeric_limits
#include <utility>                                  // for make_pair, pair

#include "modle/common/common.hpp"       // for bp_t, PRNG_t, MODLE_UNLIKELY, fwd, rev
#include "modle/common/utils.hpp"        // for ndebug_defined
#include "modle/dna.hpp"                 // for Chromosome
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier, NOT_OCCUPIED
#include "modle/extrusion_factors.hpp"   // for ExtrusionUnit, Lef
#include "modle/simulation.hpp"

namespace modle {

std::pair<size_t, size_t> Simulation::detect_units_at_chrom_boundaries(
    const Chromosome* chrom, absl::Span<const Lef> lefs, absl::Span<const size_t> rev_lef_ranks,
    absl::Span<const size_t> fwd_lef_ranks, absl::Span<const bp_t> rev_moves,
    absl::Span<const bp_t> fwd_moves, absl::Span<collision_t> rev_collisions,
    absl::Span<collision_t> fwd_collisions) {
  {
    assert(lefs.size() == fwd_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == rev_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == fwd_moves.size());                           // NOLINT
    assert(lefs.size() == rev_moves.size());                           // NOLINT
    assert(lefs.size() == fwd_collisions.size());                      // NOLINT
    assert(lefs.size() == rev_collisions.size());                      // NOLINT
    assert(std::is_sorted(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),  // NOLINT
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                          }));
    assert(std::is_sorted(rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                          }));
    assert(std::all_of(rev_collisions.begin(), rev_collisions.end(),  // NOLINT
                       [](const auto c) { return c == NO_COLLISION; }));
    assert(std::all_of(fwd_collisions.begin(), fwd_collisions.end(),  // NOLINT
                       [](const auto c) { return c == NO_COLLISION; }));
  }
  // Detect if the first rev unit or last fwd unit are about to fall off chrom. boundaries
  // Also detect extr. units that are already at chrom boundaries

  auto num_rev_units_at_5prime = 0UL;
  auto num_fwd_units_at_3prime = 0UL;

  const auto& first_active_fwd_unit = lefs[fwd_lef_ranks[0]].fwd_unit;
  const auto& last_active_rev_unit =
      lefs[*std::find_if(rev_lef_ranks.rbegin(), rev_lef_ranks.rend(), [&](const auto i) {
        return lefs[i].is_bound();
      })].rev_unit;

  // Detect and count the number of rev units located at the 5'-end
  for (auto i = 0UL; i < lefs.size(); ++i) {
    const auto& rev_idx = rev_lef_ranks[i];
    const auto& rev_unit = lefs[rev_idx].rev_unit;
    auto& rev_move = rev_moves[rev_idx];

    assert(lefs[rev_idx].is_bound());                         // NOLINT
    assert(chrom->start_pos() + rev_move <= rev_unit.pos());  // NOLINT

    if (rev_unit.pos() == chrom->start_pos()) {  // Unit already at the 5'-end
      assert(rev_moves[rev_idx] == 0);           // NOLINT
      ++num_rev_units_at_5prime;
      rev_collisions[rev_idx] = REACHED_CHROM_BOUNDARY;

    } else if (rev_unit.pos() > first_active_fwd_unit.pos()) {
      // As the rev_unit is located downstream of a fwd unit, a LEF-LEF collision might take place
      // before the rev unit reaches the 5'-end
      break;

    } else if (rev_unit.pos() - rev_move == chrom->start_pos()) {
      // Unit will reach the 5'-end by the end of the current epoch
      rev_collisions[rev_idx] = REACHED_CHROM_BOUNDARY;
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

    assert(fwd_unit.pos() + fwd_move < chrom->end_pos());  // NOLINT
    if (fwd_unit.pos() == chrom->end_pos() - 1) {
      assert(fwd_moves[fwd_idx] == 0);  // NOLINT
      ++num_fwd_units_at_3prime;
      fwd_collisions[fwd_idx] = REACHED_CHROM_BOUNDARY;

    } else if (fwd_unit.pos() < last_active_rev_unit.pos()) {
      break;

    } else if (fwd_unit.pos() + fwd_move == chrom->end_pos() - 1) {
      fwd_collisions[fwd_idx] = REACHED_CHROM_BOUNDARY;
      ++num_fwd_units_at_3prime;
      break;
    }
  }

  return std::make_pair(num_rev_units_at_5prime, num_fwd_units_at_3prime);
}

void Simulation::detect_lef_bar_collisions(
    const absl::Span<const Lef> lefs, const absl::Span<const size_t> rev_lef_ranks,
    const absl::Span<const size_t> fwd_lef_ranks, const absl::Span<const bp_t> rev_moves,
    const absl::Span<const bp_t> fwd_moves, const absl::Span<const ExtrusionBarrier> extr_barriers,
    const boost::dynamic_bitset<>& barrier_mask, const absl::Span<collision_t> rev_collisions,
    const absl::Span<collision_t> fwd_collisions, modle::PRNG_t& rand_eng,
    size_t num_rev_units_at_5prime, size_t num_fwd_units_at_3prime) const
    noexcept(utils::ndebug_defined()) {
  {
    assert(lefs.size() == fwd_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == rev_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == fwd_moves.size());                           // NOLINT
    assert(lefs.size() == rev_moves.size());                           // NOLINT
    assert(lefs.size() == fwd_collisions.size());                      // NOLINT
    assert(lefs.size() == rev_collisions.size());                      // NOLINT
    assert(barrier_mask.size() == extr_barriers.size());               // NOLINT
    assert(lefs.size() >= num_rev_units_at_5prime);                    // NOLINT
    assert(lefs.size() >= num_fwd_units_at_3prime);                    // NOLINT
    assert(std::is_sorted(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),  // NOLINT
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                          }));
    assert(std::is_sorted(rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
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
  auto j1 = num_rev_units_at_5prime == 0UL ? 0 : num_rev_units_at_5prime - 1;
  size_t j2 = 0;
  const auto j2_end = lefs.size() - std::min(num_fwd_units_at_3prime, num_fwd_units_at_3prime - 1);

  auto rev_idx = rev_lef_ranks[j1];
  auto fwd_idx = fwd_lef_ranks[j2];
  auto rev_unit_pos = lefs[rev_idx].rev_unit.pos();
  auto fwd_unit_pos = lefs[fwd_idx].fwd_unit.pos();

  // Loop over extr. barriers and find the first, possibly colliding extr. unit
  for (auto i = 0UL; i < extr_barriers.size(); ++i) {
    if (barrier_mask[i] == CTCF::NOT_OCCUPIED) {  // Extrusion barriers that are not occupied are
      continue;                                   // transparent to extr. units
    }

    const auto& barrier = extr_barriers[i];

    if (j1 < lefs.size()) {  // Process rev unit
      // Probability of block is set based on the extr. barrier blocking direction
      const auto& pblock = barrier.blocking_direction_major() == dna::rev
                               ? this->lef_hard_collision_pblock
                               : this->lef_soft_collision_pblock;

      // Look for the first rev extr. unit that comes after the current barrier
      while (rev_unit_pos <= barrier.pos()) {
        if (++j1 >= lefs.size()) {  // All rev units have been processed
          goto process_fwd_unit;    // Move to the next section
        }

        // Update idx and position with those corresponding to the rev extr. unit that comes next
        rev_idx = rev_lef_ranks[j1];  // in 5'-3' order
        rev_unit_pos = lefs[rev_idx].rev_unit.pos();
      }

      if (lefs[rev_idx].is_bound()) {
        // We have a LEF-BAR collision event if the distance between the rev. unit and the extr.
        // barrier is less or equal than the distance that the rev extr. unit is set to move in
        // the current iteration. If pblock != 1, then we also require a successful bernoulli
        // trial before calling a collision
        const auto delta = rev_unit_pos - barrier.pos();
        if (delta > 0 && delta <= rev_moves[rev_idx] &&
            (pblock == 1.0 || modle::bernoulli_trial{pblock}(rand_eng))) {
          // Collision detected. Assign barrier idx to the respective entry in the collision mask
          rev_collisions[rev_idx] = i;
          // Move LEF close to the extr. barrier (i.e 1bp upstream of the extr. barrier)
          // rev_move = delta > 1 ? delta - 1 : 0;
        }
      }
    }

  // Look in the previous section for detailed comments
  process_fwd_unit:
    if (j2 < j2_end) {
      const auto& pblock = barrier.blocking_direction_major() == dna::fwd
                               ? this->lef_hard_collision_pblock
                               : this->lef_soft_collision_pblock;
      // Look for the next fwd unit that comes strictly before the current extr. barrier
      while (fwd_unit_pos < barrier.pos()) {
        if (++j2 >= j2_end) {
          goto end_of_loop;
        }

        fwd_idx = fwd_lef_ranks[j2];
        fwd_unit_pos = lefs[fwd_idx].fwd_unit.pos();
      }

      // Decrement j2 by one (if it is legal to do so), so that j2 corresponds to the index of the
      // fwd extr. unit that is located as close as possible to the extr. barrier that is being
      // processed
      fwd_idx = fwd_lef_ranks[std::min(j2, j2 - 1)];
      fwd_unit_pos = lefs[fwd_idx].fwd_unit.pos();

      if (lefs[fwd_idx].is_bound()) {
        const auto delta = barrier.pos() - fwd_unit_pos;
        if (delta > 0 && delta <= fwd_moves[fwd_idx] &&
            (pblock == 1.0 || modle::bernoulli_trial{pblock}(rand_eng))) {
          fwd_collisions[fwd_idx] = i;
          // fwd_move = delta > 1 ? delta - 1 : 0;
        }
      }
    }

  end_of_loop:
    // Return immediately if all extr. units have been processed (regardless of whether there are
    // still extr. barriers to be processed)
    if (j1 == lefs.size() && j2 == j2_end) {
      return;
    }
  }
}

void Simulation::detect_primary_lef_lef_collisions(
    const absl::Span<const Lef> lefs, const absl::Span<const ExtrusionBarrier> barriers,
    const absl::Span<const size_t> rev_lef_ranks, const absl::Span<const size_t> fwd_lef_ranks,
    const absl::Span<const bp_t> rev_moves, const absl::Span<const bp_t> fwd_moves,
    const absl::Span<collision_t> rev_collisions, const absl::Span<collision_t> fwd_collisions,
    PRNG_t& rand_eng, size_t num_rev_units_at_5prime, size_t num_fwd_units_at_3prime) const
    noexcept(utils::ndebug_defined()) {
  {
    assert(lefs.size() == fwd_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == rev_lef_ranks.size());                       // NOLINT
    assert(lefs.size() == fwd_moves.size());                           // NOLINT
    assert(lefs.size() == rev_moves.size());                           // NOLINT
    assert(lefs.size() == fwd_collisions.size());                      // NOLINT
    assert(lefs.size() == rev_collisions.size());                      // NOLINT
    assert(lefs.size() >= num_rev_units_at_5prime);                    // NOLINT
    assert(lefs.size() >= num_fwd_units_at_3prime);                    // NOLINT
    assert(std::is_sorted(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),  // NOLINT
                          [&](const auto r1, const auto r2) {
                            return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                          }));
    assert(std::is_sorted(rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
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
  //      outcome of a bernoulli trial with probability of success equal to 1 - the prob. of bypass
  //    - If a LEF-LEF collision caused by two extr. unit moving in opposite direction is detected,
  //      encode the index of the extr. unit that caused the collision in the appropriate collision
  //      mask. This kind of collisions are encoded as offset + i, where offset = nbarriers and i =
  //      the index of the unit that is colliding.

  if (num_rev_units_at_5prime == lefs.size() || num_fwd_units_at_3prime == lefs.size())
    MODLE_UNLIKELY { return; }

  //    Initialize indexes so that we skip over rev units at the 5' and fwd units at the 3' (if any)
  auto i1 = 0UL;
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
      if (++j1 == j2) MODLE_UNLIKELY {
          return;  // all rev units have been processed
        }
      rev_idx = rev_lef_ranks[j1];
      rev_pos = lefs[rev_idx].rev_unit.pos();
    }

    // Find the last fwd unit that comes right before the jth rev unit
    // This is necessary in order to handle ties in fwd units ranking, as well as the case where
    // there are many fwd units between the pos of the jth-1 and jth rev extr. units
    while (fwd_pos < rev_pos) {
      if (++i1 == i2) MODLE_UNLIKELY {
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
    if (const auto delta = rev_pos - fwd_pos;
        delta > 0 && delta < rev_moves[rev_idx] + fwd_moves[fwd_idx] &&
        (this->probability_of_extrusion_unit_bypass == 0 ||
         boost::random::bernoulli_distribution{1.0 - this->probability_of_extrusion_unit_bypass}(
             rand_eng))) {
      // Declare few aliases to reduce code verbosity later on
      const auto& rev_move = rev_moves[rev_idx];
      const auto& fwd_move = fwd_moves[fwd_idx];
      auto& cause_of_collision_rev = rev_collisions[rev_idx];
      auto& cause_of_collision_fwd = fwd_collisions[fwd_idx];

      // Note that the rev collision pos will always be 1bp upstream of the fwd collision pos
      auto [collision_pos_rev, collision_pos_fwd] = compute_lef_lef_collision_pos(
          lefs[rev_idx].rev_unit, lefs[fwd_idx].fwd_unit, rev_move, fwd_move);

      if (cause_of_collision_rev == NO_COLLISION && cause_of_collision_fwd == NO_COLLISION) {
        // In the simplest case both units are free to move.
        // Thus we update their respective moves so that after calling Simulate::extrude, the two
        // units will be locate at their respective collision sites
        cause_of_collision_rev = barriers.size() + fwd_idx;
        cause_of_collision_fwd = barriers.size() + rev_idx;

      } else if (cause_of_collision_rev != NO_COLLISION && cause_of_collision_fwd == NO_COLLISION) {
        // In this case only the fwd unit is free to move.
        // This case is a bit more complicated than the first one, because we have to handle the
        // edge case where we have mistakenly predicted a LEF-BAR collision. This usually happens
        // when two extr. units are moving in opposite directions, and the unit extruding in rev
        // direction is located downstream of the unit extruding in fwd direction. If located
        // between the two units there's an extr. barrier, and one of the two units manages to
        // bypass the extr. barrier, then there's a chance that the LEF-LEF collision will take
        // place before the LEF-BAR collision (whether this happens or not depends on the distance
        // between the extr. units and the extr. barriers, as well as their respective extr. speed

        assert(cause_of_collision_rev < barriers.size() /* is LEF-BAR collision*/);  // NOLINT
        assert(collision_pos_rev != 0 && collision_pos_fwd != 0);                    // NOLINT
        const auto& barrier_pos = barriers[cause_of_collision_rev].pos();
        if (collision_pos_fwd > barrier_pos) {
          // Detected the mis-prediction mentioned above: make the LEF-BAR collision a LEF-LEF
          // collision
          cause_of_collision_rev = barriers.size() + fwd_idx;
          cause_of_collision_fwd = barriers.size() + rev_idx;
        } else {
          // fwd extr unit is being blocked by a rev unit that is itself being blocked by an extr.
          // barrier
          cause_of_collision_fwd = barriers.size() + rev_idx;
        }
        // This branch follows the same logic as the previous one. In this case the rev unit is free
        // to move, while the fwd unit has been predicted to be stalled by an extr. barrier
      } else if (cause_of_collision_rev == NO_COLLISION && cause_of_collision_fwd != NO_COLLISION) {
        assert(cause_of_collision_fwd < barriers.size() /* is LEF-BAR collision*/);  // NOLINT
        assert(collision_pos_rev != 0 && collision_pos_fwd != 0);                    // NOLINT
        const auto& barrier_pos = barriers[cause_of_collision_fwd].pos();
        if (collision_pos_rev < barrier_pos) {
          cause_of_collision_rev = barriers.size() + fwd_idx;
          cause_of_collision_fwd = barriers.size() + rev_idx;
        } else {
          cause_of_collision_rev = barriers.size() + fwd_idx;
        }
      }
    }
  }
}

void Simulation::process_secondary_lef_lef_collisions(
    const Chromosome* chrom, const absl::Span<const Lef> lefs, const size_t nbarriers,
    const absl::Span<const size_t> rev_lef_ranks, const absl::Span<const size_t> fwd_lef_ranks,
    const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
    const absl::Span<collision_t> rev_collisions, const absl::Span<collision_t> fwd_collisions,
    PRNG_t& rand_eng, size_t num_rev_units_at_5prime, size_t num_fwd_units_at_3prime) const
    noexcept(utils::ndebug_defined()) {
  {
    assert(lefs.size() == fwd_lef_ranks.size());     // NOLINT
    assert(lefs.size() == rev_lef_ranks.size());     // NOLINT
    assert(lefs.size() == fwd_moves.size());         // NOLINT
    assert(lefs.size() == rev_moves.size());         // NOLINT
    assert(lefs.size() == fwd_collisions.size());    // NOLINT
    assert(lefs.size() == rev_collisions.size());    // NOLINT
    assert(lefs.size() >= num_rev_units_at_5prime);  // NOLINT
    assert(lefs.size() >= num_fwd_units_at_3prime);  // NOLINT
    (void)chrom;
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

  const auto offset = nbarriers + lefs.size();
  const auto num_active_rev_units = std::max(1UL, num_rev_units_at_5prime);
  const auto num_active_fwd_units =
      lefs.size() - std::min(num_fwd_units_at_3prime, num_fwd_units_at_3prime - 1);

  // Loop over rev units in 5'-3' direction, and skip over units that are located at the 5'-end
  for (auto i = num_active_rev_units; i < lefs.size(); ++i) {
    // index of the ith-1 rev unit in 5'-3' order (i.e. U1)
    const auto& rev_idx1 = rev_lef_ranks[i - 1];
    const auto& rev_pos1 = lefs[rev_idx1].rev_unit.pos();
    if (rev_collisions[rev_idx1] == NO_COLLISION) {
      // If U1 is not stalled, then it is not possible to have a secondary LEF-LEF collision
      continue;
    }

    // index of the ith rev unit in 5'-3' order (i.e. U2)
    const auto& rev_idx2 = rev_lef_ranks[i];
    const auto& rev_pos2 = lefs[rev_idx2].rev_unit.pos();
    if (rev_collisions[rev_idx2] != NO_COLLISION) {
      // LEF-BAR and primary LEF-LEF collisions have precedence over secondary LEF-LEF collisions
      continue;
    }

    const auto& move1 = rev_moves[rev_idx1];
    auto& move2 = rev_moves[rev_idx2];

    assert(rev_pos2 >= rev_pos1);                    // NOLINT
    assert(rev_pos1 >= move1);                       // NOLINT
    assert(rev_pos2 >= move2);                       // NOLINT
    assert(rev_pos1 - move1 >= chrom->start_pos());  // NOLINT
    assert(rev_pos2 - move2 >= chrom->start_pos());  // NOLINT

    if (rev_pos2 - move2 <= rev_pos1 - move1 &&
        (this->probability_of_extrusion_unit_bypass == 0 ||
         boost::random::bernoulli_distribution{1.0 - this->probability_of_extrusion_unit_bypass}(
             rand_eng))) {
      rev_collisions[rev_idx2] = offset + rev_idx1;
      const auto move = rev_pos2 - (rev_pos1 - move1);
      move2 = move > 0UL ? move - 1UL : 0UL;
    }
  }

  // Loop over fwd units in 3'-5' direction, and skip over units that are located at the 3'-end.
  // Look in the section above for detailed comments
  for (auto i = num_active_fwd_units - 1; i > 0; --i) {
    // index of the ith fwd unit in 3'-5' order (i.e. U2)
    const auto& fwd_idx2 = fwd_lef_ranks[i];
    const auto& fwd_pos2 = lefs[fwd_idx2].fwd_unit.pos();
    if (fwd_collisions[fwd_idx2] == NO_COLLISION) {
      // If U2 is not stalled, then it is not possible to have a secondary LEF-LEF collision
      continue;
    }

    const auto& fwd_idx1 = fwd_lef_ranks[i - 1];
    const auto& fwd_pos1 = lefs[fwd_idx1].fwd_unit.pos();
    if (fwd_collisions[fwd_idx1] != NO_COLLISION) {
      continue;
    }

    auto& move1 = fwd_moves[fwd_idx1];
    const auto& move2 = fwd_moves[fwd_idx2];

    assert(fwd_pos2 >= fwd_pos1);                 // NOLINT
    assert(fwd_pos1 + move1 < chrom->end_pos());  // NOLINT
    assert(fwd_pos2 + move2 < chrom->end_pos());  // NOLINT

    if (fwd_pos1 + move1 >= fwd_pos2 + move2 &&
        (this->probability_of_extrusion_unit_bypass == 0 ||
         boost::random::bernoulli_distribution{1.0 - this->probability_of_extrusion_unit_bypass}(
             rand_eng))) {
      fwd_collisions[fwd_idx1] = offset + fwd_idx2;
      const auto move = (fwd_pos2 + move2) - fwd_pos1;
      move1 = move > 0UL ? move - 1UL : 0UL;
    }
  }
}

}  // namespace modle
