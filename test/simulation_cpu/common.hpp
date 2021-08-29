#pragma once

#include <absl/types/span.h>
#include <fmt/format.h>

#include <algorithm>
#include <catch2/catch.hpp>
#include <cstddef>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/extrusion_factors.hpp"
#include "modle/simulation.hpp"

namespace modle::test::libmodle {

inline const auto& REACHED_CHROM_BOUNDARY = Simulation::REACHED_CHROM_BOUNDARY;
inline const auto& NO_COLLISION = Simulation::NO_COLLISION;
template <typename I>
[[nodiscard]] inline Lef construct_lef(I p1, I p2, size_t binding_epoch = 0) {
  return Lef{binding_epoch, ExtrusionUnit{static_cast<bp_t>(p1)},
             ExtrusionUnit{static_cast<bp_t>(p2)}};
}

[[maybe_unused]] inline void print_debug_info(
    size_t i, const std::vector<bp_t>& rev_moves, const std::vector<bp_t>& fwd_moves,
    const std::vector<bp_t>& rev_moves_expected, const std::vector<bp_t>& fwd_moves_expected,
    const std::vector<collision_t>& rev_collision_mask,
    const std::vector<collision_t>& rev_collision_mask_expected,
    const std::vector<collision_t>& fwd_collision_mask,
    const std::vector<collision_t>& fwd_collision_mask_expected) {
  fmt::print(stderr, FMT_STRING("i={}; rev_move={}/{}; fwd_move={}/{};\n"), i, rev_moves[i],
             rev_moves_expected[i], fwd_moves[i], fwd_moves_expected[i]);
  fmt::print(stderr, FMT_STRING("i={}; rev_status={}/{}; fwd_status={}/{};\n"), i,
             rev_collision_mask[i], rev_collision_mask_expected[i], fwd_collision_mask[i],
             fwd_collision_mask_expected[i]);
}

[[maybe_unused]] inline void check_simulation_result(
    const std::vector<Lef>& lefs, const std::vector<bp_t>& rev_moves,
    const std::vector<bp_t>& fwd_moves, const std::vector<bp_t>& rev_moves_expected,
    const std::vector<bp_t>& fwd_moves_expected, const std::vector<collision_t>& rev_collision_mask,
    const std::vector<collision_t>& rev_collision_mask_expected,
    const std::vector<collision_t>& fwd_collision_mask,
    const std::vector<collision_t>& fwd_collision_mask_expected, bool print_debug_info_ = false) {
  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
    CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
    CHECK(rev_moves[i] == rev_moves_expected[i]);
    CHECK(fwd_moves[i] == fwd_moves_expected[i]);

    if (print_debug_info_) {
      print_debug_info(i, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                       rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                       fwd_collision_mask_expected);
    }
  }
}

inline void check_that_lefs_are_sorted_by_idx(const std::vector<Lef>& lefs,
                                              const std::vector<size_t>& rev_ranks,
                                              const std::vector<size_t>& fwd_ranks) {
  CHECK(std::is_sorted(fwd_ranks.begin(), fwd_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
  }));
  CHECK(std::is_sorted(rev_ranks.begin(), rev_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
  }));
}

inline void require_that_lefs_are_sorted_by_idx(const std::vector<Lef>& lefs,
                                                const std::vector<size_t>& rev_ranks,
                                                const std::vector<size_t>& fwd_ranks) {
  REQUIRE(std::is_sorted(fwd_ranks.begin(), fwd_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
  }));
  REQUIRE(std::is_sorted(rev_ranks.begin(), rev_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
  }));
}

}  // namespace modle::test::libmodle
