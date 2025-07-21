// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>
#include <fmt/format.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/extrusion_factors.hpp"
#include "modle/simulation.hpp"

namespace modle::libmodle::test {

inline constexpr auto DEFAULT_PRNG{random::PRNG(10556020843759504871ULL)};

using CollisionT = Collision<u32f>;
constexpr auto NO_COLLISION = CollisionT::NO_COLLISION;
constexpr auto COLLISION = CollisionT::COLLISION;
constexpr auto CHROM_BOUNDARY = COLLISION | CollisionT::CHROM_BOUNDARY;
constexpr auto LEF_BAR = COLLISION | CollisionT::LEF_BAR;
constexpr auto LEF_LEF_PRIMARY = COLLISION | CollisionT::LEF_LEF_PRIMARY;
constexpr auto LEF_LEF_SECONDARY = COLLISION | CollisionT::LEF_LEF_SECONDARY;

constexpr auto C = COLLISION;

template <class I>
[[nodiscard]] inline Lef construct_lef(I p1, I p2, usize binding_epoch = 0) {
  return Lef{binding_epoch, ExtrusionUnit{static_cast<bp_t>(p1)},
             ExtrusionUnit{static_cast<bp_t>(p2)}};
}

template <class BPCollection, class CollisionCollection>
[[maybe_unused]] inline void print_debug_info(
    usize i, const BPCollection& rev_moves, const BPCollection& fwd_moves,
    const BPCollection& rev_moves_expected, const BPCollection& fwd_moves_expected,
    const CollisionCollection& rev_collisions, const CollisionCollection& rev_collisions_expected,
    const CollisionCollection& fwd_collisions, const CollisionCollection& fwd_collisions_expected) {
  static_assert(std::is_same_v<typename BPCollection::value_type, bp_t>);
  fmt::print(stderr, "i={}; rev_move={}/{}; fwd_move={}/{};\n", i, rev_moves[i],
             rev_moves_expected[i], fwd_moves[i], fwd_moves_expected[i]);
  fmt::print(
      stderr, "i={}; rev_status: expected {:l} got {:l}]; fwd_status: expected {:l} got {:l}\n", i,
      rev_collisions[i], rev_collisions_expected[i], fwd_collisions[i], fwd_collisions_expected[i]);
}

template <class CollisionCollection>
[[maybe_unused]] inline void print_debug_info(usize i, const CollisionCollection& rev_collisions,
                                              const CollisionCollection& rev_collisions_expected,
                                              const CollisionCollection& fwd_collisions,
                                              const CollisionCollection& fwd_collisions_expected) {
  fmt::print(stderr, "i={}; rev_status=[{:l}\t{:l}]; fwd_status=[{:l}\t{:l}];\n", i,
             rev_collisions[i], rev_collisions_expected[i], fwd_collisions[i],
             fwd_collisions_expected[i]);
}

template <class LefCollection, class BPCollection, class CollisionCollection>
[[maybe_unused]] inline void check_simulation_result(
    const LefCollection& lefs, const BPCollection& rev_moves, const BPCollection& fwd_moves,
    const BPCollection& rev_moves_expected, const BPCollection& fwd_moves_expected,
    const CollisionCollection& rev_collisions, const CollisionCollection& rev_collisions_expected,
    const CollisionCollection& fwd_collisions, const CollisionCollection& fwd_collisions_expected,
    bool print_debug_info_ = false) {
  static_assert(std::is_same_v<typename LefCollection::value_type, Lef>);
  static_assert(std::is_same_v<typename BPCollection::value_type, bp_t>);
  for (usize i = 0; i < lefs.size(); ++i) {
    CHECK(rev_collisions[i] == rev_collisions_expected[i]);
    CHECK(fwd_collisions[i] == fwd_collisions_expected[i]);
    CHECK(rev_moves[i] == rev_moves_expected[i]);
    CHECK(fwd_moves[i] == fwd_moves_expected[i]);

    if (print_debug_info_) {
      print_debug_info(i, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                       rev_collisions, rev_collisions_expected, fwd_collisions,
                       fwd_collisions_expected);
    }
  }
}

template <class LefCollection, class CollisionCollection>
[[maybe_unused]] inline void check_collisions(const LefCollection& lefs,
                                              const CollisionCollection& rev_collisions,
                                              const CollisionCollection& rev_collisions_expected,
                                              const CollisionCollection& fwd_collisions,
                                              const CollisionCollection& fwd_collisions_expected,
                                              bool print_debug_info_ = false) {
  static_assert(std::is_same_v<typename LefCollection::value_type, Lef>);
  for (usize i = 0; i < lefs.size(); ++i) {
    CHECK(rev_collisions[i] == rev_collisions_expected[i]);
    CHECK(fwd_collisions[i] == fwd_collisions_expected[i]);

    if (print_debug_info_) {
      print_debug_info(i, rev_collisions, rev_collisions_expected, fwd_collisions,
                       fwd_collisions_expected);
    }
  }
}

template <class LefCollection, class UsizeCollection>
inline void check_that_lefs_are_sorted_by_idx(const LefCollection& lefs,
                                              const UsizeCollection& rev_ranks,
                                              const UsizeCollection& fwd_ranks) {
  static_assert(std::is_same_v<typename LefCollection::value_type, Lef>);
  static_assert(std::is_same_v<typename UsizeCollection::value_type, usize>);
  CHECK(std::is_sorted(fwd_ranks.begin(), fwd_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
  }));
  CHECK(std::is_sorted(rev_ranks.begin(), rev_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
  }));
}

template <class LefCollection, class UsizeCollection>
inline void require_that_lefs_are_sorted_by_idx(const LefCollection& lefs,
                                                const UsizeCollection& rev_ranks,
                                                const UsizeCollection& fwd_ranks) {
  static_assert(std::is_same_v<typename LefCollection::value_type, Lef>);
  static_assert(std::is_same_v<typename UsizeCollection::value_type, usize>);
  REQUIRE(std::is_sorted(fwd_ranks.begin(), fwd_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
  }));
  REQUIRE(std::is_sorted(rev_ranks.begin(), rev_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
  }));
}

[[nodiscard]] inline modle::Config init_config(usize rev_extrusion_speed, usize fwd_extrusion_speed,
                                               double rev_extrusion_speed_std = 0,
                                               double fwd_extrusion_speed_std = 0) noexcept {
  modle::Config c;

  c.rev_extrusion_speed = rev_extrusion_speed;
  c.rev_extrusion_speed_std = rev_extrusion_speed_std;
  c.fwd_extrusion_speed = fwd_extrusion_speed;
  c.fwd_extrusion_speed_std = fwd_extrusion_speed_std;
  c.probability_of_extrusion_unit_bypass = 0;
  c.lef_bar_major_collision_pblock = 1;
  c.lef_bar_minor_collision_pblock = 0;

  return c;
}

[[nodiscard]] inline GenomicInterval init_interval(
    std::string name, bp_t chrom_size, bp_t chrom_start = 0,
    bp_t chrom_end = (std::numeric_limits<bp_t>::max)(), bp_t resolution = 1, bp_t diag_width = 1) {
  chrom_end = std::min(chrom_end, chrom_size);
  assert(chrom_start < chrom_end);
  auto chrom = std::make_shared<Chromosome>(0, std::move(name), chrom_size);
  return {0, chrom, chrom_start, chrom_end, resolution, diag_width};
}

// This function is here as a compatibility layer between the old and new way to construct extrusion
// barriers, and should be removed at some point in the future.
template <ExtrusionBarriers::State state = ExtrusionBarriers::State::ACTIVE, class... Args>
[[nodiscard]] inline ExtrusionBarriers construct_barriers(Args... barriers) {
  ExtrusionBarriers barriers_;
  for (const auto& barrier : {barriers...}) {
    barriers_.push_back(barrier, state);
  }
  return barriers_;
}

}  // namespace modle::libmodle::test
