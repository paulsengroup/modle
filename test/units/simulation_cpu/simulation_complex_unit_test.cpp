// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/types/span.h>  // for MakeSpan

#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cassert>                                  // for assert
#include <catch2/catch_test_macros.hpp>
#include <memory>       // for allocator, allocator_traits<>::value_...
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "./common.hpp"                        // for construct_lef, NO_COLLISION, check_si...
#include "modle/common/common.hpp"             // for bp_t, usize
#include "modle/common/random.hpp"             // for PRNG
#include "modle/common/simulation_config.hpp"  // for Config
#include "modle/extrusion_barriers.hpp"        // for ExtrusionBarrier, State, OCCUPIED
#include "modle/extrusion_factors.hpp"         // for Lef
#include "modle/genome.hpp"                    // for Chromosome
#include "modle/simulation.hpp"                // for Simulation

namespace modle::test::libmodle {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 001", "[simulation][short]") {
  const auto c = init_config(75, 75);
  const usize nlefs = 7;
  const usize nbarriers = 5;
  const Chromosome chrom{0, "chr1", 0, 1000, 1000};
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(25, 30, 0),
      construct_lef(150, 150, 1),
      construct_lef(200, 350, 2),
      construct_lef(230, 399, 3),
      construct_lef(425, 425, 4),
      construct_lef(625, 800, 5),
      construct_lef(650, 650, 6)
  };

  const auto barriers = construct_barriers(ExtrusionBarrier{100, 1.0, 0.0, '-'},
                                           ExtrusionBarrier{105, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{400, 1.0, 0.0, '-'},
                                           ExtrusionBarrier{600, 1.0, 0.0, '-'},
                                           ExtrusionBarrier{850, 1.0, 0.0, '+'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<usize, nlefs> rev_ranks{0, 1, 2, 3, 4, 5, 6};
  const std::array<usize, nlefs> fwd_ranks{0, 1, 2, 3, 4, 6, 5};

  std::array<bp_t, nlefs> rev_moves{25, 75, 75, 75, 75, 75, 75};
  std::array<bp_t, nlefs> fwd_moves{75, 75, 75, 75, 75, 75, 75};

  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{5, CHROM_BOUNDARY},
                                                              CollisionT{1, LEF_BAR},
                                                              CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{2, LEF_LEF_SECONDARY},
                                                              CollisionT{3, LEF_LEF_PRIMARY},
                                                              CollisionT{},
                                                              CollisionT{}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{0, LEF_BAR},
                                                              CollisionT{2, LEF_LEF_PRIMARY},
                                                              CollisionT{3, LEF_LEF_SECONDARY},
                                                              CollisionT{2, LEF_BAR},
                                                              CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{}};
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{25, 44, 25, 54, 25, 75, 75};
  const std::array<bp_t, nlefs> fwd_moves_expected{69, 24, 48,  0, 75, 75, 75};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_and_clamp_moves(chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 002", "[simulation][short]") {
  const auto c = init_config(75, 75);
  const usize nlefs = 7;
  const usize nbarriers = 5;
  const Chromosome chrom{0, "chr1", 0, 1000, 1000};
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(200, 375, 0),
      construct_lef(350, 350, 1),
      construct_lef(575, 575, 2),
      construct_lef(601, 770, 3),
      construct_lef(650, 800, 4),
      construct_lef(850, 850, 5),
      construct_lef(970, 975, 6)
  };

  const auto barriers = construct_barriers(ExtrusionBarrier{150, 1.0, 0.0, '-'},
                                           ExtrusionBarrier{400, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{600, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{895, 1.0, 0.0, '-'},
                                           ExtrusionBarrier{900, 1.0, 0.0, '+'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<usize, nlefs> rev_ranks{0, 1, 2, 3, 4, 5, 6};
  const std::array<usize, nlefs> fwd_ranks{1, 0, 2, 3, 4, 5, 6};

  std::array<bp_t, nlefs> rev_moves{75, 75, 75, 75, 75, 75, 75};
  std::array<bp_t, nlefs> fwd_moves{75, 75, 75, 75, 75, 75, 24};

  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{2, LEF_BAR},
                                                              CollisionT{3, LEF_LEF_SECONDARY},
                                                              CollisionT{4, LEF_LEF_PRIMARY},
                                                              CollisionT{4, LEF_BAR}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{3, LEF_LEF_PRIMARY},
                                                              CollisionT{4, LEF_LEF_SECONDARY},
                                                              CollisionT{5, LEF_LEF_PRIMARY},
                                                              CollisionT{3, LEF_BAR},
                                                              CollisionT{3, CHROM_BOUNDARY}};
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  // clang-format off
  const std::array<bp_t, nlefs> rev_moves_expected{75, 75, 75,  0, 48, 25, 69};
  const std::array<bp_t, nlefs> fwd_moves_expected{75, 75, 25, 53, 24, 44, 24};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_and_clamp_moves(chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 003 - Soft collisions on", "[simulation][short]") {
  auto c = init_config(75, 75);
  c.lef_bar_major_collision_pblock = 1;
  c.lef_bar_minor_collision_pblock = 1;
  const usize nlefs = 7;
  const usize nbarriers = 5;
  const Chromosome chrom{0, "chr1", 0, 1000, 1000};
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(200, 375, 0),
      construct_lef(350, 350, 1),
      construct_lef(575, 575, 2),
      construct_lef(601, 770, 3),
      construct_lef(650, 800, 4),
      construct_lef(850, 850, 5),
      construct_lef(970, 975, 6)
  };

  const auto barriers = construct_barriers(ExtrusionBarrier{150, 1.0, 0.0, '-'},
                                           ExtrusionBarrier{400, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{600, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{895, 1.0, 0.0, '-'},
                                           ExtrusionBarrier{900, 1.0, 0.0, '+'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<usize, nlefs> rev_ranks{0, 1, 2, 3, 4, 5, 6};
  const std::array<usize, nlefs> fwd_ranks{1, 0, 2, 3, 4, 5, 6};

  std::array<bp_t, nlefs> rev_moves{75, 75, 75, 75, 75, 75, 75};
  std::array<bp_t, nlefs> fwd_moves{75, 75, 75, 75, 75, 75, 24};

  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{0, LEF_BAR},
                                                              CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{2, LEF_BAR},
                                                              CollisionT{3, LEF_LEF_SECONDARY},
                                                              CollisionT{4, LEF_LEF_PRIMARY},
                                                              CollisionT{4, LEF_BAR}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{1, LEF_BAR},
                                                              CollisionT{0, LEF_LEF_SECONDARY},
                                                              CollisionT{2, LEF_BAR},
                                                              CollisionT{4, LEF_LEF_SECONDARY},
                                                              CollisionT{5, LEF_LEF_PRIMARY},
                                                              CollisionT{3, LEF_BAR},
                                                              CollisionT{3, CHROM_BOUNDARY}};
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{49, 75, 75,  0, 48, 25, 69};
  const std::array<bp_t, nlefs> fwd_moves_expected{24, 48, 24, 53, 24, 44, 24};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_and_clamp_moves(chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 004 - Inactive barriers", "[simulation][short]") {
  const auto c = init_config(75, 75);
  const usize nlefs = 7;
  const usize nbarriers = 5;
  const Chromosome chrom{0, "chr1", 0, 1000, 1000};
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(200, 375, 0),
      construct_lef(350, 350, 1),
      construct_lef(575, 575, 2),
      construct_lef(601, 770, 3),
      construct_lef(650, 800, 4),
      construct_lef(850, 850, 5),
      construct_lef(970, 975, 6)
  };

  constexpr auto s =ExtrusionBarriers::State::INACTIVE;
  const auto barriers = construct_barriers<s>(ExtrusionBarrier{150, 1.0, 0.0, '-'},
                                              ExtrusionBarrier{400, 1.0, 0.0, '+'},
                                              ExtrusionBarrier{600, 1.0, 0.0, '+'},
                                              ExtrusionBarrier{895, 1.0, 0.0, '-'},
                                              ExtrusionBarrier{900, 1.0, 0.0, '+'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<usize, nlefs> rev_ranks{0, 1, 2, 3, 4, 5, 6};
  const std::array<usize, nlefs> fwd_ranks{1, 0, 2, 3, 4, 5, 6};

  std::array<bp_t, nlefs> rev_moves{75, 75, 75, 75, 75, 75, 75};
  std::array<bp_t, nlefs> fwd_moves{75, 75, 75, 75, 75, 75, 24};

  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{2, LEF_LEF_PRIMARY},
                                                              CollisionT{3, LEF_LEF_SECONDARY},
                                                              CollisionT{4, LEF_LEF_PRIMARY},
                                                              CollisionT{5, LEF_LEF_PRIMARY}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{3, LEF_LEF_PRIMARY},
                                                              CollisionT{4, LEF_LEF_SECONDARY},
                                                              CollisionT{5, LEF_LEF_PRIMARY},
                                                              CollisionT{6, LEF_LEF_PRIMARY},
                                                              CollisionT{3, CHROM_BOUNDARY}};
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{75, 75, 75, 13, 61, 25, 60};
  const std::array<bp_t, nlefs> fwd_moves_expected{75, 75, 12, 53, 24, 59, 24};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_and_clamp_moves(chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 005 - Multiple LEFs located at the same site", "[simulation][short]") {
  const auto c = init_config(25, 25);
  constexpr usize nlefs = 6;
  constexpr usize nbarriers = 1;
  const Chromosome chrom{0, "chr1", 0, 150, 150};
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(30, 50, 0),
      construct_lef(60, 80, 1),
      construct_lef(60, 80, 2),
      construct_lef(65, 125, 3),
      construct_lef(140, 140, 4),
      construct_lef(140, 140, 5)
  };

  const auto barriers = construct_barriers(ExtrusionBarrier{100, 1.0, 0.0, '-'});
  REQUIRE(barriers.size() == nbarriers);

  const std::array<usize, nlefs> rev_ranks{0, 1, 2, 3, 4, 5};
  const std::array<usize, nlefs> fwd_ranks{0, 1, 2, 3, 4, 5};

  std::array<bp_t, nlefs> rev_moves{25, 25, 25, 25, 25, 25};
  std::array<bp_t, nlefs> fwd_moves{25, 25, 25, 24,  8,  9};

  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{0, LEF_LEF_PRIMARY},
                                                              CollisionT{1, LEF_LEF_SECONDARY},
                                                              CollisionT{2, LEF_LEF_SECONDARY},
                                                              CollisionT{3, LEF_LEF_PRIMARY},
                                                              CollisionT{4, LEF_LEF_SECONDARY}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{2, LEF_LEF_SECONDARY},
                                                              CollisionT{0, LEF_BAR},
                                                              CollisionT{4, LEF_LEF_PRIMARY},
                                                              CollisionT{5, LEF_LEF_SECONDARY},
                                                              CollisionT{3, CHROM_BOUNDARY}};
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{25,  5,  4, 8, 8, 7};
  const std::array<bp_t, nlefs> fwd_moves_expected{ 4, 18, 19, 6, 8, 9};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_and_clamp_moves(chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 006 - Few inactive LEFs", "[simulation][short]") {
  const auto c = init_config(25, 25);
  constexpr usize nlefs = 6;
  constexpr usize nbarriers = 1;
  const Chromosome chrom{0, "chr1", 0, 150, 150};
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(30, 50, 0),
      construct_lef(60, 80, 1),
      construct_lef(60, 80, 2),
      construct_lef(65, 125, 3),
      construct_lef(140, 140, 4),
      construct_lef(140, 140, 5)
  };

  lefs[2].release();
  lefs[5].release();

  const auto barriers = construct_barriers(ExtrusionBarrier{100, 1.0, 0.0, '-'});
  REQUIRE(barriers.size() == nbarriers);

  const std::array<usize, nlefs> rev_ranks{0, 1, 3, 4, 2, 5};
  const std::array<usize, nlefs> fwd_ranks{0, 1, 3, 4, 2, 5};

  std::array<bp_t, nlefs> rev_moves{25, 25, 0, 25, 25, 0};
  std::array<bp_t, nlefs> fwd_moves{25, 25, 0, 24,  9, 0};
  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{0, LEF_LEF_PRIMARY},
                                                              CollisionT{},
                                                              CollisionT{1, LEF_LEF_SECONDARY},
                                                              CollisionT{3, LEF_LEF_PRIMARY},
                                                              CollisionT{}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{0, LEF_BAR},
                                                              CollisionT{},
                                                              CollisionT{4, LEF_LEF_PRIMARY},
                                                              CollisionT{3, CHROM_BOUNDARY},
                                                              CollisionT{}};
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{25,  5, 0, 9, 8, 0};
  const std::array<bp_t, nlefs> fwd_moves_expected{ 4, 19, 0, 6, 9, 0};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_and_clamp_moves(chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 007 - LEF-LEF collision overrides LEF-BAR collision 1",
          "[simulation][short]") {
  const auto c = init_config(20, 20);
  constexpr usize nlefs = 2;
  constexpr usize nbarriers = 1;
  const Chromosome chrom{0, "chr1", 0, 200, 200};
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(50, 95, 0),
      construct_lef(110, 150, 1)
  };

  const auto barriers = construct_barriers(ExtrusionBarrier{100, 1.0, 0.0, '+'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<usize, nlefs> rev_ranks{0, 1};
  const std::array<usize, nlefs> fwd_ranks{0, 1};

  std::array<bp_t, nlefs> rev_moves{20, 20};
  std::array<bp_t, nlefs> fwd_moves{20, 20};

  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{0, LEF_LEF_PRIMARY}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{}};
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  // clang-format off
  const std::array<bp_t, nlefs> rev_moves_expected{20,  7};
  const std::array<bp_t, nlefs> fwd_moves_expected{ 7, 20};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 008 - LEF-LEF collision overrides LEF-BAR collision 2",
          "[simulation][short]") {
  const auto c = init_config(20, 20);
  constexpr usize nlefs = 2;
  constexpr usize nbarriers = 1;
  const Chromosome chrom{0, "chr1", 0, 200, 200};
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(50, 90, 0),
      construct_lef(105, 150, 1)
  };
  // clang-format on

  const auto barriers = construct_barriers(ExtrusionBarrier{100, 1.0, 0.0, '-'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<usize, nlefs> rev_ranks{0, 1};
  const std::array<usize, nlefs> fwd_ranks{0, 1};

  std::array<bp_t, nlefs> rev_moves{20, 20};
  std::array<bp_t, nlefs> fwd_moves{20, 20};

  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{0, LEF_LEF_PRIMARY}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{}};
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  // clang-format off
  const std::array<bp_t, nlefs> rev_moves_expected{20,  7};
  const std::array<bp_t, nlefs> fwd_moves_expected{ 7, 20};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 009 - Ensure stacked LEFs do not interfere with surrounding extr. barriers",
          "[simulation][short]") {
  const auto c = init_config(10, 10);
  constexpr usize nlefs = 5;
  constexpr usize nbarriers = 2;
  const Chromosome chrom{0, "chr1", 0, 200, 200};
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(95, 100, 0),
      construct_lef(101, 103, 1),
      construct_lef(102, 110, 2),
      construct_lef(104, 111, 3),
      construct_lef(105, 112, 4)
  };

  const auto barriers = construct_barriers(ExtrusionBarrier{100, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{105, 1.0, 0.0, '-'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<usize, nlefs> rev_ranks{0, 1, 2, 3, 4};
  const std::array<usize, nlefs> fwd_ranks{0, 1, 2, 3, 4};

  std::array<bp_t, nlefs> rev_moves{10, 10, 10, 10, 10};
  std::array<bp_t, nlefs> fwd_moves{10, 10, 10, 10, 10};
  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{0, LEF_BAR},
                                                              CollisionT{1, LEF_LEF_SECONDARY},
                                                              CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{3, LEF_LEF_SECONDARY}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{3, LEF_LEF_PRIMARY},
                                                              CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{}};
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{10, 0,  0,  0,  0};
  const std::array<bp_t, nlefs> fwd_moves_expected{0,  0, 10, 10, 10};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_and_clamp_moves(chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 010 - Ensure stacked LEFs do not interfere with surrounding extr. barriers",
          "[simulation][short]") {
  const auto c = init_config(10, 10);
  constexpr usize nlefs = 6;
  constexpr usize nbarriers = 2;
  const Chromosome chrom{0, "chr1", 0, 200, 200};
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(95, 100, 0),
      construct_lef(101, 103, 1),
      construct_lef(102, 110, 2),
      construct_lef(104, 111, 3),
      construct_lef(105, 112, 4),
      construct_lef(102, 102, 5)
  };

  const auto barriers = construct_barriers(ExtrusionBarrier{100, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{105, 1.0, 0.0, '-'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<usize, nlefs> rev_ranks{0, 1, 2, 5, 3, 4};
  const std::array<usize, nlefs> fwd_ranks{0, 5, 1, 2, 3, 4};

  std::array<bp_t, nlefs> rev_moves{10, 10, 10, 10, 10, 10};
  std::array<bp_t, nlefs> fwd_moves{10, 10, 10, 10, 10, 10};
  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{0, LEF_BAR},
                                                              CollisionT{1, LEF_LEF_SECONDARY},
                                                              CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{3, LEF_LEF_SECONDARY},
                                                              CollisionT{2, LEF_LEF_SECONDARY}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{3, LEF_LEF_PRIMARY},
                                                              CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{1, LEF_LEF_SECONDARY}};
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{10, 0,  0,  0,  0, 0};
  const std::array<bp_t, nlefs> fwd_moves_expected{ 0, 0, 10, 10, 10, 0};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_and_clamp_moves(chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 011 - fix_secondary_lef_lef_collisions", "[simulation][short]") {
  auto c = init_config(3, 2);
  c.probability_of_extrusion_unit_bypass = 0.25;
  constexpr usize nlefs = 2;
  constexpr usize nbarriers = 1;
  const Chromosome chrom{0, "chr1", 0, 200, 200};
  // DO NOT CHANGE SEED!
  // The outcome of this test depends on the seed
  auto rand_eng = random::PRNG(752741483ULL);

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(25, 95, 0),
      construct_lef(50, 99, 0)
  };
  // clang-format on

  const auto barriers = construct_barriers(ExtrusionBarrier{100, 1.0, 0.0, '-'});
  REQUIRE(barriers.size() == nbarriers);

  std::array<usize, nlefs> rev_ranks{0, 1};
  std::array<usize, nlefs> fwd_ranks{0, 1};

  std::array<bp_t, nlefs> rev_moves{10, 10};
  std::array<bp_t, nlefs> fwd_moves{10, 10};
  // clang-format off
  std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                        CollisionT{}};
  std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{0, LEF_BAR},
                                                        CollisionT{1, LEF_LEF_SECONDARY}};
  // clang-format on
  // Important! The expected collisions and moves look the opposite of one would expect
  // because fix_secondary_lef_lef_collisions is swapping ranks
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  // clang-format off
  const std::array<bp_t, nlefs> rev_moves_expected{10, 10};
  const std::array<bp_t, nlefs> fwd_moves_expected{ 0,  3};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_and_clamp_moves(chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  modle::Simulation::test_fix_secondary_lef_lef_collisions(
      chrom, absl::MakeSpan(lefs), absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks),
      absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions));

  CHECK(rev_ranks[0] == 0);
  CHECK(rev_ranks[1] == 1);
  CHECK(fwd_ranks[0] == 1);
  CHECK(fwd_ranks[1] == 0);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 012 - fix_secondary_lef_lef_collisions", "[simulation][short]") {
  auto c = init_config(3, 2);
  c.probability_of_extrusion_unit_bypass = 0.25;
  constexpr usize nlefs = 2;
  constexpr usize nbarriers = 1;
  const Chromosome chrom{0, "chr1", 0, 200, 200};
  // DO NOT CHANGE SEED!
  // The outcome of this test depends on the seed
  auto rand_eng = random::PRNG(752741483ULL);

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(26, 75, 0),
      construct_lef(30, 80, 0)
  };

  const auto barriers = construct_barriers(ExtrusionBarrier{25, 1.0, 0.0, '+'});
  REQUIRE(barriers.size() == nbarriers);

  std::array<usize, nlefs> rev_ranks{0, 1};
  std::array<usize, nlefs> fwd_ranks{0, 1};

  std::array<bp_t, nlefs> rev_moves{10, 10};
  std::array<bp_t, nlefs> fwd_moves{10, 10};
  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{0, LEF_LEF_SECONDARY},
                                                              CollisionT{0, LEF_BAR}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{},
                                                              CollisionT{}};
  // Important! The expecte collisions and moves look the opposite of one would expect
  // because fix_secondary_lef_lef_collisions is swapping ranks
  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{3, 0};
  const std::array<bp_t, nlefs> fwd_moves_expected{10, 10};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_and_clamp_moves(chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, absl::MakeSpan(lefs), rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);

  modle::Simulation::test_fix_secondary_lef_lef_collisions(
      chrom, absl::MakeSpan(lefs), absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks),
      absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions));

  CHECK(rev_ranks[0] == 1);
  CHECK(rev_ranks[1] == 0);
  CHECK(fwd_ranks[0] == 0);
  CHECK(fwd_ranks[1] == 1);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

}  // namespace modle::test::libmodle
