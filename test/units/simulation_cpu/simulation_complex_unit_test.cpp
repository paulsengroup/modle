// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/types/span.h>  // for MakeSpan

#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cassert>                                  // for assert
#include <catch2/catch.hpp>                         // for SourceLineInfo, StringRef, TEST_CASE
#include <memory>                                   // for allocator, allocator_traits<>::value_...
#include <string_view>                              // for string_view
#include <vector>                                   // for vector

#include "./common.hpp"                  // for construct_lef, NO_COLLISION, check_si...
#include "modle/common/common.hpp"       // for bp_t, usize
#include "modle/common/config.hpp"       // for Config
#include "modle/common/random.hpp"       // for PRNG
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier, State, OCCUPIED
#include "modle/extrusion_factors.hpp"   // for Lef
#include "modle/genome.hpp"              // for Chromosome
#include "modle/simulation.hpp"          // for Simulation

namespace modle::test::libmodle {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 001", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 75;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 75;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_bar_major_collision_pblock = 1;        // NOLINT
  c.lef_bar_minor_collision_pblock = 0;        // NOLINT
  const usize nlefs = 7;
  const usize nbarriers = 5;
  const Chromosome chrom{0, "chr1", 0, 1000, 1000};
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = random::PRNG(10556020843759504871ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      construct_lef(25, 30, 0),    // NOLINT
      construct_lef(150, 150, 1),  // NOLINT
      construct_lef(200, 350, 2),  // NOLINT
      construct_lef(230, 399, 3),  // NOLINT
      construct_lef(425, 425, 4),  // NOLINT
      construct_lef(625, 800, 5),  // NOLINT
      construct_lef(650, 650, 6)   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{105, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{400, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{600, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{850, 1.0, 0.0, '+'}};
  // clang-format off
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<usize> rev_ranks{0, 1, 2, 3, 4, 5, 6};  // NOLINT
  const std::vector<usize> fwd_ranks{0, 1, 2, 3, 4, 6, 5};  // NOLINT

  std::vector<bp_t> rev_moves{25, 75, 75, 75, 75, 75, 75};  // NOLINT
  std::vector<bp_t> fwd_moves{75, 75, 75, 75, 75, 75, 75};  // NOLINT

  const std::vector<CollisionT> rev_collisions_expected{{5, CHROM_BOUNDARY},
                                                        {1, LEF_BAR},
                                                        {1, LEF_LEF_PRIMARY},
                                                        {2, LEF_LEF_SECONDARY},
                                                        {3, LEF_LEF_PRIMARY},
                                                        {},
                                                        {}};
  const std::vector<CollisionT> fwd_collisions_expected{{0, LEF_BAR},
                                                        {2, LEF_LEF_PRIMARY},
                                                        {3, LEF_LEF_SECONDARY},
                                                        {2, LEF_BAR},
                                                        {},
                                                        {},
                                                        {}};
  std::vector<CollisionT> rev_collisions(nlefs);
  std::vector<CollisionT> fwd_collisions(nlefs);

  const std::vector<bp_t> rev_moves_expected{25, 44, 25, 54, 25, 75, 75};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{69, 24, 48, 0, 75, 75, 75};   // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                     // NOLINT
  assert(barriers.size() == nbarriers);             // NOLINT
  assert(rev_collisions_expected.size() == nlefs);  // NOLINT
  assert(fwd_collisions_expected.size() == nlefs);  // NOLINT
  assert(rev_collisions.size() == nlefs);           // NOLINT
  assert(fwd_collisions.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_moves(chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
                                absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      barriers, barrier_mask, absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 002", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 75;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 75;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_bar_major_collision_pblock = 1;        // NOLINT
  c.lef_bar_minor_collision_pblock = 0;        // NOLINT
  const usize nlefs = 7;
  const usize nbarriers = 5;
  const Chromosome chrom{0, "chr1", 0, 1000, 1000};
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = random::PRNG(9668495374385482848ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      construct_lef(200, 375, 0),  // NOLINT
      construct_lef(350, 350, 1),  // NOLINT
      construct_lef(575, 575, 2),  // NOLINT
      construct_lef(601, 770, 3),  // NOLINT
      construct_lef(650, 800, 4),  // NOLINT
      construct_lef(850, 850, 5),  // NOLINT
      construct_lef(970, 975, 6)   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{150, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{400, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{600, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{895, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{900, 1.0, 0.0, '+'}};
  // clang-format off
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<usize> rev_ranks{0, 1, 2, 3, 4, 5, 6};  // NOLINT
  const std::vector<usize> fwd_ranks{1, 0, 2, 3, 4, 5, 6};  // NOLINT

  std::vector<bp_t> rev_moves{75, 75, 75, 75, 75, 75, 75};  // NOLINT
  std::vector<bp_t> fwd_moves{75, 75, 75, 75, 75, 75, 24};  // NOLINT

  const std::vector<CollisionT> rev_collisions_expected{{},
                                                        {},
                                                        {},
                                                        {2, LEF_BAR},
                                                        {3, LEF_LEF_SECONDARY},
                                                        {4, LEF_LEF_PRIMARY},
                                                        {4, LEF_BAR}};
  const std::vector<CollisionT> fwd_collisions_expected{{},
                                                        {},
                                                        {3, LEF_LEF_PRIMARY},
                                                        {4, LEF_LEF_SECONDARY},
                                                        {5, LEF_LEF_PRIMARY},
                                                        {3, LEF_BAR},
                                                        {3, CHROM_BOUNDARY}};
  std::vector<CollisionT> rev_collisions(nlefs);
  std::vector<CollisionT> fwd_collisions(nlefs);

  const std::vector<bp_t> rev_moves_expected{75, 75, 75, 0, 48, 25, 69};   // NOLINT
  const std::vector<bp_t> fwd_moves_expected{75, 75, 25, 53, 24, 44, 24};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                     // NOLINT
  assert(barriers.size() == nbarriers);             // NOLINT
  assert(rev_collisions_expected.size() == nlefs);  // NOLINT
  assert(fwd_collisions_expected.size() == nlefs);  // NOLINT
  assert(rev_collisions.size() == nlefs);           // NOLINT
  assert(fwd_collisions.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_moves(chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
                                absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      barriers, barrier_mask, absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 003 - Soft collisions on", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 75;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 75;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_bar_major_collision_pblock = 1;        // NOLINT
  c.lef_bar_minor_collision_pblock = 1;        // NOLINT
  const usize nlefs = 7;
  const usize nbarriers = 5;
  const Chromosome chrom{0, "chr1", 0, 1000, 1000};
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = random::PRNG(919341715542527390ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      construct_lef(200, 375, 0),  // NOLINT
      construct_lef(350, 350, 1),  // NOLINT
      construct_lef(575, 575, 2),  // NOLINT
      construct_lef(601, 770, 3),  // NOLINT
      construct_lef(650, 800, 4),  // NOLINT
      construct_lef(850, 850, 5),  // NOLINT
      construct_lef(970, 975, 6)   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{150, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{400, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{600, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{895, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{900, 1.0, 0.0, '+'}};
  // clang-format off
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<usize> rev_ranks{0, 1, 2, 3, 4, 5, 6};  // NOLINT
  const std::vector<usize> fwd_ranks{1, 0, 2, 3, 4, 5, 6};  // NOLINT

  std::vector<bp_t> rev_moves{75, 75, 75, 75, 75, 75, 75};  // NOLINT
  std::vector<bp_t> fwd_moves{75, 75, 75, 75, 75, 75, 24};  // NOLINT

  const std::vector<CollisionT> rev_collisions_expected{{0, LEF_BAR},
                                                        {},
                                                        {},
                                                        {2, LEF_BAR},
                                                        {3, LEF_LEF_SECONDARY},
                                                        {4, LEF_LEF_PRIMARY},
                                                        {4, LEF_BAR}};
  const std::vector<CollisionT> fwd_collisions_expected{{1, LEF_BAR},
                                                        {0, LEF_LEF_SECONDARY},
                                                        {2, LEF_BAR},
                                                        {4, LEF_LEF_SECONDARY},
                                                        {5, LEF_LEF_PRIMARY},
                                                        {3, LEF_BAR},
                                                        {3, CHROM_BOUNDARY}};
  std::vector<CollisionT> rev_collisions(nlefs);
  std::vector<CollisionT> fwd_collisions(nlefs);

  const std::vector<bp_t> rev_moves_expected{49, 75, 75, 0, 48, 25, 69};   // NOLINT
  const std::vector<bp_t> fwd_moves_expected{24, 48, 24, 53, 24, 44, 24};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                     // NOLINT
  assert(barriers.size() == nbarriers);             // NOLINT
  assert(rev_collisions_expected.size() == nlefs);  // NOLINT
  assert(fwd_collisions_expected.size() == nlefs);  // NOLINT
  assert(rev_collisions.size() == nlefs);           // NOLINT
  assert(fwd_collisions.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_moves(chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
                                absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      barriers, barrier_mask, absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 004 - Inactive barriers", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 75;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 75;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_bar_major_collision_pblock = 1;        // NOLINT
  c.lef_bar_minor_collision_pblock = 0;        // NOLINT
  const usize nlefs = 7;
  const usize nbarriers = 5;
  const Chromosome chrom{0, "chr1", 0, 1000, 1000};
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = random::PRNG(138440880292496584ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      construct_lef(200, 375, 0),  // NOLINT
      construct_lef(350, 350, 1),  // NOLINT
      construct_lef(575, 575, 2),  // NOLINT
      construct_lef(601, 770, 3),  // NOLINT
      construct_lef(650, 800, 4),  // NOLINT
      construct_lef(850, 850, 5),  // NOLINT
      construct_lef(970, 975, 6)   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{150, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{400, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{600, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{895, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{900, 1.0, 0.0, '+'}};
  // clang-format off
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::NOT_OCCUPIED));

  const std::vector<usize> rev_ranks{0, 1, 2, 3, 4, 5, 6};  // NOLINT
  const std::vector<usize> fwd_ranks{1, 0, 2, 3, 4, 5, 6};  // NOLINT

  std::vector<bp_t> rev_moves{75, 75, 75, 75, 75, 75, 75};  // NOLINT
  std::vector<bp_t> fwd_moves{75, 75, 75, 75, 75, 75, 24};  // NOLINT

  const std::vector<CollisionT> rev_collisions_expected{{},
                                                        {},
                                                        {},
                                                        {2, LEF_LEF_PRIMARY},
                                                        {3, LEF_LEF_SECONDARY},
                                                        {4, LEF_LEF_PRIMARY},
                                                        {5, LEF_LEF_PRIMARY}};
  const std::vector<CollisionT> fwd_collisions_expected{{},
                                                        {},
                                                        {3, LEF_LEF_PRIMARY},
                                                        {4, LEF_LEF_SECONDARY},
                                                        {5, LEF_LEF_PRIMARY},
                                                        {6, LEF_LEF_PRIMARY},
                                                        {3, CHROM_BOUNDARY}};
  std::vector<CollisionT> rev_collisions(nlefs);
  std::vector<CollisionT> fwd_collisions(nlefs);

  const std::vector<bp_t> rev_moves_expected{75, 75, 75, 13, 61, 25, 60};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{75, 75, 12, 53, 24, 59, 24};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                     // NOLINT
  assert(barriers.size() == nbarriers);             // NOLINT
  assert(rev_collisions_expected.size() == nlefs);  // NOLINT
  assert(fwd_collisions_expected.size() == nlefs);  // NOLINT
  assert(rev_collisions.size() == nlefs);           // NOLINT
  assert(fwd_collisions.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_moves(chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
                                absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      barriers, barrier_mask, absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 005 - Multiple LEFs located at the same site", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 25;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 25;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_bar_major_collision_pblock = 1;        // NOLINT
  c.lef_bar_minor_collision_pblock = 0;        // NOLINT
  constexpr usize nlefs = 6;
  constexpr usize nbarriers = 1;
  const Chromosome chrom{0, "chr1", 0, 150, 150};  // NOLINT
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = random::PRNG(9509334554644345044ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      construct_lef(30, 50, 0),    // NOLINT
      construct_lef(60, 80, 1),    // NOLINT
      construct_lef(60, 80, 2),    // NOLINT
      construct_lef(65, 125, 3),   // NOLINT
      construct_lef(140, 140, 4),  // NOLINT
      construct_lef(140, 140, 5)   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '-'}};

  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<usize> rev_ranks{0, 1, 2, 3, 4, 5};  // NOLINT
  const std::vector<usize> fwd_ranks{0, 1, 2, 3, 4, 5};  // NOLINT

  std::vector<bp_t> rev_moves{25, 25, 25, 25, 25, 25};  // NOLINT
  std::vector<bp_t> fwd_moves{25, 25, 25, 24,  8,  9};  // NOLINT
  // clang-format off
  const std::vector<CollisionT> rev_collisions_expected{{},
                                                        {0, LEF_LEF_PRIMARY},
                                                        {1, LEF_LEF_SECONDARY},
                                                        {2, LEF_LEF_SECONDARY},
                                                        {3, LEF_LEF_PRIMARY},
                                                        {4, LEF_LEF_SECONDARY}};
  const std::vector<CollisionT> fwd_collisions_expected{{1, LEF_LEF_PRIMARY},
                                                        {2, LEF_LEF_SECONDARY},
                                                        {0, LEF_BAR},
                                                        {4, LEF_LEF_PRIMARY},
                                                        {5, LEF_LEF_SECONDARY},
                                                        {3, CHROM_BOUNDARY}};
  std::vector<CollisionT> rev_collisions(nlefs);
  std::vector<CollisionT> fwd_collisions(nlefs);

  const std::vector<bp_t> rev_moves_expected{25, 5, 4, 8, 8, 7};   // NOLINT
  const std::vector<bp_t> fwd_moves_expected{4, 18, 19, 6, 8, 9};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                     // NOLINT
  assert(barriers.size() == nbarriers);             // NOLINT
  assert(rev_collisions_expected.size() == nlefs);  // NOLINT
  assert(fwd_collisions_expected.size() == nlefs);  // NOLINT
  assert(rev_collisions.size() == nlefs);           // NOLINT
  assert(fwd_collisions.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_moves(chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
                                absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      barriers, barrier_mask, absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 006 - Few inactive LEFs", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 25;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 25;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_bar_major_collision_pblock = 1;        // NOLINT
  c.lef_bar_minor_collision_pblock = 0;        // NOLINT
  constexpr usize nlefs = 6;
  constexpr usize nbarriers = 1;
  const Chromosome chrom{0, "chr1", 0, 150, 150};  // NOLINT
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = random::PRNG(8213068426516999476ULL);

  // clang-format off
  std::vector<Lef> lefs{
      construct_lef(30, 50, 0),    // NOLINT
      construct_lef(60, 80, 1),    // NOLINT
      construct_lef(60, 80, 2),    // NOLINT
      construct_lef(65, 125, 3),   // NOLINT
      construct_lef(140, 140, 4),  // NOLINT
      construct_lef(140, 140, 5)   // NOLINT
  };

  lefs[2].release(); // NOLINT
  lefs[5].release(); // NOLINT

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '-'}};

  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<usize> rev_ranks{0, 1, 3, 4, 2, 5};  // NOLINT
  const std::vector<usize> fwd_ranks{0, 1, 3, 4, 2, 5};  // NOLINT

  std::vector<bp_t> rev_moves{25, 25,  0, 25, 25,  0};  // NOLINT
  std::vector<bp_t> fwd_moves{25, 25,  0, 24,  9,  0};  // NOLINT
  // clang-format off
  const std::vector<CollisionT> rev_collisions_expected{{},
                                                        {0, LEF_LEF_PRIMARY},
                                                        {},
                                                        {1, LEF_LEF_SECONDARY},
                                                        {3, LEF_LEF_PRIMARY},
                                                        {}};
  const std::vector<CollisionT> fwd_collisions_expected{{1, LEF_LEF_PRIMARY},
                                                        {0, LEF_BAR},
                                                        {},
                                                        {4, LEF_LEF_PRIMARY},
                                                        {3, CHROM_BOUNDARY},
                                                        {}};
  std::vector<CollisionT> rev_collisions(nlefs);
  std::vector<CollisionT> fwd_collisions(nlefs);

  const std::vector<bp_t> rev_moves_expected{25, 5, 0, 9, 8, 0};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{4, 19, 0, 6, 9, 0};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                     // NOLINT
  assert(barriers.size() == nbarriers);             // NOLINT
  assert(rev_collisions_expected.size() == nlefs);  // NOLINT
  assert(fwd_collisions_expected.size() == nlefs);  // NOLINT
  assert(rev_collisions.size() == nlefs);           // NOLINT
  assert(fwd_collisions.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_moves(chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
                                absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      barriers, barrier_mask, absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 007 - LEF-LEF collision overrides LEF-BAR collision 1",
          "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 20;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 20;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_bar_major_collision_pblock = 1;        // NOLINT
  c.lef_bar_minor_collision_pblock = 0;        // NOLINT
  constexpr usize nlefs = 2;
  constexpr usize nbarriers = 1;
  const Chromosome chrom{0, "chr1", 0, 200, 200};  // NOLINT
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = random::PRNG(4003028292286930594ULL);

  // clang-format off
  std::vector<Lef> lefs{
      construct_lef(50, 95, 0),    // NOLINT
      construct_lef(110, 150, 1)   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '+'}};

  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<usize> rev_ranks{0, 1};  // NOLINT
  const std::vector<usize> fwd_ranks{0, 1};  // NOLINT

  std::vector<bp_t> rev_moves{20, 20};  // NOLINT
  std::vector<bp_t> fwd_moves{20, 20};  // NOLINT
  // clang-format on

  const std::vector<CollisionT> rev_collisions_expected{{}, {0, LEF_LEF_PRIMARY}};
  const std::vector<CollisionT> fwd_collisions_expected{{1, LEF_LEF_PRIMARY}, {}};
  std::vector<CollisionT> rev_collisions(nlefs);
  std::vector<CollisionT> fwd_collisions(nlefs);

  const std::vector<bp_t> rev_moves_expected{20, 7};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{7, 20};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                     // NOLINT
  assert(barriers.size() == nbarriers);             // NOLINT
  assert(rev_collisions_expected.size() == nlefs);  // NOLINT
  assert(fwd_collisions_expected.size() == nlefs);  // NOLINT
  assert(rev_collisions.size() == nlefs);           // NOLINT
  assert(fwd_collisions.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  modle::Simulation{c, false}.test_process_collisions(
      chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      barriers, barrier_mask, absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 008 - LEF-LEF collision overrides LEF-BAR collision 2",
          "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 20;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 20;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_bar_major_collision_pblock = 1;        // NOLINT
  c.lef_bar_minor_collision_pblock = 0;        // NOLINT
  constexpr usize nlefs = 2;
  constexpr usize nbarriers = 1;
  const Chromosome chrom{0, "chr1", 0, 200, 200};  // NOLINT
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = random::PRNG(12931236264094635263ULL);

  // clang-format off
  std::vector<Lef> lefs{
      construct_lef(50, 90, 0),    // NOLINT
      construct_lef(105, 150, 1)   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '-'}};

  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<usize> rev_ranks{0, 1};  // NOLINT
  const std::vector<usize> fwd_ranks{0, 1};  // NOLINT

  std::vector<bp_t> rev_moves{20, 20};  // NOLINT
  std::vector<bp_t> fwd_moves{20, 20};  // NOLINT
  // clang-format on

  const std::vector<CollisionT> rev_collisions_expected{{}, {0, LEF_LEF_PRIMARY}};
  const std::vector<CollisionT> fwd_collisions_expected{{1, LEF_LEF_PRIMARY}, {}};
  std::vector<CollisionT> rev_collisions(nlefs);
  std::vector<CollisionT> fwd_collisions(nlefs);

  const std::vector<bp_t> rev_moves_expected{20, 7};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{7, 20};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                     // NOLINT
  assert(barriers.size() == nbarriers);             // NOLINT
  assert(rev_collisions_expected.size() == nlefs);  // NOLINT
  assert(fwd_collisions_expected.size() == nlefs);  // NOLINT
  assert(rev_collisions.size() == nlefs);           // NOLINT
  assert(fwd_collisions.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  modle::Simulation{c, false}.test_process_collisions(
      chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      barriers, barrier_mask, absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 009 - Ensure stacked LEFs do not interfere with surrounding extr. barriers",
          "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 10;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 10;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_bar_major_collision_pblock = 1;        // NOLINT
  c.lef_bar_minor_collision_pblock = 0;        // NOLINT
  constexpr usize nlefs = 5;
  constexpr usize nbarriers = 2;
  const Chromosome chrom{0, "chr1", 0, 200, 200};  // NOLINT
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = random::PRNG(8668729858575330735ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      construct_lef(95, 100, 0),   // NOLINT
      construct_lef(101, 103, 1),  // NOLINT
      construct_lef(102, 110, 2),  // NOLINT
      construct_lef(104, 111, 3),  // NOLINT
      construct_lef(105, 112, 4)   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{105, 1.0, 0.0, '-'}};

  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<usize> rev_ranks{0, 1, 2, 3, 4};  // NOLINT
  const std::vector<usize> fwd_ranks{0, 1, 2, 3, 4};  // NOLINT

  std::vector<bp_t> rev_moves{10, 10, 10, 10, 10};  // NOLINT
  std::vector<bp_t> fwd_moves{10, 10, 10, 10, 10};  // NOLINT
  // clang-format off
  const std::vector<CollisionT> rev_collisions_expected{{},
                                                        {0, LEF_BAR},
                                                        {1, LEF_LEF_SECONDARY},
                                                        {1, LEF_LEF_PRIMARY},
                                                        {3, LEF_LEF_SECONDARY}};
  const std::vector<CollisionT> fwd_collisions_expected{{1, LEF_LEF_PRIMARY},
                                                        {3, LEF_LEF_PRIMARY},
                                                        {},
                                                        {},
                                                        {}};
  std::vector<CollisionT> rev_collisions(nlefs);
  std::vector<CollisionT> fwd_collisions(nlefs);

  const std::vector<bp_t> rev_moves_expected{10, 0, 0, 0, 0};    // NOLINT
  const std::vector<bp_t> fwd_moves_expected{0, 0, 10, 10, 10};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                     // NOLINT
  assert(barriers.size() == nbarriers);             // NOLINT
  assert(rev_collisions_expected.size() == nlefs);  // NOLINT
  assert(fwd_collisions_expected.size() == nlefs);  // NOLINT
  assert(rev_collisions.size() == nlefs);           // NOLINT
  assert(fwd_collisions.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_moves(chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
                                absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      barriers, barrier_mask, absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Simulation 010 - Ensure stacked LEFs do not interfere with surrounding extr. barriers",
          "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 10;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 10;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_bar_major_collision_pblock = 1;        // NOLINT
  c.lef_bar_minor_collision_pblock = 0;        // NOLINT
  constexpr usize nlefs = 6;
  constexpr usize nbarriers = 2;
  const Chromosome chrom{0, "chr1", 0, 200, 200};  // NOLINT
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = random::PRNG(4320101097510038983ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      construct_lef(95, 100, 0),   // NOLINT
      construct_lef(101, 103, 1),  // NOLINT
      construct_lef(102, 110, 2),  // NOLINT
      construct_lef(104, 111, 3),  // NOLINT
      construct_lef(105, 112, 4),  // NOLINT
      construct_lef(102, 102, 5)   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{105, 1.0, 0.0, '-'}};

  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<usize> rev_ranks{0, 1, 2, 5, 3, 4};  // NOLINT
  const std::vector<usize> fwd_ranks{0, 5, 1, 2, 3, 4};  // NOLINT

  std::vector<bp_t> rev_moves{10, 10, 10, 10, 10, 10};  // NOLINT
  std::vector<bp_t> fwd_moves{10, 10, 10, 10, 10, 10};  // NOLINT
  // clang-format off
  const std::vector<CollisionT> rev_collisions_expected{{},
                                                        {0, LEF_BAR},
                                                        {1, LEF_LEF_SECONDARY},
                                                        {1, LEF_LEF_PRIMARY},
                                                        {3, LEF_LEF_SECONDARY},
                                                        {2, LEF_LEF_SECONDARY}};
  const std::vector<CollisionT> fwd_collisions_expected{{1, LEF_LEF_PRIMARY},
                                                        {3, LEF_LEF_PRIMARY},
                                                        {},
                                                        {},
                                                        {},
                                                        {1, LEF_LEF_SECONDARY}};
  std::vector<CollisionT> rev_collisions(nlefs);
  std::vector<CollisionT> fwd_collisions(nlefs);

  const std::vector<bp_t> rev_moves_expected{10, 0, 0, 0, 0, 0};    // NOLINT
  const std::vector<bp_t> fwd_moves_expected{0, 0, 10, 10, 10, 0};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                     // NOLINT
  assert(barriers.size() == nbarriers);             // NOLINT
  assert(rev_collisions_expected.size() == nlefs);  // NOLINT
  assert(fwd_collisions_expected.size() == nlefs);  // NOLINT
  assert(rev_collisions.size() == nlefs);           // NOLINT
  assert(fwd_collisions.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation::test_adjust_moves(chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
                                absl::MakeSpan(fwd_moves));

  modle::Simulation{c, false}.test_process_collisions(
      chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      barriers, barrier_mask, absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

}  // namespace modle::test::libmodle
