// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/types/span.h>

#include <algorithm>
#include <cassert>
#include <catch2/catch_test_macros.hpp>
#include <iterator>
#include <memory>
#include <numeric>
#include <string_view>
#include <vector>

#include "./common.hpp"
#include "modle/common/common.hpp"
#include "modle/common/random.hpp"
#include "modle/common/simulation_config.hpp"
#include "modle/extrusion_barriers.hpp"
#include "modle/extrusion_factors.hpp"
#include "modle/genome.hpp"
#include "modle/simulation.hpp"

namespace modle::libmodle::test {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Bind LEFs 001", "[bind-lefs][simulation][short]") {
  const auto chrom = init_interval("chr1", 1000);
  const std::size_t nlefs = 10;
  std::array<Lef, nlefs> lefs;
  std::array<std::size_t, nlefs> rank1;
  std::array<std::size_t, nlefs> rank2;
  std::iota(rank1.begin(), rank1.end(), 0);
  std::copy(rank1.begin(), rank1.end(), rank2.begin());
  std::array<std::size_t, nlefs> mask{};
  std::fill(mask.begin(), mask.end(), 1);
  auto rand_eng = DEFAULT_PRNG;

  for (std::size_t i = 0; i < nlefs; ++i) {  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
    mask[i] = random::bernoulli_trial{0.50}(rand_eng);
  }

  const auto c = Config{};
  Simulation{c, false}.test_bind_lefs(chrom, lefs, rank1, rank2, mask, rand_eng, 0);

  CHECK(lefs.size() == nlefs);
  CHECK(rank1.size() == nlefs);
  CHECK(rank2.size() == nlefs);
  CHECK(mask.size() == nlefs);
  check_that_lefs_are_sorted_by_idx(lefs, rank1, rank2);

  for (std::size_t i = 0; i < nlefs; ++i) {
    const auto& lef = lefs[i];
    if (mask[i]) {  // NOLINT(readability-implicit-bool-conversion)
      CHECK(lef.is_bound());
      CHECK(lef.rev_unit.pos() == lef.fwd_unit.pos());
      CHECK(lef.rev_unit.pos() >= chrom.start());
      CHECK(lef.rev_unit.pos() < chrom.end());
    } else {
      CHECK(!lef.is_bound());
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Bind LEFs 002 - No LEFs to bind", "[bind-lefs][simulation][short]") {
  const auto chrom = init_interval("chr1", 1000);
  std::vector<Lef> lefs;
  std::vector<std::size_t> rank1, rank2;
  boost::dynamic_bitset<> mask1;
  std::vector<int> mask2;
  auto rand_eng = DEFAULT_PRNG;
  const auto c = Config{};

  Simulation{c, false}.test_bind_lefs(chrom, lefs, rank1, rank2, mask1, rand_eng, 1);

  CHECK(lefs.empty());
  CHECK(rank1.empty());
  CHECK(rank2.empty());
  CHECK(mask1.empty());

  Simulation{c, false}.test_bind_lefs(chrom, lefs, rank1, rank2, mask2, rand_eng, 2);

  CHECK(lefs.empty());
  CHECK(rank1.empty());
  CHECK(rank2.empty());
  CHECK(mask2.empty());

  const std::size_t nlefs = 10;
  lefs = std::vector<Lef>(nlefs, Lef{});
  rank1.resize(nlefs);
  std::iota(rank1.begin(), rank1.end(), 0);
  rank2 = rank1;
  mask1.resize(nlefs, false);

  Simulation{c, false}.test_bind_lefs(chrom, lefs, rank1, rank2, mask1, rand_eng, 0);

  CHECK(lefs.size() == nlefs);
  CHECK(rank1.size() == nlefs);
  CHECK(rank2.size() == nlefs);
  CHECK(mask1.size() == nlefs);
  check_that_lefs_are_sorted_by_idx(lefs, rank1, rank2);
  CHECK(mask1.count() == 0);

  CHECK(std::all_of(lefs.begin(), lefs.end(), [](const auto& lef) { return !lef.is_bound(); }));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Bind LEFs 003 - Empty mask (i.e. bind all LEFs)", "[bind-lefs][simulation][short]") {
  const auto chrom = init_interval("chr1", 1000);
  const std::size_t nlefs = 10;
  std::array<Lef, nlefs> lefs{};
  std::array<std::size_t, nlefs> rank1{};
  std::array<std::size_t, nlefs> rank2{};
  std::iota(rank1.begin(), rank1.end(), 0);
  std::copy(rank1.begin(), rank1.end(), rank2.begin());
  boost::dynamic_bitset<> mask;
  auto rand_eng = DEFAULT_PRNG;
  const auto c = Config{};

  Simulation{c, false}.test_bind_lefs(chrom, lefs, rank1, rank2, mask, rand_eng, 0);

  CHECK(mask.empty());

  CHECK(std::all_of(lefs.begin(), lefs.end(), [](const auto& lef) { return lef.is_bound(); }));
  CHECK(std::all_of(lefs.begin(), lefs.end(),
                    [](const auto& lef) { return lef.rev_unit.pos() == lef.fwd_unit.pos(); }));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Adjust LEF moves 001", "[adjust-lef-moves][simulation][short]") {
  const auto chrom = init_interval("chr1", 101);
  const std::size_t nlefs = 3;
  // clang-format off
  const std::array<Lef, nlefs> lefs{construct_lef(5, 25, 1),
                                    construct_lef(10, 20, 2),
                                    construct_lef(90, 90, 3)};

  // clang-format on
  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 2};
  const std::array<std::size_t, nlefs> fwd_ranks{1, 0, 2};

  std::array<bp_t, nlefs> rev_moves{5, 10, 15};
  std::array<bp_t, nlefs> fwd_moves{10, 20, 10};

  const std::array<bp_t, nlefs> rev_moves_adjusted{5, 10, 15};
  const std::array<bp_t, nlefs> fwd_moves_adjusted{16, 20, 10};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  auto c = Config{};
  c.rev_extrusion_speed_std = 1;
  c.fwd_extrusion_speed_std = 1;
  Simulation::test_adjust_and_clamp_moves(chrom, lefs, rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  for (std::size_t i = 0; i < lefs.size(); ++i) {
    CHECK(rev_moves[i] == rev_moves_adjusted[i]);
    CHECK(fwd_moves[i] == fwd_moves_adjusted[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Adjust LEF moves 002", "[adjust-lef-moves][simulation][short]") {
  const auto chrom = init_interval("chr1", 400, 10);
  const std::size_t nlefs = 6;
  // clang-format off
  const std::array<Lef, nlefs> lefs{construct_lef(20, 50, 0),
                                    construct_lef(60, 60, 1),
                                    construct_lef(200, 310, 2),
                                    construct_lef(220, 300, 3),
                                    construct_lef(240, 250, 4),
                                    construct_lef(125, 305, 5)};

  // clang-format on
  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 5, 2, 3, 4};
  const std::array<std::size_t, nlefs> fwd_ranks{0, 1, 4, 3, 5, 2};

  std::array<bp_t, nlefs> rev_moves{10, 10, 5, 25, 50, 10};
  std::array<bp_t, nlefs> fwd_moves{25, 10, 5, 20, 20, 0};

  const std::array<bp_t, nlefs> rev_moves_adjusted{10, 10, 12, 31, 50, 10};
  const std::array<bp_t, nlefs> fwd_moves_adjusted{25, 16, 12, 20, 20, 16};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  auto c = Config{};
  c.rev_extrusion_speed_std = 1;
  c.fwd_extrusion_speed_std = 1;
  Simulation::test_adjust_and_clamp_moves(chrom, lefs, rev_ranks, fwd_ranks,
                                          absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));

  for (std::size_t i = 0; i < nlefs; ++i) {
    CHECK(rev_moves[i] == rev_moves_adjusted[i]);
    CHECK(fwd_moves[i] == fwd_moves_adjusted[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Generate LEF moves 001", "[generate-lef-moves][simulation][medium]") {
  const auto chrom = init_interval("chr1", 2000, 1000);
  const std::size_t nlefs = 100;
  const std::size_t iters = 1000;
  auto rand_eng = DEFAULT_PRNG;
  std::array<Lef, nlefs> lefs{};
  std::array<std::size_t, nlefs> rev_ranks{};
  std::array<std::size_t, nlefs> fwd_ranks{};
  std::array<bp_t, nlefs> rev_moves{};
  std::array<bp_t, nlefs> fwd_moves{};
  auto c = Config{};
  c.bin_size = 500;
  c.rev_extrusion_speed = c.bin_size;
  c.fwd_extrusion_speed = c.bin_size;
  c.rev_extrusion_speed_std = static_cast<double>(c.rev_extrusion_speed) * 0.2;
  c.fwd_extrusion_speed_std = static_cast<double>(c.fwd_extrusion_speed) * 0.2;

  Simulation simulation{c, false};

  for (std::size_t i = 0; i < iters; ++i) {
    std::generate(lefs.begin(), lefs.end(), [&]() {
      const auto pos =
          random::uniform_int_distribution<bp_t>{chrom.start(), chrom.end() - 1}(rand_eng);
      return construct_lef(pos, pos, i);
    });
    Simulation::test_rank_lefs(lefs, absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks), false,
                               true);

    simulation.test_generate_moves(chrom, absl::MakeConstSpan(lefs), absl::MakeConstSpan(rev_ranks),
                                   absl::MakeConstSpan(fwd_ranks), absl::MakeSpan(rev_moves),
                                   absl::MakeSpan(fwd_moves), rand_eng);

    std::for_each(rev_moves.begin(), rev_moves.end(), [&, j = 0UL](const auto n) mutable {
      CHECK(lefs[j].rev_unit.pos() >= chrom.start() + n);
      CHECK(lefs[j++].rev_unit.pos() < chrom.end());
    });

    std::for_each(fwd_moves.begin(), fwd_moves.end(), [&, j = 0UL](const auto n) mutable {
      CHECK(lefs[j].fwd_unit.pos() + n < chrom.end());
      CHECK(lefs[j++].fwd_unit.pos() >= chrom.start());
    });
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Detect LEF-LEF collisions 001", "[lef-lef-collisions][simulation][short]") {
  const auto c = init_config(3, 2);
  const auto chrom = init_interval("chr1", 30);
  [[maybe_unused]] const std::size_t nlefs = 4;
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  /*       __0__      __1__        2        __3__
  *        |   |      |   |        |        |   |
  *      <-0- -2->  <-4- -8->  <--14-->  <-18- -23->
  */
  std::array<Lef, nlefs> lefs{
    construct_lef(0, 2, 0),
    construct_lef(4, 8, 1),
    construct_lef(14, 14, 2),
    construct_lef(18, 23, 3)
  };
  // clang-format on
  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 2, 3};
  const std::array<std::size_t, nlefs> fwd_ranks{0, 1, 2, 3};

  const ExtrusionBarriers barriers{};

  std::array<CollisionT, nlefs> fwd_collisions;
  std::array<CollisionT, nlefs> rev_collisions;

  std::array<bp_t, nlefs> rev_moves{0, 3, 3, 3};
  std::array<bp_t, nlefs> fwd_moves{2, 2, 2, 2};
  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{5, CHROM_BOUNDARY},
                                                              CollisionT{0, LEF_LEF_PRIMARY},
                                                              CollisionT{},
                                                              CollisionT{2, LEF_LEF_PRIMARY}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{},
                                                              CollisionT{3, LEF_LEF_PRIMARY},
                                                              CollisionT{}};
  // clang-format on

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  modle::Simulation::test_detect_units_at_interval_boundaries(
      chrom, lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions));

  modle::Simulation{c, false}.test_detect_primary_lef_lef_collisions(
      lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions), rand_eng);

  check_collisions(lefs, rev_collisions, rev_collisions_expected, fwd_collisions,
                   fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Process LEF-LEF collisions 001", "[lef-lef-collisions][simulation][short]") {
  const auto c = init_config(3, 2);
  const auto chrom = init_interval("chr1", 30);
  [[maybe_unused]] const std::size_t nlefs = 4;
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
/*       __0__      __1__        2        __3__
*        |   |      |   |        |        |   |
*      <-0- -2->  <-4- -8->  <--14-->  <-18- -23->
*/
  std::array<Lef, nlefs> lefs{
    construct_lef(0, 2, 0),
    construct_lef(4, 8, 1),
    construct_lef(14, 14, 2),
    construct_lef(18, 23, 3)
  };
  // clang-format on
  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 2, 3};
  const std::array<std::size_t, nlefs> fwd_ranks{0, 1, 2, 3};

  const ExtrusionBarriers barriers{};

  std::array<CollisionT, nlefs> fwd_collisions;
  std::array<CollisionT, nlefs> rev_collisions;

  std::array<bp_t, nlefs> rev_moves{0, 3, 3, 3};
  std::array<bp_t, nlefs> fwd_moves{2, 2, 2, 2};
  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{5, CHROM_BOUNDARY},
                                                              CollisionT{0, LEF_LEF_PRIMARY},
                                                              CollisionT{},
                                                              CollisionT{2, LEF_LEF_PRIMARY}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{},
                                                              CollisionT{3, LEF_LEF_PRIMARY},
                                                              CollisionT{}};
  // clang-format on

  const std::array<bp_t, nlefs> rev_moves_expected{0, 1, 3, 2};
  const std::array<bp_t, nlefs> fwd_moves_expected{0, 2, 1, 2};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  modle::Simulation::test_detect_units_at_interval_boundaries(
      chrom, lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions));

  modle::Simulation{c, false}.test_process_lef_lef_collisions(
      chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Detect LEF-LEF collisions 002", "[lef-lef-collisions][simulation][short]") {
  const auto c = init_config(3, 2);
  const auto chrom = init_interval("chr1", 16);
  [[maybe_unused]] const std::size_t nlefs = 4;
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
/*
 *                 _______1_______                _______3_______
 *        _______0_|_____        |      ________2_|______       |
 *        |        |    |        |      |         |     |       |
 *      <-0-     <-4-  -5->     -6->  <-9-     <-11-   -14->   -15->
 */

std::vector<Lef> lefs{
    construct_lef(0, 5, 0),
    construct_lef(4, 6, 1),
    construct_lef(9, 14, 2),
    construct_lef(11, 15, 3)
};
  // clang-format on
  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 2, 3};
  const std::array<std::size_t, nlefs> fwd_ranks{0, 1, 2, 3};

  const ExtrusionBarriers barriers{};

  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  std::array<bp_t, nlefs> rev_moves{0, 3, 3, 4};
  std::array<bp_t, nlefs> fwd_moves{3, 2, 1, 0};
  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{5, CHROM_BOUNDARY},
                                                              CollisionT{},
                                                              CollisionT{1, LEF_LEF_PRIMARY},
                                                              CollisionT{2, LEF_LEF_SECONDARY}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{1, LEF_LEF_SECONDARY},
                                                              CollisionT{2, LEF_LEF_PRIMARY},
                                                              CollisionT{3, CHROM_BOUNDARY},
                                                              CollisionT{3, CHROM_BOUNDARY}};
  // clang-format on

  const std::array<bp_t, nlefs> rev_moves_expected{0, 3, 2, 3};
  const std::array<bp_t, nlefs> fwd_moves_expected{0, 0, 1, 0};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  modle::Simulation::test_detect_units_at_interval_boundaries(
      chrom, lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions));

  modle::Simulation{c, false}.test_process_lef_lef_collisions(
      chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Detect LEF-LEF collisions 003", "[lef-lef-collisions][simulation][short]") {
  const auto c = init_config(3, 2);
  const auto chrom = init_interval("chr1", 201, 100);
  const std::size_t nlefs = 3;
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
/*
 *                ___________________0___________________
 *                |       ___________1___________       |
 *                |       |       ___2___       |       |
 *                |       |       |     |       |       |
 *             <-120-  <-130-  <-140- -141->  -160->  -180->
 */

std::vector<Lef> lefs{
    construct_lef(120, 180, 0),
    construct_lef(130, 160, 1),
    construct_lef(140, 141, 2)
};
  // clang-format on
  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 2};
  const std::array<std::size_t, nlefs> fwd_ranks{2, 1, 0};

  const ExtrusionBarriers barriers{};

  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  std::array<bp_t, nlefs> rev_moves{20, 30, 40};
  std::array<bp_t, nlefs> fwd_moves{20, 40, 59};
  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{5, CHROM_BOUNDARY},
                                                              CollisionT{0, LEF_LEF_SECONDARY},
                                                              CollisionT{1, LEF_LEF_SECONDARY}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{3, CHROM_BOUNDARY},
                                                              CollisionT{0, LEF_LEF_SECONDARY},
                                                              CollisionT{1, LEF_LEF_SECONDARY}};
  // clang-format on

  const std::array<bp_t, nlefs> rev_moves_expected{20, 29, 38};
  const std::array<bp_t, nlefs> fwd_moves_expected{20, 39, 57};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  modle::Simulation::test_detect_units_at_interval_boundaries(
      chrom, lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions));

  modle::Simulation{c, false}.test_process_lef_lef_collisions(
      chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions),
      rand_eng);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Detect LEF-BAR collisions 001 - wo soft collisions fwd CTCFs",
          "[lef-bar-collisions][simulation][short]") {
  const auto c = init_config(2, 2);
  constexpr std::size_t nlefs = 3;
  constexpr std::size_t nbarriers = 3;
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(0, 1, 0),
      construct_lef(3, 4, 1),
      construct_lef(5, 5, 2)
  };
  const auto barriers = construct_barriers(ExtrusionBarrier{2, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{4, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{8, 1.0, 0.0, '+'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 2};
  const std::array<std::size_t, nlefs> fwd_ranks{0, 1, 2};

  std::array<bp_t, nlefs> rev_moves{0, 2, 2};
  std::array<bp_t, nlefs> fwd_moves{2, 2, 2};

  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{0, LEF_BAR},
                                                              CollisionT{1, LEF_BAR}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{}};
  // clang-format on

  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{0, 0, 0};
  const std::array<bp_t, nlefs> fwd_moves_expected{2, 2, 2};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation{c, false}.test_detect_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, barriers, absl::MakeSpan(rev_collisions),
      absl::MakeSpan(fwd_collisions), rand_eng);
  Simulation::test_correct_moves_for_lef_bar_collisions(lefs, barriers, absl::MakeSpan(rev_moves),
                                                        absl::MakeSpan(fwd_moves), rev_collisions,
                                                        fwd_collisions);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Detect LEF-BAR collisions 002 - wo soft-collisions rev CTCFs",
          "[lef-bar-collisions][simulation][short]") {
  const auto c = init_config(2, 2);
  constexpr std::size_t nlefs = 3;
  constexpr std::size_t nbarriers = 3;
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(0, 1, 0),
      construct_lef(3, 4, 1),
      construct_lef(5, 5, 2)
  };
  const auto barriers = construct_barriers(ExtrusionBarrier{2, 1.0, 0.0, '-'},
                                           ExtrusionBarrier{4, 1.0, 0.0, '-'},
                                           ExtrusionBarrier{8, 1.0, 0.0, '-'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 2};
  const std::array<std::size_t, nlefs> fwd_ranks{0, 1, 2};

  std::array<bp_t, nlefs> rev_moves{0, 2, 2};
  std::array<bp_t, nlefs> fwd_moves{2, 2, 2};

  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{0, LEF_BAR},
                                                              CollisionT{},
                                                              CollisionT{}};
  // clang-format on

  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{0, 2, 2};
  const std::array<bp_t, nlefs> fwd_moves_expected{0, 2, 2};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation{c, false}.test_detect_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions), rand_eng);
  Simulation::test_correct_moves_for_lef_bar_collisions(lefs, barriers, absl::MakeSpan(rev_moves),
                                                        absl::MakeSpan(fwd_moves), rev_collisions,
                                                        fwd_collisions);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Detect LEF-BAR collisions 003 - w soft-collisions fwd CTCFs",
          "[lef-bar-collisions][simulation][short]") {
  auto c = init_config(2, 2);
  c.lef_bar_major_collision_pblock = 1;
  c.lef_bar_minor_collision_pblock = 1;
  constexpr std::size_t nlefs = 3;
  constexpr std::size_t nbarriers = 3;
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(0, 1, 0),
      construct_lef(3, 4, 1),
      construct_lef(5, 5, 2)
  };
  const auto barriers = construct_barriers(ExtrusionBarrier{2, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{4, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{8, 1.0, 0.0, '+'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 2};
  const std::array<std::size_t, nlefs> fwd_ranks{0, 1, 2};

  std::array<bp_t, nlefs> rev_moves{0, 2, 2};
  std::array<bp_t, nlefs> fwd_moves{2, 2, 2};

  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{0, LEF_BAR},
                                                              CollisionT{1, LEF_BAR}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{0, LEF_BAR},
                                                              CollisionT{},
                                                              CollisionT{}};
  // clang-format on

  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{0, 0, 0};
  const std::array<bp_t, nlefs> fwd_moves_expected{0, 2, 2};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation{c, false}.test_detect_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions), rand_eng);
  Simulation::test_correct_moves_for_lef_bar_collisions(lefs, barriers, absl::MakeSpan(rev_moves),
                                                        absl::MakeSpan(fwd_moves), rev_collisions,
                                                        fwd_collisions);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Detect LEF-BAR collisions 004 - wo soft-collisions mixed CTCFs",
          "[lef-bar-collisions][simulation][short]") {
  const auto c = init_config(5, 5);
  const std::size_t nlefs = 5;
  const std::size_t nbarriers = 4;
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(10, 20, 0),
      construct_lef(26, 26, 1),
      construct_lef(30, 35, 2),
      construct_lef(42, 43, 3),
      construct_lef(44, 60, 4)
  };

  const auto barriers = construct_barriers(ExtrusionBarrier{25, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{27, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{40, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{46, 1.0, 0.0, '-'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 2, 3, 4};
  const std::array<std::size_t, nlefs> fwd_ranks{0, 1, 2, 3, 4};

  std::array<bp_t, nlefs> rev_moves{5, 5, 5, 5, 5};
  std::array<bp_t, nlefs> fwd_moves{5, 5, 5, 5, 5};

  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{0, LEF_BAR},
                                                              CollisionT{1, LEF_BAR},
                                                              CollisionT{2, LEF_BAR},
                                                              CollisionT{}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{3, LEF_BAR},
                                                              CollisionT{}};
  // clang-format on

  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{5, 0, 2, 1, 5};
  const std::array<bp_t, nlefs> fwd_moves_expected{5, 5, 5, 2, 5};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation{c, false}.test_detect_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions), rand_eng);
  Simulation::test_correct_moves_for_lef_bar_collisions(lefs, barriers, absl::MakeSpan(rev_moves),
                                                        absl::MakeSpan(fwd_moves), rev_collisions,
                                                        fwd_collisions);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Detect LEF-BAR collisions 005 - wo soft-collisions mixed CTCFs, different extr. speeds",
          "[lef-bar-collisions][simulation][short]") {
  const auto c = init_config(2, 5);
  const std::size_t nlefs = 5;
  const std::size_t nbarriers = 4;
  auto rand_eng = DEFAULT_PRNG;

  // clang-format off
  std::array<Lef, nlefs> lefs{
      construct_lef(10, 20, 0),
      construct_lef(26, 26, 1),
      construct_lef(30, 35, 2),
      construct_lef(42, 43, 3),
      construct_lef(44, 60, 4)
  };

  const auto barriers = construct_barriers(ExtrusionBarrier{25, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{27, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{40, 1.0, 0.0, '+'},
                                           ExtrusionBarrier{46, 1.0, 0.0, '-'});
  // clang-format on
  REQUIRE(barriers.size() == nbarriers);

  const std::array<std::size_t, nlefs> rev_ranks{0, 1, 2, 3, 4};
  const std::array<std::size_t, nlefs> fwd_ranks{0, 1, 2, 3, 4};

  std::array<bp_t, nlefs> rev_moves{2, 2, 2, 2, 2};
  std::array<bp_t, nlefs> fwd_moves{5, 5, 5, 5, 5};

  // clang-format off
  const std::array<CollisionT, nlefs> rev_collisions_expected{CollisionT{},
                                                              CollisionT{0, LEF_BAR},
                                                              CollisionT{},
                                                              CollisionT{2, LEF_BAR},
                                                              CollisionT{}};
  const std::array<CollisionT, nlefs> fwd_collisions_expected{CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{},
                                                              CollisionT{3, LEF_BAR},
                                                              CollisionT{}};
  // clang-format on

  std::array<CollisionT, nlefs> rev_collisions;
  std::array<CollisionT, nlefs> fwd_collisions;

  const std::array<bp_t, nlefs> rev_moves_expected{2, 0, 2, 1, 2};
  const std::array<bp_t, nlefs> fwd_moves_expected{5, 5, 5, 2, 5};

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  Simulation{c, false}.test_detect_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      absl::MakeSpan(rev_collisions), absl::MakeSpan(fwd_collisions), rand_eng);
  Simulation::test_correct_moves_for_lef_bar_collisions(lefs, barriers, absl::MakeSpan(rev_moves),
                                                        absl::MakeSpan(fwd_moves), rev_collisions,
                                                        fwd_collisions);

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collisions, rev_collisions_expected, fwd_collisions,
                          fwd_collisions_expected);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("LEFs ranking 001 - Rev extr. unit tied", "[simulation][short]") {
  modle::Config c{};
  constexpr std::size_t nlefs = 6;

  // clang-format off
  const std::array<Lef, nlefs> lefs1{
      construct_lef(95, 100, 0),
      construct_lef(101, 103, 0),
      construct_lef(102, 110, 0),
      construct_lef(104, 111, 0),
      construct_lef(105, 112, 0),
      construct_lef(102, 102, 1)
  };
  const std::array<Lef, nlefs> lefs2{
      construct_lef(95, 100, 0),
      construct_lef(101, 103, 0),
      construct_lef(102, 102, 1),
      construct_lef(102, 110, 0),
      construct_lef(104, 111, 0),
      construct_lef(105, 112, 0)
  };
  // clang-format on

  std::array<std::size_t, nlefs> rev_ranks{};
  std::array<std::size_t, nlefs> fwd_ranks{};

  const std::array<std::size_t, nlefs> rev_ranks_expected1{0, 1, 2, 5, 3, 4};
  const std::array<std::size_t, nlefs> fwd_ranks_expected1{0, 5, 1, 2, 3, 4};

  const std::array<std::size_t, nlefs> rev_ranks_expected2{0, 1, 3, 2, 4, 5};
  const std::array<std::size_t, nlefs> fwd_ranks_expected2{0, 2, 1, 3, 4, 5};

  Simulation::test_rank_lefs(lefs1, absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks), false,
                             true);
  for (std::size_t i = 0; i < lefs1.size(); ++i) {
    CHECK(rev_ranks[i] == rev_ranks_expected1[i]);
    CHECK(fwd_ranks[i] == fwd_ranks_expected1[i]);
  }

  Simulation::test_rank_lefs(lefs2, absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks), false,
                             true);
  for (std::size_t i = 0; i < lefs2.size(); ++i) {
    CHECK(rev_ranks[i] == rev_ranks_expected2[i]);
    CHECK(fwd_ranks[i] == fwd_ranks_expected2[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("LEFs ranking 002 - Fwd extr. unit tied", "[simulation][short]") {
  modle::Config c{};
  constexpr std::size_t nlefs = 6;

  // clang-format off
  const std::array<Lef, nlefs> lefs1{
      construct_lef(95, 100, 0),
      construct_lef(101, 104, 0),
      construct_lef(102, 110, 0),
      construct_lef(103, 111, 0),
      construct_lef(105, 112, 0),
      construct_lef(104, 104, 1)
  };
  const std::array<Lef, nlefs> lefs2{
      construct_lef(95, 100, 0),
      construct_lef(104, 104, 1),
      construct_lef(101, 104, 0),
      construct_lef(102, 110, 0),
      construct_lef(103, 111, 0),
      construct_lef(105, 112, 0)
  };
  // clang-format on

  std::array<std::size_t, nlefs> rev_ranks{};
  std::array<std::size_t, nlefs> fwd_ranks{};

  const std::array<std::size_t, nlefs> rev_ranks_expected1{0, 1, 2, 3, 5, 4};
  const std::array<std::size_t, nlefs> fwd_ranks_expected1{0, 5, 1, 2, 3, 4};

  const std::array<std::size_t, nlefs> rev_ranks_expected2{0, 2, 3, 4, 1, 5};
  const std::array<std::size_t, nlefs> fwd_ranks_expected2{0, 1, 2, 3, 4, 5};

  Simulation::test_rank_lefs(lefs1, absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks), false,
                             true);
  for (std::size_t i = 0; i < lefs1.size(); ++i) {
    CHECK(rev_ranks[i] == rev_ranks_expected1[i]);
    CHECK(fwd_ranks[i] == fwd_ranks_expected1[i]);
  }

  Simulation::test_rank_lefs(lefs2, absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks), false,
                             true);
  for (std::size_t i = 0; i < lefs2.size(); ++i) {
    CHECK(rev_ranks[i] == rev_ranks_expected2[i]);
    CHECK(fwd_ranks[i] == fwd_ranks_expected2[i]);
  }
}

}  // namespace modle::libmodle::test
