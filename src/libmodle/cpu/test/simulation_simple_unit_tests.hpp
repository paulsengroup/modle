#pragma once

#include <absl/types/span.h>  // for MakeSpan

#include <algorithm>         // for fill
#include <cassert>           // for assert
#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr
#include <cstddef>           // IWYU pragma: keep for size_t
#include <memory>            // for allocator_traits<>::value_type
#include <vector>            // for vector

#include "modle/chrom_sizes.hpp"
#include "modle/common/common.hpp"                      // for bp_t, PRNG, PRNG_t
#include "modle/common/config.hpp"                      // for Config
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNING_PUSH
#include "modle/extrusion_barriers.hpp"                 // for ExtrusionBarrier
#include "modle/extrusion_factors.hpp"                  // for Lef, ExtrusionUnit
#include "modle/simulation.hpp"  // for Simulation, Simulation::collision_t, seeder
#include "modle/test/common.hpp"

namespace modle::test::libmodle {

TEST_CASE("Bind LEFs 001", "[bind-lefs][simulation][short]") {
  const Chromosome chrom{{"chr1", 0, 1000}};
  constexpr auto nlefs = 10UL;
  std::vector<Lef> lefs(nlefs, Lef{});
  std::vector<std::size_t> rank1(nlefs), rank2(nlefs);  // NOLINT
  std::iota(rank1.begin(), rank1.end(), 0);
  std::copy(rank1.begin(), rank1.end(), rank2.begin());
  std::vector<size_t> mask(nlefs, 1);
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(8865403569783063175ULL);

  for (auto i = 0UL; i < nlefs; ++i) {
    mask[i] = boost::random::bernoulli_distribution{0.50}(rand_eng);
  }

  const auto c = Config{};
  CHECK_NOTHROW(Simulation{c, false}.test_bind_lefs(
      &chrom, absl::MakeSpan(lefs), absl::MakeSpan(rank1), absl::MakeSpan(rank2), mask, rand_eng));

  CHECK(lefs.size() == nlefs);
  CHECK(rank1.size() == nlefs);
  CHECK(rank2.size() == nlefs);
  CHECK(mask.size() == nlefs);
  check_that_lefs_are_sorted_by_idx(lefs, rank1, rank2);

  for (auto i = 0UL; i < nlefs; ++i) {
    const auto& lef = lefs[i];
    if (mask[i]) {
      CHECK(lef.is_bound());
      CHECK(lef.rev_unit.pos() == lef.fwd_unit.pos());
      CHECK(lef.rev_unit.pos() >= chrom.start_pos());
      CHECK(lef.rev_unit.pos() < chrom.end_pos());
    } else {
      CHECK(!lef.is_bound());
    }
  }
}

TEST_CASE("Bind LEFs 002 - No LEFs to bind", "[bind-lefs][simulation][short]") {
  const Chromosome chrom{{"chr1", 0, 1000}};
  std::vector<Lef> lefs;
  std::vector<std::size_t> rank1, rank2;  // NOLINT
  boost::dynamic_bitset<> mask1;
  std::vector<int> mask2;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(16044114709020280409ULL);
  const auto c = Config{};

  CHECK_NOTHROW(Simulation{c, false}.test_bind_lefs(
      &chrom, absl::MakeSpan(lefs), absl::MakeSpan(rank1), absl::MakeSpan(rank2), mask1, rand_eng));

  CHECK(lefs.empty());
  CHECK(rank1.empty());
  CHECK(rank2.empty());
  CHECK(mask1.empty());

  CHECK_NOTHROW(Simulation{c, false}.test_bind_lefs(
      &chrom, absl::MakeSpan(lefs), absl::MakeSpan(rank1), absl::MakeSpan(rank2), mask2, rand_eng));

  CHECK(lefs.empty());
  CHECK(rank1.empty());
  CHECK(rank2.empty());
  CHECK(mask2.empty());

  constexpr auto nlefs = 10UL;
  std::generate_n(std::back_inserter(lefs), nlefs, []() { return Lef{}; });
  rank1.resize(nlefs);
  std::iota(rank1.begin(), rank1.end(), 0);
  rank2 = rank1;
  mask1.resize(nlefs, 0);

  CHECK_NOTHROW(Simulation{c, false}.test_bind_lefs(
      &chrom, absl::MakeSpan(lefs), absl::MakeSpan(rank1), absl::MakeSpan(rank2), mask1, rand_eng));

  CHECK(lefs.size() == nlefs);
  CHECK(rank1.size() == nlefs);
  CHECK(rank2.size() == nlefs);
  CHECK(mask1.size() == nlefs);
  check_that_lefs_are_sorted_by_idx(lefs, rank1, rank2);
  CHECK(mask1.count() == 0);

  CHECK(std::all_of(lefs.begin(), lefs.end(), [](const auto& lef) { return !lef.is_bound(); }));
}

TEST_CASE("Bind LEFs 003 - Empty mask (i.e. bind all LEFs)", "[bind-lefs][simulation][short]") {
  const Chromosome chrom{{"chr1", 0, 1000}};
  constexpr auto nlefs = 10UL;
  std::vector<Lef> lefs(nlefs, Lef{});
  std::vector<std::size_t> rank1(nlefs), rank2(nlefs);
  std::iota(rank1.begin(), rank1.end(), 0);
  std::copy(rank1.begin(), rank1.end(), rank2.begin());
  boost::dynamic_bitset<> mask;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(17550568853244630438ULL);
  const auto c = Config{};

  CHECK_NOTHROW(Simulation{c, false}.test_bind_lefs(
      &chrom, absl::MakeSpan(lefs), absl::MakeSpan(rank1), absl::MakeSpan(rank2), mask, rand_eng));

  CHECK(lefs.size() == nlefs);
  CHECK(rank1.size() == nlefs);
  CHECK(rank2.size() == nlefs);
  CHECK(mask.empty());

  CHECK(std::all_of(lefs.begin(), lefs.end(), [](const auto& lef) { return lef.is_bound(); }));
  CHECK(std::all_of(lefs.begin(), lefs.end(),
                    [](const auto& lef) { return lef.rev_unit.pos() == lef.fwd_unit.pos(); }));
}

TEST_CASE("Adjust LEF moves 001", "[adjust-lef-moves][simulation][short]") {
  const Chromosome chrom{{"chr1", 0, 101}};
  // clang-format off
const std::vector<Lef> lefs{Lef{EU{5},  EU{25}, true},
                            Lef{EU{10}, EU{20}, true},
                            Lef{EU{90}, EU{90}, true}};

  // clang-format on
  const std::vector<std::size_t> rev_ranks{0, 1, 2};  // NOLINT
  const std::vector<std::size_t> fwd_ranks{1, 0, 2};  // NOLINT

  std::vector<bp_t> rev_moves{5, 10, 15};   // NOLINT
  std::vector<bp_t> fwd_moves{10, 20, 10};  // NOLINT

  const std::vector<bp_t> rev_moves_adjusted{5, 10, 15};   // NOLINT
  const std::vector<bp_t> fwd_moves_adjusted{15, 20, 10};  // NOLINT

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  auto c = Config{};
  c.rev_extrusion_speed_std = 1;
  c.fwd_extrusion_speed_std = 1;
  CHECK_NOTHROW(Simulation{c, false}.test_adjust_moves(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  CHECK(std::equal(rev_moves.begin(), rev_moves.end(), rev_moves_adjusted.begin()));
  CHECK(std::equal(fwd_moves.begin(), fwd_moves.end(), fwd_moves_adjusted.begin()));
}

TEST_CASE("Adjust LEF moves 002", "[adjust-lef-moves][simulation][short]") {
  const Chromosome chrom{{"chr1", 10, 400}};
  // clang-format off
const std::vector<Lef> lefs{Lef{EU{20},  EU{50}, true},
                            Lef{EU{60},  EU{60}, true},
                            Lef{EU{200}, EU{310}, true},
                            Lef{EU{220}, EU{300}, true},
                            Lef{EU{240}, EU{250}, true},
                            Lef{EU{125}, EU{305}, true}};

  // clang-format on
  const std::vector<std::size_t> rev_ranks{0, 1, 5, 2, 3, 4};  // NOLINT
  const std::vector<std::size_t> fwd_ranks{0, 1, 4, 3, 5, 2};  // NOLINT

  std::vector<bp_t> rev_moves{10, 10, 5, 25, 50, 10};  // NOLINT
  std::vector<bp_t> fwd_moves{25, 10, 5, 20, 20, 0};   // NOLINT

  const std::vector<bp_t> rev_moves_adjusted{10, 10, 10, 30, 50, 10};  // NOLINT
  const std::vector<bp_t> fwd_moves_adjusted{25, 15, 10, 20, 20, 15};  // NOLINT

  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  auto c = Config{};
  c.rev_extrusion_speed_std = 1;
  c.fwd_extrusion_speed_std = 1;
  CHECK_NOTHROW(Simulation{c, false}.test_adjust_moves(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  for (auto i = 0UL; i < rev_moves.size(); ++i) {
    CHECK(rev_moves[i] == rev_moves_adjusted[i]);
    CHECK(fwd_moves[i] == fwd_moves_adjusted[i]);
  }
}

TEST_CASE("Generate LEF moves 001", "[generate-lef-moves][simulation][short]") {
  const Chromosome chrom{{"chr1", 1000, 2000}};
  constexpr auto nlefs = 100UL;
  constexpr auto iters = 1000UL;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(8312545934532053745ULL);
  std::vector<Lef> lefs(nlefs);
  std::vector<std::size_t> rev_ranks(nlefs);
  std::vector<std::size_t> fwd_ranks(nlefs);
  std::vector<bp_t> rev_moves(nlefs, 0), fwd_moves(nlefs, 0);  // NOLINT
  auto c = Config{};
  c.bin_size = 500;  // NOLINT
  c.rev_extrusion_speed = c.bin_size;
  c.fwd_extrusion_speed = c.bin_size;
  c.rev_extrusion_speed_std = static_cast<double>(c.rev_extrusion_speed) * 0.2;  // NOLINT
  c.fwd_extrusion_speed_std = static_cast<double>(c.fwd_extrusion_speed) * 0.2;  // NOLINT

  for (auto i = 0UL; i < iters; ++i) {
    std::generate(lefs.begin(), lefs.end(), [&]() {
      const auto pos = boost::random::uniform_int_distribution<bp_t>{chrom.start_pos(),
                                                                     chrom.end_pos() - 1}(rand_eng);
      return Lef{EU{pos}, EU{pos}, true};
    });
    Simulation::test_rank_lefs(lefs, absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks), false,
                               true);

    CHECK_NOTHROW(Simulation{c, false}.test_generate_moves(
        &chrom, absl::MakeConstSpan(lefs), absl::MakeConstSpan(rev_ranks),
        absl::MakeConstSpan(fwd_ranks), absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
        rand_eng));

    std::for_each(rev_moves.begin(), rev_moves.end(), [&, j = 0UL](const auto n) mutable {
      CHECK(lefs[j].rev_unit.pos() >= chrom.start_pos() + n);
      CHECK(lefs[j++].rev_unit.pos() < chrom.end_pos());
    });

    std::for_each(fwd_moves.begin(), fwd_moves.end(), [&, j = 0UL](const auto n) mutable {
      CHECK(lefs[j].fwd_unit.pos() + n < chrom.end_pos());
      CHECK(lefs[j++].fwd_unit.pos() >= chrom.start_pos());
    });
  }
}

TEST_CASE("Detect LEF-LEF collisions 001", "[lef-lef-collisions][simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 3;                   // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 2;                   // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  const Chromosome chrom{{"chr1", 0, 30}};
  constexpr auto nlefs = 4UL;
  (void)nlefs;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(9589236971571381257ULL);

  // clang-format off
/*       __0__      __1__        2        __3__
*        |   |      |   |        |        |   |
*      <-0- -2->  <-4- -8->  <--14-->  <-18- -23->
*/
std::vector<Lef> lefs{
    Lef{EU{0},  EU{2},  true},  // NOLINT
    Lef{EU{4},  EU{8},  true},  // NOLINT
    Lef{EU{14}, EU{14}, true},  // NOLINT
    Lef{EU{18}, EU{23}, true}   // NOLINT
};
  // clang-format on
  const std::vector<size_t> rev_ranks{0, 1, 2, 3};  // NOLINT
  const std::vector<size_t> fwd_ranks{0, 1, 2, 3};  // NOLINT

  const std::vector<ExtrusionBarrier> barriers;
  const boost::dynamic_bitset<> barrier_mask;

  std::vector<collision_t> fwd_collision_mask(nlefs, NO_COLLISION);
  std::vector<collision_t> rev_collision_mask(nlefs, NO_COLLISION);

  std::vector<bp_t> rev_moves{0, 3, 3, 3};  // NOLINT
  std::vector<bp_t> fwd_moves{2, 2, 2, 2};  // NOLINT
  // clang-format off
const std::vector<collision_t> rev_collision_mask_expected{REACHED_CHROM_BOUNDARY,
                                                           0UL,
                                                           NO_COLLISION,
                                                           2UL};
const std::vector<collision_t> fwd_collision_mask_expected{1UL,
                                                           NO_COLLISION,
                                                           3UL,
                                                           NO_COLLISION};
  // clang-format on

  const std::vector<bp_t> rev_moves_expected{0, 1, 3, 2};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{0, 2, 1, 2};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(modle::Simulation{c, false}.test_detect_units_at_chrom_boundaries(
      &chrom, lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask)));

  CHECK_NOTHROW(modle::Simulation{c, false}.test_process_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Detect LEF-LEF collisions 002", "[lef-lef-collisions][simulation][short]") {
  modle::Config c;
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  const Chromosome chrom{{"chr1", 0, 16}};
  constexpr auto nlefs = 4UL;
  (void)nlefs;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(9589236971571381257ULL);

  // clang-format off
/*
 *                 _______1_______                _______3_______
 *        _______0_|_____        |      ________2_|______       |
 *        |        |    |        |      |         |     |       |
 *      <-0-     <-4-  -5->     -6->  <-9-     <-11-   -14->   -15->
 */

std::vector<Lef> lefs{
    Lef{EU{0},  EU{5},  true},  // NOLINT
    Lef{EU{4},  EU{6},  true},  // NOLINT
    Lef{EU{9},  EU{14}, true},  // NOLINT
    Lef{EU{11}, EU{15}, true}   // NOLINT
};
  // clang-format on
  const std::vector<size_t> rev_ranks{0, 1, 2, 3};  // NOLINT
  const std::vector<size_t> fwd_ranks{0, 1, 2, 3};  // NOLINT

  const std::vector<ExtrusionBarrier> barriers;
  const boost::dynamic_bitset<> barrier_mask;

  std::vector<collision_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<collision_t> fwd_collision_mask(nlefs, NO_COLLISION);

  std::vector<bp_t> rev_moves{0, 3, 3, 4};  // NOLINT
  std::vector<bp_t> fwd_moves{3, 2, 1, 0};  // NOLINT
  // clang-format off
const std::vector<collision_t> rev_collision_mask_expected{REACHED_CHROM_BOUNDARY,
                                                           NO_COLLISION,
                                                           1UL,
                                                           6UL};
const std::vector<collision_t> fwd_collision_mask_expected{5UL,
                                                           2UL,
                                                           REACHED_CHROM_BOUNDARY,
                                                           REACHED_CHROM_BOUNDARY};
  // clang-format on

  const std::vector<bp_t> rev_moves_expected{0, 3, 2, 3};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{0, 0, 1, 0};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(modle::Simulation{c, false}.test_detect_units_at_chrom_boundaries(
      &chrom, lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask)));

  CHECK_NOTHROW(modle::Simulation{c, false}.test_process_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Detect LEF-LEF collisions 003", "[lef-lef-collisions][simulation][short]") {
  modle::Config c;
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  const Chromosome chrom{{"chr1", 100, 201}};
  constexpr auto nlefs = 3UL;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(3778138500607566040ULL);

  // clang-format off
/*
 *                ___________________0___________________
 *                |       ___________1___________       |
 *                |       |       ___2___       |       |
 *                |       |       |     |       |       |
 *             <-120-  <-130-  <-140- -141->  -160->  -180->
 */

std::vector<Lef> lefs{
    Lef{EU{120}, EU{180}, true},  // NOLINT
    Lef{EU{130}, EU{160}, true},  // NOLINT
    Lef{EU{140}, EU{141}, true}   // NOLINT
};
  // clang-format on
  const std::vector<size_t> rev_ranks{0, 1, 2};  // NOLINT
  const std::vector<size_t> fwd_ranks{2, 1, 0};  // NOLINT

  const std::vector<ExtrusionBarrier> barriers;
  const boost::dynamic_bitset<> barrier_mask;

  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  std::vector<bp_t> rev_moves{20, 30, 40};  // NOLINT
  std::vector<bp_t> fwd_moves{20, 40, 59};  // NOLINT
  // clang-format off
const std::vector<collision_t> rev_collision_mask_expected{REACHED_CHROM_BOUNDARY,
                                                           3UL,
                                                           4UL};
const std::vector<collision_t> fwd_collision_mask_expected{REACHED_CHROM_BOUNDARY,
                                                           3UL,
                                                           4UL};
  // clang-format on

  const std::vector<bp_t> rev_moves_expected{20, 29, 38};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{20, 39, 57};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(modle::Simulation{c, false}.test_detect_units_at_chrom_boundaries(
      &chrom, lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask)));

  CHECK_NOTHROW(modle::Simulation{c, false}.test_process_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Detect LEF-BAR collisions 001 - wo soft collisions fwd CTCFs",
          "[lef-bar-collisions][simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 2;        // NOLINT
  c.rev_extrusion_speed_std = 0;    // NOLINT
  c.fwd_extrusion_speed = 2;        // NOLINT
  c.fwd_extrusion_speed_std = 0;    // NOLINT
  c.lef_hard_collision_pblock = 1;  // NOLINT
  c.lef_soft_collision_pblock = 0;  // NOLINT
  constexpr std::size_t nlefs = 3;
  constexpr std::size_t nbarriers = 3;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(245021575020225667ULL);

  // clang-format off
const std::vector<Lef> lefs{
    Lef{EU{0}, EU{1}, 1},  // NOLINT
    Lef{EU{3}, EU{4}, 1},  // NOLINT
    Lef{EU{5}, EU{5}, 1}   // NOLINT
};
const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{2, 1.0, 0.0, '+'},
                                             ExtrusionBarrier{4, 1.0, 0.0, '+'},
                                             ExtrusionBarrier{8, 1.0, 0.0, '+'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<size_t> rev_ranks{0, 1, 2};  // NOLINT
  const std::vector<size_t> fwd_ranks{0, 1, 2};  // NOLINT

  std::vector<bp_t> rev_moves{0, 2, 2};  // NOLINT
  std::vector<bp_t> fwd_moves{2, 2, 2};  // NOLINT

  const std::vector<collision_t> rev_collision_mask_expected{NO_COLLISION, 0UL, 1UL};
  const std::vector<collision_t> fwd_collision_mask_expected{NO_COLLISION, NO_COLLISION,
                                                             NO_COLLISION};
  // clang-format on

  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{0, 0, 0};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{2, 2, 2};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(Simulation{c, false}.test_detect_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, barriers, barrier_mask,
      absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask), rand_eng));
  CHECK_NOTHROW(Simulation::test_correct_moves_for_lef_bar_collisions(
      lefs, barriers, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), rev_collision_mask,
      fwd_collision_mask));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Detect LEF-BAR collisions 002 - wo soft-collisions rev CTCFs",
          "[lef-bar-collisions][simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 2;        // NOLINT
  c.rev_extrusion_speed_std = 0;    // NOLINT
  c.fwd_extrusion_speed = 2;        // NOLINT
  c.fwd_extrusion_speed_std = 0;    // NOLINT
  c.lef_hard_collision_pblock = 1;  // NOLINT
  c.lef_soft_collision_pblock = 0;  // NOLINT
  constexpr std::size_t nlefs = 3;
  constexpr std::size_t nbarriers = 3;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(502428610956409790ULL);

  // clang-format off
const std::vector<Lef> lefs{
    Lef{EU{0}, EU{1}, 1},  // NOLINT
    Lef{EU{3}, EU{4}, 1},  // NOLINT
    Lef{EU{5}, EU{5}, 1}   // NOLINT
};
const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{2, 1.0, 0.0, '-'},
                                             ExtrusionBarrier{4, 1.0, 0.0, '-'},
                                             ExtrusionBarrier{8, 1.0, 0.0, '-'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<size_t> rev_ranks{0, 1, 2};  // NOLINT
  const std::vector<size_t> fwd_ranks{0, 1, 2};  // NOLINT

  std::vector<bp_t> rev_moves{0, 2, 2};  // NOLINT
  std::vector<bp_t> fwd_moves{2, 2, 2};  // NOLINT

  const std::vector<collision_t> rev_collision_mask_expected{NO_COLLISION, NO_COLLISION,
                                                             NO_COLLISION};
  const std::vector<collision_t> fwd_collision_mask_expected{0UL, NO_COLLISION, NO_COLLISION};
  // clang-format on

  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{0, 2, 2};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{0, 2, 2};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(Simulation{c, false}.test_detect_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));
  CHECK_NOTHROW(Simulation::test_correct_moves_for_lef_bar_collisions(
      lefs, barriers, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), rev_collision_mask,
      fwd_collision_mask));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Detect LEF-BAR collisions 003 - w soft-collisions fwd CTCFs",
          "[lef-bar-collisions][simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 2;        // NOLINT
  c.rev_extrusion_speed_std = 0;    // NOLINT
  c.fwd_extrusion_speed = 2;        // NOLINT
  c.fwd_extrusion_speed_std = 0;    // NOLINT
  c.lef_hard_collision_pblock = 1;  // NOLINT
  c.lef_soft_collision_pblock = 1;  // NOLINT
  constexpr std::size_t nlefs = 3;
  constexpr std::size_t nbarriers = 3;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(17304119060106843975ULL);

  // clang-format off
const std::vector<Lef> lefs{
    Lef{EU{0}, EU{1}, 1},  // NOLINT
    Lef{EU{3}, EU{4}, 1},  // NOLINT
    Lef{EU{5}, EU{5}, 1}   // NOLINT
};
const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{2, 1.0, 0.0, '+'},
                                             ExtrusionBarrier{4, 1.0, 0.0, '+'},
                                             ExtrusionBarrier{8, 1.0, 0.0, '+'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<size_t> rev_ranks{0, 1, 2};  // NOLINT
  const std::vector<size_t> fwd_ranks{0, 1, 2};  // NOLINT

  std::vector<bp_t> rev_moves{0, 2, 2};  // NOLINT
  std::vector<bp_t> fwd_moves{2, 2, 2};  // NOLINT

  const std::vector<collision_t> rev_collision_mask_expected{NO_COLLISION, 0UL, 1UL};
  const std::vector<collision_t> fwd_collision_mask_expected{0UL, NO_COLLISION, NO_COLLISION};
  // clang-format on

  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{0, 0, 0};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{0, 2, 2};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(Simulation{c, false}.test_detect_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));
  CHECK_NOTHROW(Simulation::test_correct_moves_for_lef_bar_collisions(
      lefs, barriers, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), rev_collision_mask,
      fwd_collision_mask));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Detect LEF-BAR collisions 004 - wo soft-collisions mixed CTCFs",
          "[lef-bar-collisions][simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 5;        // NOLINT
  c.rev_extrusion_speed_std = 0;    // NOLINT
  c.fwd_extrusion_speed = 5;        // NOLINT
  c.fwd_extrusion_speed_std = 0;    // NOLINT
  c.lef_hard_collision_pblock = 1;  // NOLINT
  c.lef_soft_collision_pblock = 0;  // NOLINT
  const std::size_t nlefs = 5;
  const std::size_t nbarriers = 4;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(4870652027555157985ULL);

  // clang-format off
const std::vector<Lef> lefs{
    Lef{EU{10}, EU{20}, 1},  // NOLINT
    Lef{EU{26}, EU{26}, 1},  // NOLINT
    Lef{EU{30}, EU{35}, 1},  // NOLINT
    Lef{EU{42}, EU{43}, 1},  // NOLINT
    Lef{EU{44}, EU{60}, 1}   // NOLINT
};

const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{25, 1.0, 0.0, '+'},
                                             ExtrusionBarrier{27, 1.0, 0.0, '+'},
                                             ExtrusionBarrier{40, 1.0, 0.0, '+'},
                                             ExtrusionBarrier{46, 1.0, 0.0, '-'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<size_t> rev_ranks{0, 1, 2, 3, 4};  // NOLINT
  const std::vector<size_t> fwd_ranks{0, 1, 2, 3, 4};  // NOLINT

  std::vector<bp_t> rev_moves{5, 5, 5, 5, 5};  // NOLINT
  std::vector<bp_t> fwd_moves{5, 5, 5, 5, 5};  // NOLINT

  const std::vector<collision_t> rev_collision_mask_expected{NO_COLLISION, 0UL, 1UL, 2UL,
                                                             NO_COLLISION};
  const std::vector<collision_t> fwd_collision_mask_expected{NO_COLLISION, NO_COLLISION,
                                                             NO_COLLISION, 3UL, NO_COLLISION};
  // clang-format on

  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{5, 0, 2, 1, 5};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{5, 5, 5, 2, 5};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(Simulation{c, false}.test_detect_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));
  CHECK_NOTHROW(Simulation::test_correct_moves_for_lef_bar_collisions(
      lefs, barriers, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), rev_collision_mask,
      fwd_collision_mask));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Detect LEF-BAR collisions 005 - wo soft-collisions mixed CTCFs, different extr. speeds",
          "[lef-bar-collisions][simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 2;        // NOLINT
  c.rev_extrusion_speed_std = 0;    // NOLINT
  c.fwd_extrusion_speed = 5;        // NOLINT
  c.fwd_extrusion_speed_std = 0;    // NOLINT
  c.lef_hard_collision_pblock = 1;  // NOLINT
  c.lef_soft_collision_pblock = 0;  // NOLINT
  const std::size_t nlefs = 5;
  const std::size_t nbarriers = 4;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(4870652027555157985ULL);

  // clang-format off
const std::vector<Lef> lefs{
    Lef{EU{10}, EU{20}, 1},  // NOLINT
    Lef{EU{26}, EU{26}, 1},  // NOLINT
    Lef{EU{30}, EU{35}, 1},  // NOLINT
    Lef{EU{42}, EU{43}, 1},  // NOLINT
    Lef{EU{44}, EU{60}, 1}   // NOLINT
};

const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{25, 1.0, 0.0, '+'},
                                             ExtrusionBarrier{27, 1.0, 0.0, '+'},
                                             ExtrusionBarrier{40, 1.0, 0.0, '+'},
                                             ExtrusionBarrier{46, 1.0, 0.0, '-'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<size_t> rev_ranks{0, 1, 2, 3, 4};  // NOLINT
  const std::vector<size_t> fwd_ranks{0, 1, 2, 3, 4};  // NOLINT

  std::vector<bp_t> rev_moves{2, 2, 2, 2, 2};  // NOLINT
  std::vector<bp_t> fwd_moves{5, 5, 5, 5, 5};  // NOLINT

  const std::vector<collision_t> rev_collision_mask_expected{NO_COLLISION, 0UL, NO_COLLISION, 2UL,
                                                             NO_COLLISION};
  const std::vector<collision_t> fwd_collision_mask_expected{NO_COLLISION, NO_COLLISION,
                                                             NO_COLLISION, 3UL, NO_COLLISION};
  // clang-format on

  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{2, 0, 2, 1, 2};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{5, 5, 5, 2, 5};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(Simulation{c, false}.test_detect_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));
  CHECK_NOTHROW(Simulation::test_correct_moves_for_lef_bar_collisions(
      lefs, barriers, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), rev_collision_mask,
      fwd_collision_mask));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("LEFs ranking 001 - Rev extr. unit tied", "[simulation][short]") {
  modle::Config c{};
  constexpr size_t nlefs = 6;

  // clang-format off
  const std::vector<Lef> lefs1{
      Lef{EU{ 95},  EU{100},  1},  // NOLINT
      Lef{EU{101},  EU{103},  1},  // NOLINT
      Lef{EU{102},  EU{110},  1},  // NOLINT
      Lef{EU{104},  EU{111},  1},  // NOLINT
      Lef{EU{105},  EU{112},  1},  // NOLINT
      Lef{EU{102},  EU{102},  1}   // NOLINT
  };
  const std::vector<Lef> lefs2{
      Lef{EU{ 95},  EU{100},  1},  // NOLINT
      Lef{EU{101},  EU{103},  1},  // NOLINT
      Lef{EU{102},  EU{102},  1},  // NOLINT
      Lef{EU{102},  EU{110},  1},  // NOLINT
      Lef{EU{104},  EU{111},  1},  // NOLINT
      Lef{EU{105},  EU{112},  1}   // NOLINT
  };
  // clang-format on

  std::vector<size_t> rev_ranks(nlefs);
  std::vector<size_t> fwd_ranks(nlefs);

  const std::vector<size_t> rev_ranks_expected1{0, 1, 2, 5, 3, 4};  // NOLINT
  const std::vector<size_t> fwd_ranks_expected1{0, 5, 1, 2, 3, 4};  // NOLINT

  const std::vector<size_t> rev_ranks_expected2{0, 1, 3, 2, 4, 5};  // NOLINT
  const std::vector<size_t> fwd_ranks_expected2{0, 2, 1, 3, 4, 5};  // NOLINT

  assert(lefs1.size() == nlefs);                // NOLINT
  assert(lefs2.size() == nlefs);                // NOLINT
  assert(rev_ranks.size() == nlefs);            // NOLINT
  assert(fwd_ranks.size() == nlefs);            // NOLINT
  assert(rev_ranks_expected1.size() == nlefs);  // NOLINT
  assert(rev_ranks_expected1.size() == nlefs);  // NOLINT
  assert(fwd_ranks_expected2.size() == nlefs);  // NOLINT
  assert(fwd_ranks_expected2.size() == nlefs);  // NOLINT

  REQUIRE_NOTHROW(Simulation{c, false}.test_rank_lefs(lefs1, absl::MakeSpan(rev_ranks),
                                                      absl::MakeSpan(fwd_ranks), false, true));
  for (auto i = 0UL; i < lefs1.size(); ++i) {
    CHECK(rev_ranks[i] == rev_ranks_expected1[i]);
    CHECK(fwd_ranks[i] == fwd_ranks_expected1[i]);
  }

  REQUIRE_NOTHROW(Simulation{c, false}.test_rank_lefs(lefs2, absl::MakeSpan(rev_ranks),
                                                      absl::MakeSpan(fwd_ranks), false, true));
  for (auto i = 0UL; i < lefs2.size(); ++i) {
    CHECK(rev_ranks[i] == rev_ranks_expected2[i]);
    CHECK(fwd_ranks[i] == fwd_ranks_expected2[i]);
  }
}

TEST_CASE("LEFs ranking 002 - Fwd extr. unit tied", "[simulation][short]") {
  modle::Config c{};
  constexpr size_t nlefs = 6;

  // clang-format off
  const std::vector<Lef> lefs1{
      Lef{EU{ 95},  EU{100},  1},  // NOLINT
      Lef{EU{101},  EU{104},  1},  // NOLINT
      Lef{EU{102},  EU{110},  1},  // NOLINT
      Lef{EU{103},  EU{111},  1},  // NOLINT
      Lef{EU{105},  EU{112},  1},  // NOLINT
      Lef{EU{104},  EU{104},  1}   // NOLINT
  };
  const std::vector<Lef> lefs2{
      Lef{EU{ 95},  EU{100},  1},  // NOLINT
      Lef{EU{104},  EU{104},  1},  // NOLINT
      Lef{EU{101},  EU{104},  1},  // NOLINT
      Lef{EU{102},  EU{110},  1},  // NOLINT
      Lef{EU{103},  EU{111},  1},  // NOLINT
      Lef{EU{105},  EU{112},  1}   // NOLINT
  };
  // clang-format on

  std::vector<size_t> rev_ranks(nlefs);
  std::vector<size_t> fwd_ranks(nlefs);

  const std::vector<size_t> rev_ranks_expected1{0, 1, 2, 3, 5, 4};  // NOLINT
  const std::vector<size_t> fwd_ranks_expected1{0, 5, 1, 2, 3, 4};  // NOLINT

  const std::vector<size_t> rev_ranks_expected2{0, 2, 3, 4, 1, 5};  // NOLINT
  const std::vector<size_t> fwd_ranks_expected2{0, 1, 2, 3, 4, 5};  // NOLINT

  {
    assert(lefs1.size() == nlefs);                // NOLINT
    assert(lefs2.size() == nlefs);                // NOLINT
    assert(rev_ranks.size() == nlefs);            // NOLINT
    assert(fwd_ranks.size() == nlefs);            // NOLINT
    assert(rev_ranks_expected1.size() == nlefs);  // NOLINT
    assert(rev_ranks_expected1.size() == nlefs);  // NOLINT
    assert(fwd_ranks_expected2.size() == nlefs);  // NOLINT
    assert(fwd_ranks_expected2.size() == nlefs);  // NOLINT
  }

  REQUIRE_NOTHROW(Simulation{c, false}.test_rank_lefs(lefs1, absl::MakeSpan(rev_ranks),
                                                      absl::MakeSpan(fwd_ranks), false, true));
  for (auto i = 0UL; i < lefs1.size(); ++i) {
    CHECK(rev_ranks[i] == rev_ranks_expected1[i]);
    CHECK(fwd_ranks[i] == fwd_ranks_expected1[i]);
  }

  REQUIRE_NOTHROW(Simulation{c, false}.test_rank_lefs(lefs2, absl::MakeSpan(rev_ranks),
                                                      absl::MakeSpan(fwd_ranks), false, true));
  for (auto i = 0UL; i < lefs2.size(); ++i) {
    CHECK(rev_ranks[i] == rev_ranks_expected2[i]);
    CHECK(fwd_ranks[i] == fwd_ranks_expected2[i]);
  }
}

}  // namespace modle::test::libmodle
