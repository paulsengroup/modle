#pragma once

#include <absl/types/span.h>  // for MakeSpan

#include <algorithm>         // for fill
#include <cassert>           // for assert
#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr
#include <cstddef>           // IWYU pragma: keep for size_t
#include <memory>            // for allocator_traits<>::value_type
#include <vector>            // for vector

#include "modle/common.hpp"                      // for Bp
#include "modle/config.hpp"                      // for Config
#include "modle/extrusion_barriers.hpp"          // for ExtrusionBarrier
#include "modle/extrusion_factors.hpp"           // for Lef, ExtrusionUnit
#include "modle/simulation.hpp"                  // for Simulation, Simulation::collision_t, seeder
#include "modle/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNING_PUSH

namespace modle::test::libmodle {

[[nodiscard]] inline modle::PRNG init_rand_eng(uint64_t seed) {
#ifdef USE_XOSHIRO
  seeder s(seed);
  PRNG rand_eng(s.generateSeedSequence<4>());
#else
  seeder s{seed};
  PRNG rand_eng{s};
#endif
  return rand_eng;
}

inline void check_that_lefs_are_sorted_by_idx(const std::vector<Lef>& lefs,
                                              const std::vector<std::size_t>& rev_ranks,
                                              const std::vector<std::size_t>& fwd_ranks) {
  CHECK(std::is_sorted(fwd_ranks.begin(), fwd_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
  }));
  CHECK(std::is_sorted(rev_ranks.begin(), rev_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
  }));
}

inline void require_that_lefs_are_sorted_by_idx(const std::vector<Lef>& lefs,
                                                const std::vector<std::size_t>& rev_ranks,
                                                const std::vector<std::size_t>& fwd_ranks) {
  REQUIRE(std::is_sorted(fwd_ranks.begin(), fwd_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
  }));
  REQUIRE(std::is_sorted(rev_ranks.begin(), rev_ranks.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
  }));
}

TEST_CASE("Bind LEFs 001", "[bind-lefs][simulation][short]") {
  const Chromosome chrom{{"chr1", 0, 1000}};
  constexpr auto nlefs = 10UL;
  std::vector<Lef> lefs(nlefs, Lef{});
  std::vector<std::size_t> rank1(nlefs), rank2(nlefs);
  std::iota(rank1.begin(), rank1.end(), 0);
  std::copy(rank1.begin(), rank1.end(), rank2.begin());
  boost::dynamic_bitset<> mask(nlefs);
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = init_rand_eng(8865403569783063175ULL);

  for (auto i = 0UL; i < nlefs; ++i) {
    mask[i] = std::bernoulli_distribution{0.50}(rand_eng);
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
  auto rand_eng = init_rand_eng(16044114709020280409ULL);
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
  auto rand_eng = init_rand_eng(17550568853244630438ULL);
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
  const std::vector<Lef> lefs{Lef{{5},  {25}, true},
                              Lef{{10}, {20}, true},
                              Lef{{90}, {90}, true}};

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
  const std::vector<Lef> lefs{Lef{{20},  {50}, true},
                              Lef{{60},  {60}, true},
                              Lef{{200}, {310}, true},
                              Lef{{220}, {300}, true},
                              Lef{{240}, {250}, true},
                              Lef{{125}, {305}, true}};

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

  CHECK(std::equal(rev_moves.begin(), rev_moves.end(), rev_moves_adjusted.begin()));
  CHECK(std::equal(fwd_moves.begin(), fwd_moves.end(), fwd_moves_adjusted.begin()));
}

TEST_CASE("Generate LEF moves 001", "[generate-lef-moves][simulation][short]") {
  const Chromosome chrom{{"chr1", 1000, 2000}};
  constexpr auto nlefs = 100UL;
  constexpr auto iters = 1000UL;
  const auto pos = (chrom.start_pos() + chrom.end_pos() + 1) / 2;
  std::vector<Lef> lefs(nlefs, Lef{{pos}, {pos}, true});
  std::vector<std::size_t> rev_ranks(nlefs);
  std::iota(rev_ranks.begin(), rev_ranks.end(), 0);
  std::vector<bp_t> rev_moves(nlefs, 0), fwd_moves(nlefs, 0);  // NOLINT
  auto fwd_ranks = rev_ranks;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = init_rand_eng(8312545934532053745ULL);
  auto c = Config{};
  c.bin_size = 500;  // NOLINT
  c.rev_extrusion_speed = c.bin_size;
  c.fwd_extrusion_speed = c.bin_size;
  c.rev_extrusion_speed_std = static_cast<double>(c.rev_extrusion_speed) * 0.2;  // NOLINT
  c.fwd_extrusion_speed_std = static_cast<double>(c.fwd_extrusion_speed) * 0.2;  // NOLINT

  for (auto i = 0UL; i < iters; ++i) {
    CHECK_NOTHROW(Simulation{c, false}.test_generate_moves(
        &chrom, absl::MakeConstSpan(lefs), absl::MakeConstSpan(rev_ranks),
        absl::MakeConstSpan(fwd_ranks), absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
        rand_eng));

    Simulation::test_rank_lefs(lefs, absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks), true);

    std::for_each(rev_moves.begin(), rev_moves.end(), [&, j = 0UL](const auto n) mutable {
      CHECK(lefs[j].rev_unit.pos() >= chrom.start_pos() + n);
      CHECK(lefs[j++].rev_unit.pos() < chrom.end_pos());
    });

    std::for_each(fwd_moves.begin(), fwd_moves.end(), [&, j = 0UL](const auto n) mutable {
      CHECK(lefs[j].fwd_unit.pos() + n < chrom.end_pos());
      CHECK(lefs[j++].fwd_unit.pos() >= chrom.start_pos());
    });

    /*
    fmt::print(
        stderr, "Prevented ~{} rev overflows\n",
        std::count_if(rev_moves.begin(), rev_moves.end(), [&, i = 0UL](const auto move) mutable {
          return lefs[i++].rev_unit.pos() - move == chrom.start_pos();
        }));

    fmt::print(
        stderr, "Prevented ~{} fwd overflows\n",
        std::count_if(fwd_moves.begin(), fwd_moves.end(), [&, i = 0UL](const auto move) mutable {
          return lefs[i++].fwd_unit.pos() + move == chrom.end_pos() - 1;
        }));
    */
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
  auto rand_eng = init_rand_eng(9589236971571381257ULL);

  // clang-format off
  /*       __0__      __1__        2        __3__
  *        |   |      |   |        |        |   |
  *      <-0- -2->  <-4- -8->  <--14-->  <-18- -23->
  */
  std::vector<Lef> lefs{
     Lef{{0},  {2},  true},  // NOLINT
     Lef{{4},  {8},  true},  // NOLINT
     Lef{{14}, {14}, true},  // NOLINT
     Lef{{18}, {23}, true}   // NOLINT
  };
  // clang-format on
  const std::vector<bp_t> rev_ranks{0, 1, 2, 3};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2, 3};  // NOLINT

  const auto& REACHED_CHROM_BOUNDARY = Simulation::REACHED_CHROM_BOUNDARY;
  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);

  std::vector<bp_t> rev_moves{0, 3, 3, 3};  // NOLINT
  std::vector<bp_t> fwd_moves{2, 2, 2, 2};  // NOLINT
  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{REACHED_CHROM_BOUNDARY,
                                                             LEF_LEF_COLLISION,
                                                             NO_COLLISION,
                                                             LEF_LEF_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{LEF_LEF_COLLISION,
                                                             NO_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             NO_COLLISION};
  // clang-format on

  const std::vector<bp_t> rev_moves_expected{0, 1, 3, 2};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{0, 2, 1, 2};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(modle::Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask), rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
    CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
    CHECK(rev_moves[i] == rev_moves_expected[i]);
    CHECK(fwd_moves[i] == fwd_moves_expected[i]);
  }
}

TEST_CASE("Detect LEF-LEF collisions 002", "[lef-lef-collisions][simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 3;                   // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 2;                   // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  const Chromosome chrom{{"chr1", 0, 16}};
  constexpr auto nlefs = 4UL;
  (void)nlefs;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = init_rand_eng(9589236971571381257ULL);

  // clang-format off
  /*
   *                 _______1_______                _______3_______
   *        _______0_|_____        |      ________2_|______       |
   *        |        |    |        |      |         |     |       |
   *      <-0-     <-4-  -5->     -6->  <-9-     <-11-   -14->   -15->
   */

  std::vector<Lef> lefs{
      Lef{{0},  {5},  true},  // NOLINT
      Lef{{4},  {6},  true},  // NOLINT
      Lef{{9},  {14}, true},  // NOLINT
      Lef{{11}, {15}, true}   // NOLINT
  };
  // clang-format on
  const std::vector<bp_t> rev_ranks{0, 1, 2, 3};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2, 3};  // NOLINT

  const auto& REACHED_CHROM_BOUNDARY = Simulation::REACHED_CHROM_BOUNDARY;
  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  std::vector<bp_t> rev_moves{0, 3, 3, 3};  // NOLINT
  std::vector<bp_t> fwd_moves{2, 2, 1, 0};  // NOLINT
  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{REACHED_CHROM_BOUNDARY,
                                                             NO_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             REACHED_CHROM_BOUNDARY};
  // clang-format on

  const std::vector<bp_t> rev_moves_expected{0, 3, 2, 3};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{0, 0, 0, 0};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  REQUIRE_NOTHROW(modle::Simulation{c, false}.test_adjust_moves(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  CHECK_NOTHROW(modle::Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask), rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
    CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
    CHECK(rev_moves[i] == rev_moves_expected[i]);
    CHECK(fwd_moves[i] == fwd_moves_expected[i]);
    /*
    fmt::print(stderr, "i={}; rev_move={}/{}; fwd_move={}/{};\n", i, rev_moves[i],
                rev_moves_expected[i], fwd_moves[i], fwd_moves_expected[i]);
    fmt::print(stderr, "i={}; rev_status={}/{}; fwd_status={}/{};\n", i, rev_collision_mask[i],
               rev_collision_mask_expected[i], fwd_collision_mask[i],
               fwd_collision_mask_expected[i]);
    */
  }
}

TEST_CASE("Detect LEF-LEF collisions 003", "[lef-lef-collisions][simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 50;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 50;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  const Chromosome chrom{{"chr1", 100, 201}};
  constexpr auto nlefs = 3UL;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = init_rand_eng(3778138500607566040ULL);

  // clang-format off
  /*
   *                ___________________0___________________
   *                |       ___________1___________       |
   *                |       |       ___2___       |       |
   *                |       |       |     |       |       |
   *             <-120-  <-130-  <-140- -141->  -160->  -180->
   */

  std::vector<Lef> lefs{
      Lef{{120}, {180}, true},  // NOLINT
      Lef{{130}, {160}, true},  // NOLINT
      Lef{{140}, {141}, true}   // NOLINT
  };
  // clang-format on
  const std::vector<bp_t> rev_ranks{0, 1, 2};  // NOLINT
  const std::vector<bp_t> fwd_ranks{2, 1, 0};  // NOLINT

  const auto& REACHED_CHROM_BOUNDARY = Simulation::REACHED_CHROM_BOUNDARY;
  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  std::vector<bp_t> rev_moves{20, 30, 40};  // NOLINT
  std::vector<bp_t> fwd_moves{20, 40, 59};  // NOLINT
  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{REACHED_CHROM_BOUNDARY,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{REACHED_CHROM_BOUNDARY,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION};
  // clang-format on

  const std::vector<bp_t> rev_moves_expected{20, 29, 38};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{20, 39, 58};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  REQUIRE_NOTHROW(modle::Simulation{c, false}.test_adjust_moves(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  CHECK_NOTHROW(modle::Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
      absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask), rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
    CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
    CHECK(rev_moves[i] == rev_moves_expected[i]);
    CHECK(fwd_moves[i] == fwd_moves_expected[i]);

    /*
    fmt::print(stderr, "i={}; rev_move={}/{}; fwd_move={}/{};\n", i, rev_moves[i],
                rev_moves_expected[i], fwd_moves[i], fwd_moves_expected[i]);
    fmt::print(stderr, "i={}; rev_status={}/{}; fwd_status={}/{};\n", i, rev_collision_mask[i],
               rev_collision_mask_expected[i], fwd_collision_mask[i],
               fwd_collision_mask_expected[i]);
    */
  }
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
  auto rand_eng = init_rand_eng(245021575020225667ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      Lef{{0}, {1}, 1},  // NOLINT
      Lef{{3}, {4}, 1},  // NOLINT
      Lef{{5}, {5}, 1}   // NOLINT
  };
  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{2, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{4, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{8, 1.0, 0.0, '+'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 2};
  const std::vector<bp_t> fwd_ranks{0, 1, 2};

  std::vector<bp_t> rev_moves{0, 2, 2};  // NOLINT
  std::vector<bp_t> fwd_moves{2, 2, 2};  // NOLINT

  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             0UL,
                                                             1UL};
  const std::vector<std::size_t> fwd_collision_mask_expected{NO_COLLISION,
                                                             NO_COLLISION,
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

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
    CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
    CHECK(rev_moves[i] == rev_moves_expected[i]);
    CHECK(fwd_moves[i] == fwd_moves_expected[i]);

    /*
    fmt::print(stderr, "i={}; rev_move={}/{}; fwd_move={}/{};\n", i, rev_moves[i],
               rev_moves_expected[i], fwd_moves[i], fwd_moves_expected[i]);
    fmt::print(stderr, "i={}; rev_status={}/{}; fwd_status={}/{};\n", i, rev_collision_mask[i],
               rev_collision_mask_expected[i], fwd_collision_mask[i],
               fwd_collision_mask_expected[i]);
    */
  }
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
  auto rand_eng = init_rand_eng(502428610956409790ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      Lef{{0}, {1}, 1},  // NOLINT
      Lef{{3}, {4}, 1},  // NOLINT
      Lef{{5}, {5}, 1}   // NOLINT
  };
  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{2, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{4, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{8, 1.0, 0.0, '-'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 2};
  const std::vector<bp_t> fwd_ranks{0, 1, 2};

  std::vector<bp_t> rev_moves{0, 2, 2};  // NOLINT
  std::vector<bp_t> fwd_moves{2, 2, 2};  // NOLINT

  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             NO_COLLISION,
                                                             NO_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{0UL,
                                                             1UL,
                                                             NO_COLLISION};
  // clang-format on

  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{0, 2, 2};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{0, 0, 2};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
    CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
    CHECK(rev_moves[i] == rev_moves_expected[i]);
    CHECK(fwd_moves[i] == fwd_moves_expected[i]);

    /*
    fmt::print(stderr, "i={}; rev_move={}/{}; fwd_move={}/{};\n", i, rev_moves[i],
               rev_moves_expected[i], fwd_moves[i], fwd_moves_expected[i]);
    fmt::print(stderr, "i={}; rev_status={}/{}; fwd_status={}/{};\n", i, rev_collision_mask[i],
               rev_collision_mask_expected[i], fwd_collision_mask[i],
               fwd_collision_mask_expected[i]);
    */
  }
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
  auto rand_eng = init_rand_eng(17304119060106843975ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      Lef{{0}, {1}, 1},  // NOLINT
      Lef{{3}, {4}, 1},  // NOLINT
      Lef{{5}, {5}, 1}   // NOLINT
  };
  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{2, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{4, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{8, 1.0, 0.0, '+'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 2};
  const std::vector<bp_t> fwd_ranks{0, 1, 2};

  std::vector<bp_t> rev_moves{0, 2, 2};  // NOLINT
  std::vector<bp_t> fwd_moves{2, 2, 2};  // NOLINT

  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             0UL,
                                                             1UL};
  const std::vector<std::size_t> fwd_collision_mask_expected{0UL,
                                                             1UL,
                                                             NO_COLLISION};
  // clang-format on

  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{0, 0, 0};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{0, 0, 2};  // NOLINT

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
    CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
    CHECK(rev_moves[i] == rev_moves_expected[i]);
    CHECK(fwd_moves[i] == fwd_moves_expected[i]);

    /*
    fmt::print(stderr, "i={}; rev_move={}/{}; fwd_move={}/{};\n", i, rev_moves[i],
               rev_moves_expected[i], fwd_moves[i], fwd_moves_expected[i]);
    fmt::print(stderr, "i={}; rev_status={}/{}; fwd_status={}/{};\n", i, rev_collision_mask[i],
               rev_collision_mask_expected[i], fwd_collision_mask[i],
               fwd_collision_mask_expected[i]);
    */
  }
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
  auto rand_eng = init_rand_eng(4870652027555157985ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      Lef{{10}, {20}, 1},  // NOLINT
      Lef{{26}, {26}, 1},  // NOLINT
      Lef{{30}, {35}, 1},  // NOLINT
      Lef{{42}, {43}, 1},  // NOLINT
      Lef{{44}, {60}, 1}   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{25, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{27, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{40, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{46, 1.0, 0.0, '-'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 2, 3, 4};
  const std::vector<bp_t> fwd_ranks{0, 1, 2, 3, 4};

  std::vector<bp_t> rev_moves{5, 5, 5, 5, 5};  // NOLINT
  std::vector<bp_t> fwd_moves{5, 5, 5, 5, 5};  // NOLINT

  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             0UL,
                                                             1UL,
                                                             2UL,
                                                             NO_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{NO_COLLISION,
                                                             NO_COLLISION,
                                                             NO_COLLISION,
                                                             3UL,
                                                             NO_COLLISION};
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

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
    CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
    CHECK(rev_moves[i] == rev_moves_expected[i]);
    CHECK(fwd_moves[i] == fwd_moves_expected[i]);

    /*
    fmt::print(stderr, "i={}; rev_move={}/{}; fwd_move={}/{};\n", i, rev_moves[i],
               rev_moves_expected[i], fwd_moves[i], fwd_moves_expected[i]);
    fmt::print(stderr, "i={}; rev_status={}/{}; fwd_status={}/{};\n", i, rev_collision_mask[i],
               rev_collision_mask_expected[i], fwd_collision_mask[i],
               fwd_collision_mask_expected[i]);
    */
  }
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
  auto rand_eng = init_rand_eng(4870652027555157985ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      Lef{{10}, {20}, 1},  // NOLINT
      Lef{{26}, {26}, 1},  // NOLINT
      Lef{{30}, {35}, 1},  // NOLINT
      Lef{{42}, {43}, 1},  // NOLINT
      Lef{{44}, {60}, 1}   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{25, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{27, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{40, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{46, 1.0, 0.0, '-'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 2, 3, 4};
  const std::vector<bp_t> fwd_ranks{0, 1, 2, 3, 4};

  std::vector<bp_t> rev_moves{2, 2, 2, 2, 2};  // NOLINT
  std::vector<bp_t> fwd_moves{5, 5, 5, 5, 5};  // NOLINT

  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             0UL,
                                                             NO_COLLISION,
                                                             2UL,
                                                             NO_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{NO_COLLISION,
                                                             NO_COLLISION,
                                                             NO_COLLISION,
                                                             3UL,
                                                             NO_COLLISION};
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

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
    CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
    CHECK(rev_moves[i] == rev_moves_expected[i]);
    CHECK(fwd_moves[i] == fwd_moves_expected[i]);

    /*
    fmt::print(stderr, "i={}; rev_move={}/{}; fwd_move={}/{};\n", i, rev_moves[i],
               rev_moves_expected[i], fwd_moves[i], fwd_moves_expected[i]);
    fmt::print(stderr, "i={}; rev_status={}/{}; fwd_status={}/{};\n", i, rev_collision_mask[i],
               rev_collision_mask_expected[i], fwd_collision_mask[i],
               fwd_collision_mask_expected[i]);
    */
  }
}

}  // namespace modle::test::libmodle
