#pragma once

#include <absl/types/span.h>  // for MakeSpan

#include <algorithm>         // for fill
#include <cassert>           // for assert
#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr
#include <cstddef>           // IWYU pragma: keep for size_t
#include <memory>            // for allocator_traits<>::value_type
#include <vector>            // for vector

#include "modle/common.hpp"                      // for bp_t, PRNG, PRNG_t
#include "modle/config.hpp"                      // for Config
#include "modle/extrusion_barriers.hpp"          // for ExtrusionBarrier
#include "modle/extrusion_factors.hpp"           // for Lef, ExtrusionUnit
#include "modle/simulation.hpp"                  // for Simulation, Simulation::collision_t, seeder
#include "modle/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNING_PUSH

namespace modle::test::libmodle {

[[maybe_unused]] inline void print_debug_info(
    std::size_t i, absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
    absl::Span<const bp_t> rev_moves_expected, absl::Span<const bp_t> fwd_moves_expected,
    absl::Span<const std::size_t> rev_collision_mask,
    absl::Span<const std::size_t> rev_collision_mask_expected,
    absl::Span<const std::size_t> fwd_collision_mask,
    absl::Span<const std::size_t> fwd_collision_mask_expected) {
  fmt::print(stderr, "i={}; rev_move={}/{}; fwd_move={}/{};\n", i, rev_moves[i],
             rev_moves_expected[i], fwd_moves[i], fwd_moves_expected[i]);
  fmt::print(stderr, "i={}; rev_status={}/{}; fwd_status={}/{};\n", i, rev_collision_mask[i],
             rev_collision_mask_expected[i], fwd_collision_mask[i], fwd_collision_mask_expected[i]);
}

[[maybe_unused]] inline void check_simulation_result(
    absl::Span<const Lef> lefs, absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
    absl::Span<const bp_t> rev_moves_expected, absl::Span<const bp_t> fwd_moves_expected,
    absl::Span<const std::size_t> rev_collision_mask,
    absl::Span<const std::size_t> rev_collision_mask_expected,
    absl::Span<const std::size_t> fwd_collision_mask,
    absl::Span<const std::size_t> fwd_collision_mask_expected, bool print_debug_info_ = false) {
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
using EU = ExtrusionUnit;

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
  std::vector<std::size_t> rank1(nlefs), rank2(nlefs);  // NOLINT
  std::iota(rank1.begin(), rank1.end(), 0);
  std::copy(rank1.begin(), rank1.end(), rank2.begin());
  boost::dynamic_bitset<> mask(nlefs);
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(8865403569783063175ULL);

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
  const auto pos = (chrom.start_pos() + chrom.end_pos() + 1) / 2;
  std::vector<Lef> lefs(nlefs, Lef{EU{pos}, EU{pos}, true});
  std::vector<std::size_t> rev_ranks(nlefs);
  std::iota(rev_ranks.begin(), rev_ranks.end(), 0);
  std::vector<bp_t> rev_moves(nlefs, 0), fwd_moves(nlefs, 0);  // NOLINT
  auto fwd_ranks = rev_ranks;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(8312545934532053745ULL);
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
  const std::vector<bp_t> rev_ranks{0, 1, 2, 3};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2, 3};  // NOLINT

  const std::vector<ExtrusionBarrier> barriers;

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
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

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
  const std::vector<bp_t> rev_ranks{0, 1, 2, 3};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2, 3};  // NOLINT

  const std::vector<ExtrusionBarrier> barriers;

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
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
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
  const std::vector<bp_t> rev_ranks{0, 1, 2};  // NOLINT
  const std::vector<bp_t> fwd_ranks{2, 1, 0};  // NOLINT

  const std::vector<ExtrusionBarrier> barriers;

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

  const std::vector<bp_t> rev_ranks{0, 1, 2};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2};  // NOLINT

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

  const std::vector<bp_t> rev_ranks{0, 1, 2};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2};  // NOLINT

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

  const std::vector<bp_t> rev_ranks{0, 1, 2};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2};  // NOLINT

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

  const std::vector<bp_t> rev_ranks{0, 1, 2, 3, 4};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2, 3, 4};  // NOLINT

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

  const std::vector<bp_t> rev_ranks{0, 1, 2, 3, 4};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2, 3, 4};  // NOLINT

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

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Simulation 001", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 75;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 75;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_hard_collision_pblock = 1;             // NOLINT
  c.lef_soft_collision_pblock = 0;             // NOLINT
  const std::size_t nlefs = 7;
  const std::size_t nbarriers = 5;
  const Chromosome chrom{{"chr1", 0, 1000}};
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(10556020843759504871ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      Lef{EU{25},  EU{30},  1},  // NOLINT
      Lef{EU{150}, EU{150}, 1},  // NOLINT
      Lef{EU{200}, EU{350}, 1},  // NOLINT
      Lef{EU{230}, EU{399}, 1},  // NOLINT
      Lef{EU{425}, EU{425}, 1},  // NOLINT
      Lef{EU{625}, EU{800}, 1},  // NOLINT
      Lef{EU{650}, EU{650}, 1}   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{105, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{400, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{600, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{850, 1.0, 0.0, '+'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 2, 3, 4, 5, 6};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2, 3, 4, 6, 5};  // NOLINT

  std::vector<bp_t> rev_moves{25, 75, 75, 75, 75, 75, 75};  // NOLINT
  std::vector<bp_t> fwd_moves{75, 75, 75, 75, 75, 75, 75};  // NOLINT

  const auto& REACHED_CHROM_BOUNDARY = Simulation::REACHED_CHROM_BOUNDARY;
  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{REACHED_CHROM_BOUNDARY,
                                                             1UL,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             NO_COLLISION,
                                                             NO_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{0UL,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             2UL,
                                                             NO_COLLISION,
                                                             NO_COLLISION,
                                                             NO_COLLISION};
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{25, 44, 25, 54, 25, 75, 75};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{69, 24, 48,  0, 75, 75, 75};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  REQUIRE_NOTHROW(Simulation{c, false}.test_adjust_moves(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (rev_collision_mask[i] < nbarriers) {
      CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
      CHECK(rev_moves[i] == rev_moves_expected[i]);
    }
    if (fwd_collision_mask[i] < nbarriers) {
      CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
      CHECK(fwd_moves[i] == fwd_moves_expected[i]);
    }
  }

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Simulation 002", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 75;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 75;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_hard_collision_pblock = 1;             // NOLINT
  c.lef_soft_collision_pblock = 0;             // NOLINT
  const std::size_t nlefs = 7;
  const std::size_t nbarriers = 5;
  const Chromosome chrom{{"chr1", 0, 1000}};
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(9668495374385482848ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      Lef{EU{200}, EU{375}, 1},  // NOLINT
      Lef{EU{350}, EU{350}, 1},  // NOLINT
      Lef{EU{575}, EU{575}, 1},  // NOLINT
      Lef{EU{601}, EU{770}, 1},  // NOLINT
      Lef{EU{650}, EU{800}, 1},  // NOLINT
      Lef{EU{850}, EU{850}, 1},  // NOLINT
      Lef{EU{970}, EU{975}, 1}   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{150, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{400, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{600, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{895, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{900, 1.0, 0.0, '+'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 2, 3, 4, 5, 6};  // NOLINT
  const std::vector<bp_t> fwd_ranks{1, 0, 2, 3, 4, 5, 6};  // NOLINT

  std::vector<bp_t> rev_moves{75, 75, 75, 75, 75, 75, 75};  // NOLINT
  std::vector<bp_t> fwd_moves{75, 75, 75, 75, 75, 75, 24};  // NOLINT

  const auto& REACHED_CHROM_BOUNDARY = Simulation::REACHED_CHROM_BOUNDARY;
  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             NO_COLLISION,
                                                             NO_COLLISION,
                                                             2UL,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             4UL};
  const std::vector<std::size_t> fwd_collision_mask_expected{NO_COLLISION,
                                                             NO_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             3UL,
                                                             REACHED_CHROM_BOUNDARY};
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{75, 75, 75,  0, 48, 25, 69};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{75, 75, 25, 53, 24, 44, 24};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  REQUIRE_NOTHROW(Simulation{c, false}.test_adjust_moves(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (rev_collision_mask[i] < nbarriers) {
      CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
      CHECK(rev_moves[i] == rev_moves_expected[i]);
    }
    if (fwd_collision_mask[i] < nbarriers) {
      CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
      CHECK(fwd_moves[i] == fwd_moves_expected[i]);
    }
  }

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Simulation 003 - Soft collisions on", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 75;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 75;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_hard_collision_pblock = 1;             // NOLINT
  c.lef_soft_collision_pblock = 1;             // NOLINT
  const std::size_t nlefs = 7;
  const std::size_t nbarriers = 5;
  const Chromosome chrom{{"chr1", 0, 1000}};
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(919341715542527390ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      Lef{EU{200}, EU{375}, 1},  // NOLINT
      Lef{EU{350}, EU{350}, 1},  // NOLINT
      Lef{EU{575}, EU{575}, 1},  // NOLINT
      Lef{EU{601}, EU{770}, 1},  // NOLINT
      Lef{EU{650}, EU{800}, 1},  // NOLINT
      Lef{EU{850}, EU{850}, 1},  // NOLINT
      Lef{EU{970}, EU{975}, 1}   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{150, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{400, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{600, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{895, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{900, 1.0, 0.0, '+'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 2, 3, 4, 5, 6};  // NOLINT
  const std::vector<bp_t> fwd_ranks{1, 0, 2, 3, 4, 5, 6};  // NOLINT

  std::vector<bp_t> rev_moves{75, 75, 75, 75, 75, 75, 75};  // NOLINT
  std::vector<bp_t> fwd_moves{75, 75, 75, 75, 75, 75, 24};  // NOLINT

  const auto& REACHED_CHROM_BOUNDARY = Simulation::REACHED_CHROM_BOUNDARY;
  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{0UL,
                                                             NO_COLLISION,
                                                             NO_COLLISION,
                                                             2UL,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             4UL};
  const std::vector<std::size_t> fwd_collision_mask_expected{1UL,
                                                             LEF_LEF_COLLISION,
                                                             2UL,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             3UL,
                                                             REACHED_CHROM_BOUNDARY};
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{49, 75, 75,  0, 48, 25, 69};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{24, 48, 24, 53, 24, 44, 24};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  REQUIRE_NOTHROW(Simulation{c, false}.test_adjust_moves(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (rev_collision_mask[i] < nbarriers) {
      CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
      CHECK(rev_moves[i] == rev_moves_expected[i]);
    }
    if (fwd_collision_mask[i] < nbarriers) {
      CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
      CHECK(fwd_moves[i] == fwd_moves_expected[i]);
    }
  }

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Simulation 004 - Inactive barriers", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 75;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 75;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_hard_collision_pblock = 1;             // NOLINT
  c.lef_soft_collision_pblock = 0;             // NOLINT
  const std::size_t nlefs = 7;
  const std::size_t nbarriers = 5;
  const Chromosome chrom{{"chr1", 0, 1000}};
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(138440880292496584ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      Lef{EU{200}, EU{375}, 1},  // NOLINT
      Lef{EU{350}, EU{350}, 1},  // NOLINT
      Lef{EU{575}, EU{575}, 1},  // NOLINT
      Lef{EU{601}, EU{770}, 1},  // NOLINT
      Lef{EU{650}, EU{800}, 1},  // NOLINT
      Lef{EU{850}, EU{850}, 1},  // NOLINT
      Lef{EU{970}, EU{975}, 1}   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{150, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{400, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{600, 1.0, 0.0, '+'},
                                               ExtrusionBarrier{895, 1.0, 0.0, '-'},
                                               ExtrusionBarrier{900, 1.0, 0.0, '+'}};
  // clang-format on
  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::NOT_OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 2, 3, 4, 5, 6};  // NOLINT
  const std::vector<bp_t> fwd_ranks{1, 0, 2, 3, 4, 5, 6};  // NOLINT

  std::vector<bp_t> rev_moves{75, 75, 75, 75, 75, 75, 75};  // NOLINT
  std::vector<bp_t> fwd_moves{75, 75, 75, 75, 75, 75, 24};  // NOLINT

  const auto& REACHED_CHROM_BOUNDARY = Simulation::REACHED_CHROM_BOUNDARY;
  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             NO_COLLISION,
                                                             NO_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{NO_COLLISION,
                                                             NO_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             REACHED_CHROM_BOUNDARY};
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{75, 75, 75, 13, 61, 25, 60};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{75, 75, 12, 53, 24, 59, 24};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  REQUIRE_NOTHROW(Simulation{c, false}.test_adjust_moves(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (rev_collision_mask[i] < nbarriers) {
      CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
      CHECK(rev_moves[i] == rev_moves_expected[i]);
    }
    if (fwd_collision_mask[i] < nbarriers) {
      CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
      CHECK(fwd_moves[i] == fwd_moves_expected[i]);
    }
  }

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Simulation 005 - Multiple LEFs located at the same site", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 25;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 25;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_hard_collision_pblock = 1;             // NOLINT
  c.lef_soft_collision_pblock = 0;             // NOLINT
  constexpr std::size_t nlefs = 6;
  constexpr std::size_t nbarriers = 1;
  const Chromosome chrom{{"chr1", 0, 150}};  // NOLINT
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(9509334554644345044ULL);

  // clang-format off
  const std::vector<Lef> lefs{
      Lef{EU{30},  EU{50},  1},  // NOLINT
      Lef{EU{60},  EU{80},  1},  // NOLINT
      Lef{EU{60},  EU{80},  1},  // NOLINT
      Lef{EU{65},  EU{125}, 1},  // NOLINT
      Lef{EU{140}, EU{140}, 1},  // NOLINT
      Lef{EU{140}, EU{140}, 1}   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '-'}};

  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 2, 3, 4, 5};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 2, 3, 4, 5};  // NOLINT

  std::vector<bp_t> rev_moves{25, 25, 25, 25, 25, 25};  // NOLINT
  std::vector<bp_t> fwd_moves{25, 25, 25, 24,  8,  9};  // NOLINT
  // clang-format on

  const auto& REACHED_CHROM_BOUNDARY = Simulation::REACHED_CHROM_BOUNDARY;
  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             0UL,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             REACHED_CHROM_BOUNDARY};
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{25,  5,  4,  8,  8,  7};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{ 4, 18, 19,  6,  8,  9};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  REQUIRE_NOTHROW(Simulation{c, false}.test_adjust_moves(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (rev_collision_mask[i] < nbarriers) {
      CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
      CHECK(rev_moves[i] == rev_moves_expected[i]);
    }
    if (fwd_collision_mask[i] < nbarriers) {
      CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
      CHECK(fwd_moves[i] == fwd_moves_expected[i]);
    }
  }

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Simulation 006 - Few inactive LEFs", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 25;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 25;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_hard_collision_pblock = 1;             // NOLINT
  c.lef_soft_collision_pblock = 0;             // NOLINT
  constexpr std::size_t nlefs = 6;
  constexpr std::size_t nbarriers = 1;
  const Chromosome chrom{{"chr1", 0, 150}};  // NOLINT
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(8213068426516999476ULL);

  // clang-format off
  std::vector<Lef> lefs{
      Lef{EU{30},  EU{50},  1},  // NOLINT
      Lef{EU{60},  EU{80},  1},  // NOLINT
      Lef{EU{60},  EU{80},  1},  // NOLINT
      Lef{EU{65},  EU{125}, 1},  // NOLINT
      Lef{EU{140}, EU{140}, 1},  // NOLINT
      Lef{EU{140}, EU{140}, 1}   // NOLINT
  };

  lefs[2].release(); // NOLINT
  lefs[5].release(); // NOLINT

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '-'}};

  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1, 3, 4, 2, 5};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1, 3, 4, 2, 5};  // NOLINT

  std::vector<bp_t> rev_moves{25, 25,  0, 25, 25,  0};  // NOLINT
  std::vector<bp_t> fwd_moves{25, 25,  0, 24,  9,  0};  // NOLINT
  // clang-format on

  const auto& REACHED_CHROM_BOUNDARY = Simulation::REACHED_CHROM_BOUNDARY;
  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             NO_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             NO_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{LEF_LEF_COLLISION,
                                                             0UL,
                                                             NO_COLLISION,
                                                             LEF_LEF_COLLISION,
                                                             REACHED_CHROM_BOUNDARY,
                                                             NO_COLLISION};
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{25,  5,  0,  9,  8,  0};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{ 4, 19,  0,  6,  9,  0};  // NOLINT
  // clang-format on

  assert(lefs.size() == nlefs);                         // NOLINT
  assert(barriers.size() == nbarriers);                 // NOLINT
  assert(rev_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(fwd_collision_mask_expected.size() == nlefs);  // NOLINT
  assert(rev_collision_mask.size() == nlefs);           // NOLINT
  assert(fwd_collision_mask.size() == nlefs);           // NOLINT
  require_that_lefs_are_sorted_by_idx(lefs, rev_ranks, fwd_ranks);

  REQUIRE_NOTHROW(Simulation{c, false}.test_adjust_moves(
      &chrom, lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_bar_collisions(
      lefs, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves), barriers,
      barrier_mask, absl::MakeSpan(rev_collision_mask), absl::MakeSpan(fwd_collision_mask),
      rand_eng));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (rev_collision_mask[i] < nbarriers) {
      CHECK(rev_collision_mask[i] == rev_collision_mask_expected[i]);
      CHECK(rev_moves[i] == rev_moves_expected[i]);
    }
    if (fwd_collision_mask[i] < nbarriers) {
      CHECK(fwd_collision_mask[i] == fwd_collision_mask_expected[i]);
      CHECK(fwd_moves[i] == fwd_moves_expected[i]);
    }
  }

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Simulation 007 - LEF-LEF collision overrides LEF-BAR collision 1",
          "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 20;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 20;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_hard_collision_pblock = 1;             // NOLINT
  c.lef_soft_collision_pblock = 0;             // NOLINT
  constexpr std::size_t nlefs = 2;
  constexpr std::size_t nbarriers = 1;
  const Chromosome chrom{{"chr1", 0, 200}};  // NOLINT
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(4003028292286930594ULL);

  // clang-format off
  std::vector<Lef> lefs{
      Lef{EU{50},  EU{95},  1},  // NOLINT
      Lef{EU{110}, EU{150}, 1}   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '+'}};

  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1};  // NOLINT

  std::vector<bp_t> rev_moves{20, 20};  // NOLINT
  std::vector<bp_t> fwd_moves{20, 20};  // NOLINT
  // clang-format on

  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             LEF_LEF_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{LEF_LEF_COLLISION,
                                                             NO_COLLISION};
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{20, 7};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{ 7, 20}; // NOLINT
  // clang-format on

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

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

TEST_CASE("Simulation 008 - LEF-LEF collision overrides LEF-BAR collision 2",
          "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 20;                  // NOLINT
  c.rev_extrusion_speed_std = 0;               // NOLINT
  c.fwd_extrusion_speed = 20;                  // NOLINT
  c.fwd_extrusion_speed_std = 0;               // NOLINT
  c.probability_of_extrusion_unit_bypass = 0;  // NOLINT
  c.lef_hard_collision_pblock = 1;             // NOLINT
  c.lef_soft_collision_pblock = 0;             // NOLINT
  constexpr std::size_t nlefs = 2;
  constexpr std::size_t nbarriers = 1;
  const Chromosome chrom{{"chr1", 0, 200}};  // NOLINT
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = modle::PRNG(12931236264094635263ULL);

  // clang-format off
  std::vector<Lef> lefs{
      Lef{EU{50},  EU{90},  1},  // NOLINT
      Lef{EU{105}, EU{150}, 1}   // NOLINT
  };

  const std::vector<ExtrusionBarrier> barriers{ExtrusionBarrier{100, 1.0, 0.0, '-'}};

  boost::dynamic_bitset<> barrier_mask;
  barrier_mask.resize(nbarriers, static_cast<bool>(CTCF::State::OCCUPIED));

  const std::vector<bp_t> rev_ranks{0, 1};  // NOLINT
  const std::vector<bp_t> fwd_ranks{0, 1};  // NOLINT

  std::vector<bp_t> rev_moves{20, 20};  // NOLINT
  std::vector<bp_t> fwd_moves{20, 20};  // NOLINT
  // clang-format on

  const auto& LEF_LEF_COLLISION = Simulation::LEF_LEF_COLLISION;
  const auto& NO_COLLISION = Simulation::NO_COLLISION;  // clang-format off
  const std::vector<std::size_t> rev_collision_mask_expected{NO_COLLISION,
                                                             LEF_LEF_COLLISION};
  const std::vector<std::size_t> fwd_collision_mask_expected{LEF_LEF_COLLISION,
                                                             NO_COLLISION};
  std::vector<std::size_t> rev_collision_mask(nlefs, NO_COLLISION);
  std::vector<std::size_t> fwd_collision_mask(nlefs, NO_COLLISION);

  const std::vector<bp_t> rev_moves_expected{20,  7};  // NOLINT
  const std::vector<bp_t> fwd_moves_expected{ 7, 20};  // NOLINT
  // clang-format on

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

  CHECK_NOTHROW(Simulation{c, false}.test_check_lef_lef_collisions(
      &chrom, lefs, barriers, rev_ranks, fwd_ranks, absl::MakeSpan(rev_moves),
      absl::MakeSpan(fwd_moves), absl::MakeSpan(rev_collision_mask),
      absl::MakeSpan(fwd_collision_mask), rand_eng));

  check_simulation_result(lefs, rev_moves, fwd_moves, rev_moves_expected, fwd_moves_expected,
                          rev_collision_mask, rev_collision_mask_expected, fwd_collision_mask,
                          fwd_collision_mask_expected);
}

}  // namespace modle::test::libmodle
