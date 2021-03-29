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
  CHECK_NOTHROW(Simulation{c, false}.test_clamp_moves(
      &chrom, absl::MakeConstSpan(lefs), absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks),
      absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

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
  CHECK_NOTHROW(Simulation{c, false}.test_clamp_moves(
      &chrom, absl::MakeConstSpan(lefs), absl::MakeSpan(rev_ranks), absl::MakeSpan(fwd_ranks),
      absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves)));

  CHECK(std::equal(rev_moves.begin(), rev_moves.end(), rev_moves_adjusted.begin()));
  CHECK(std::equal(fwd_moves.begin(), fwd_moves.end(), fwd_moves_adjusted.begin()));
}

TEST_CASE("Generate LEF moves 001", "[generate-lef-moves][simulation][short]") {
  const Chromosome chrom{{"chr1", 1000, 2000}};
  constexpr auto nlefs = 100UL;
  constexpr auto iters = 1000UL;
  std::vector<Lef> lefs(nlefs, Lef{{chrom.start_pos()}, {chrom.end_pos() - 1}, true});
  std::vector<std::size_t> rev_ranks(nlefs);
  std::iota(rev_ranks.begin(), rev_ranks.end(), 0);
  std::vector<bp_t> rev_moves(nlefs, 0), fwd_moves(nlefs, 0);  // NOLINT
  auto fwd_ranks = rev_ranks;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  auto rand_eng = init_rand_eng(8312545934532053745ULL);
  auto c = Config{};
  c.rev_extrusion_speed_std = 0.01;
  c.fwd_extrusion_speed_std = 0.01;

  for (auto i = 0UL; i < iters; ++i) {
    CHECK_NOTHROW(Simulation{c, false}.test_generate_moves(
        &chrom, absl::MakeConstSpan(lefs), absl::MakeConstSpan(rev_ranks),
        absl::MakeConstSpan(fwd_ranks), absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves),
        rand_eng));

    Simulation::test_rank_lefs(absl::MakeConstSpan(lefs), absl::MakeSpan(rev_ranks),
                               absl::MakeSpan(fwd_ranks), true);

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

/*

inline void apply_lef_lef_stalls_wrapper(
    const modle::Config& c, std::vector<Lef>& lefs,
    const std::vector<Simulation::collision_t>& rev_collision_buff,
    const std::vector<Simulation::collision_t>& fwd_collision_buff,
    const std::vector<std::size_t>& rev_lef_rank_buff,
    const std::vector<std::size_t>& fwd_lef_rank_buff) {
#ifdef USE_XOSHIRO
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  seeder s(1583769804295373920);
  PRNG rand_eng(s.generateSeedSequence<4>());
#else
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  seeder s{1583769804295373920};
  PRNG rand_eng{s};
#endif

  modle::Simulation(c, false).test_apply_lef_lef_stalls(
      absl::MakeSpan(lefs), rev_collision_buff, fwd_collision_buff,
      absl::MakeSpan(rev_lef_rank_buff), absl::MakeSpan(fwd_lef_rank_buff), rand_eng);
  for (auto i = 0UL; i < lefs.size(); ++i) {

    if (((lefs[i].rev_unit.lef_lef_stalls() > 0) != (rev_collision_buff[i] > 0))) {
      fmt::print(stderr, "i={}; rev_lef_lef_stalls={}; rev_collisions={};\n", i,
                 lefs[i].rev_unit.lef_lef_stalls(), rev_collision_buff[i]);
    }
    if (((lefs[i].fwd_unit.lef_lef_stalls() > 0) != (fwd_collision_buff[i] > 0))) {
      fmt::print(stderr, "i={}; fwd_lef_lef_stalls={}; fwd_collisions={};\n", i,
                 lefs[i].fwd_unit.lef_lef_stalls(), fwd_collision_buff[i]);
    }


    CHECK(((lefs[i].rev_unit.lef_lef_stalls() > 0) == (rev_collision_buff[i] > 0)));
    CHECK(((lefs[i].fwd_unit.lef_lef_stalls() > 0) == (fwd_collision_buff[i] > 0)));
  }
}

inline void check_lef_bar_collisions_wrapper(
    const modle::Config& c, const std::vector<Lef>& lefs,
    const std::vector<ExtrusionBarrier>& barriers,
    const std::vector<Simulation::collision_t>& expected_collisions_rev,
    const std::vector<Simulation::collision_t>& expected_collisions_fwd,
    const std::vector<std::size_t>& rev_lef_rank_buff,
    const std::vector<std::size_t>& fwd_lef_rank_buff) {
  const auto& nlefs = lefs.size();

  std::vector<std::size_t> fwd_collision_counter_buff(nlefs, 0);
  std::vector<std::size_t> rev_collision_counter_buff(nlefs, 0);
  std::vector<bp_t> rev_moves(nlefs, c.rev_extrusion_speed);
  std::vector<bp_t> fwd_moves(nlefs, c.fwd_extrusion_speed);

  modle::Simulation(c, false).test_check_lef_bar_collisions(
      lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_moves, fwd_moves, barriers,
      absl::MakeSpan(rev_collision_counter_buff), absl::MakeSpan(fwd_collision_counter_buff));

  // fmt::print(stderr, "rev_collisions={}\n", absl::StrJoin(rev_collision_counter_buff, ", "));
  // fmt::print(stderr, "fwd_collisions={}\n", absl::StrJoin(fwd_collision_counter_buff, ", "));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_counter_buff[i] == expected_collisions_rev[i]);
    CHECK(fwd_collision_counter_buff[i] == expected_collisions_fwd[i]);
  }
}

inline void apply_lef_bar_stalls_wrapper(
    const modle::Config& c, std::vector<Lef>& lefs, const std::vector<ExtrusionBarrier>& barriers,
    const std::vector<Simulation::collision_t>& rev_collision_buff,
    const std::vector<Simulation::collision_t>& fwd_collision_buff) {
#ifdef USE_XOSHIRO
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  seeder s(437181240247933207);
  PRNG rand_eng(s.generateSeedSequence<4>());
#else
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  seeder s{437181240247933207};
  PRNG rand_eng{s};
#endif

  modle::Simulation(c, false).test_apply_lef_bar_stalls(absl::MakeSpan(lefs), rev_collision_buff,
                                                        fwd_collision_buff, barriers, rand_eng);

  // TODO Figure out a way to test hard stalls
  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(((lefs[i].rev_unit.lef_bar_stalls() > 0) == (rev_collision_buff[i] != -1UL)));
    CHECK(((lefs[i].fwd_unit.lef_bar_stalls() > 0) == (fwd_collision_buff[i] != -1UL)));
  }
}

TEST_CASE("Detect LEF-LEF collision simple 001", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 3;      // NOLINT
  c.rev_extrusion_speed_std = 0;  // NOLINT
  c.fwd_extrusion_speed = 2;      // NOLINT
  c.fwd_extrusion_speed_std = 0;  // NOLINT
  const std::size_t nlefs = 4;

  std::vector<bp_t> fwd_lef_rank_buff{0, 1, 2, 3};
  std::vector<bp_t> rev_lef_rank_buff{0, 1, 2, 3};

  // clang-format off
  std::vector<Lef> lefs{
     Lef{{0}, {2}, 1},    // NOLINT
     Lef{{4}, {8}, 1},    // NOLINT
     Lef{{14}, {14}, 1},  // NOLINT
     Lef{{18}, {23}, 1}   // NOLINT
  };

  const std::vector<std::size_t> expected_collisions_rev{0, 1, 0, 1};
  const std::vector<std::size_t> expected_collisions_fwd{1, 0, 1, 0};
  (void)nlefs;
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  std::vector<std::size_t> fwd_collision_counter_buff(nlefs, 0);
  std::vector<std::size_t> rev_collision_counter_buff(nlefs, 0);
  std::vector<bp_t> rev_moves(nlefs, c.rev_extrusion_speed);
  std::vector<bp_t> fwd_moves(nlefs, c.fwd_extrusion_speed);

  modle::Simulation(c, false).test_check_lef_lef_collisions(
      lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_moves, fwd_moves,
      absl::MakeSpan(rev_collision_counter_buff), absl::MakeSpan(fwd_collision_counter_buff));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_counter_buff[i] == expected_collisions_rev[i]);
    CHECK(fwd_collision_counter_buff[i] == expected_collisions_fwd[i]);
  }
}

TEST_CASE("Detect LEF-LEF collision simple 002", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 3;      // NOLINT
  c.rev_extrusion_speed_std = 0;  // NOLINT
  c.fwd_extrusion_speed = 2;      // NOLINT
  c.fwd_extrusion_speed_std = 0;  // NOLINT
  const std::size_t nlefs = 4;

  std::vector<bp_t> fwd_lef_rank_buff{0, 1, 2, 3};
  std::vector<bp_t> rev_lef_rank_buff{0, 1, 2, 3};

  // clang-format off
  std::vector<Lef> lefs{ // Lef{pos, lifetime, rank}
      Lef{{1}, {5}, 1},   // NOLINT
      Lef{{4}, {6}, 1},   // NOLINT
      Lef{{9}, {14}, 1},  // NOLINT
      Lef{{11}, {15}, 1}  // NOLINT
  };

  const std::vector<std::size_t> expected_collisions_rev{0, 0, 2, 1};
  const std::vector<std::size_t> expected_collisions_fwd{1, 2, 0, 0};
  (void)nlefs;
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  std::vector<std::size_t> fwd_collision_counter_buff(nlefs, 0);
  std::vector<std::size_t> rev_collision_counter_buff(nlefs, 0);
  std::vector<bp_t> fwd_moves(nlefs, c.fwd_extrusion_speed);
  std::vector<bp_t> rev_moves(nlefs, c.rev_extrusion_speed);

  modle::Simulation(c, false).test_check_lef_lef_collisions(
      lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_moves, fwd_moves,
      absl::MakeSpan(rev_collision_counter_buff), absl::MakeSpan(fwd_collision_counter_buff));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_counter_buff[i] == expected_collisions_rev[i]);
    CHECK(fwd_collision_counter_buff[i] == expected_collisions_fwd[i]);
  }
}

TEST_CASE("Apply LEF-LEF stalls simple 001", "[simulation][short]") {
  modle::Config c;
  c.bin_size = 5;                                 // NOLINT
  c.probability_of_extrusion_unit_bypass = 0.05;  // NOLINT
  const std::size_t nlefs = 4;

  std::vector<bp_t> fwd_lef_rank_buff{0, 1, 2, 3};
  std::vector<bp_t> rev_lef_rank_buff{0, 1, 2, 3};

  // clang-format off
  std::vector<Lef> lefs{
      Lef{{0}, {2}, 1},    // NOLINT
      Lef{{4}, {8}, 1},    // NOLINT
      Lef{{14}, {14}, 1},  // NOLINT
      Lef{{18}, {23}, 1}   // NOLINT
  };

  const std::vector<std::size_t> expected_collisions_rev{0, 1, 0, 1};
  const std::vector<std::size_t> expected_collisions_fwd{1, 0, 1, 0};
  (void)nlefs;
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  apply_lef_lef_stalls_wrapper(c, lefs, expected_collisions_rev, expected_collisions_fwd,
                               rev_lef_rank_buff, fwd_lef_rank_buff);
}

TEST_CASE("Apply LEF-LEF stalls simple 002", "[simulation][short]") {
  modle::Config c;
  c.bin_size = 5;  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  c.probability_of_extrusion_unit_bypass =
      0.05;  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  const std::size_t nlefs = 4;

  std::vector<bp_t> fwd_lef_rank_buff{0, 1, 2, 3};
  std::vector<bp_t> rev_lef_rank_buff{0, 1, 2, 3};

  // clang-format off
  std::vector<Lef> lefs{
      Lef{{1}, {5}, 1},   // NOLINT
      Lef{{4}, {6}, 1},   // NOLINT
      Lef{{9}, {14}, 1},  // NOLINT
      Lef{{11}, {15}, 1}  // NOLINT
  };

  const std::vector<std::size_t> expected_collisions_rev{0, 0, 2, 1};
  const std::vector<std::size_t> expected_collisions_fwd{1, 2, 0, 0};
  (void)nlefs;
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  apply_lef_lef_stalls_wrapper(c, lefs, expected_collisions_rev, expected_collisions_fwd,
                               rev_lef_rank_buff, fwd_lef_rank_buff);
}

TEST_CASE("Detect LEF-BAR collisions simple 001", "[simulation][short]") {
  modle::Config c;
  c.rev_extrusion_speed = 1;      // NOLINT
  c.rev_extrusion_speed_std = 0;  // NOLINT
  c.fwd_extrusion_speed = 1;      // NOLINT
  c.fwd_extrusion_speed_std = 0;  // NOLINT
  const std::size_t nlefs = 3;

  std::vector<bp_t> fwd_lef_rank_buff{0, 1, 2};
  std::vector<bp_t> rev_lef_rank_buff{0, 1, 2};

  // clang-format off
  std::vector<Lef> lefs{
      Lef{{0}, {1}, 1},  // NOLINT
      Lef{{3}, {4}, 1},  // NOLINT
      Lef{{5}, {5}, 1}   // NOLINT
  };

  std::vector<ExtrusionBarrier> barriers{ // ExtrusionBarrier{pos, prob_of_block, motif_direction}
      ExtrusionBarrier{2, 1.0, '+'},
      ExtrusionBarrier{4, 1.0, '+'},
      ExtrusionBarrier{8, 1.0, '+'}
  };

DISABLE_WARNING_PUSH
DISABLE_WARNING_SIGN_CONVERSION
  const std::vector<std::size_t> expected_collisions_rev{-1UL, 0, 1};
  const std::vector<std::size_t> expected_collisions_fwd{0, 1, -1UL};
DISABLE_WARNING_POP
  (void)nlefs;
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  check_lef_bar_collisions_wrapper(c, lefs, barriers, expected_collisions_rev,
                                   expected_collisions_fwd, rev_lef_rank_buff, fwd_lef_rank_buff);
}

TEST_CASE("Detect LEF-BAR collisions simple 002", "[simulation][short]") {
  modle::Config c;
  c.bin_size = 2;
  const std::size_t nlefs = 3;

  std::vector<bp_t> fwd_lef_rank_buff{0, 1, 2};
  std::vector<bp_t> rev_lef_rank_buff{0, 1, 2};

  // clang-format off
  std::vector<Lef> lefs{
      Lef{{0}, {1}, 1},  // NOLINT
      Lef{{3}, {4}, 1},  // NOLINT
      Lef{{5}, {8}, 1}   // NOLINT
  };

  std::vector<ExtrusionBarrier> barriers{ // ExtrusionBarrier{pos, prob_of_block, motif_direction}
      ExtrusionBarrier{2, 1.0, '+'},
      ExtrusionBarrier{4, 1.0, '+'},
      ExtrusionBarrier{8, 1.0, '+'}
  };

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_CONVERSION
  const std::vector<std::size_t> expected_collisions_rev{-1UL, 0, 1};
  const std::vector<std::size_t> expected_collisions_fwd{0, 1, 2};
  DISABLE_WARNING_POP
  (void)nlefs;
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  check_lef_bar_collisions_wrapper(c, lefs, barriers, expected_collisions_rev,
                                   expected_collisions_fwd, rev_lef_rank_buff, fwd_lef_rank_buff);
}

TEST_CASE("Apply LEF-BAR stalls simple 001", "[simulation][short]") {
  modle::Config c;
  c.bin_size = 2;
  const std::size_t nlefs = 3;

  // clang-format off
  std::vector<Lef> lefs{
      Lef{{0}, {1}, 1},  // NOLINT
      Lef{{3}, {4}, 1},  // NOLINT
      Lef{{5}, {5}, 1}   // NOLINT
  };

  std::vector<ExtrusionBarrier> barriers{ // ExtrusionBarrier{pos, prob_of_block, motif_direction}
      ExtrusionBarrier{2, 0.99, '+'},  // NOLINT
      ExtrusionBarrier{4, 0.99, '+'},  // NOLINT
      ExtrusionBarrier{8, 0.99, '+'}   // NOLINT
  };

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_CONVERSION
  (void)nlefs;
  const std::vector<std::size_t> collisions_rev{-1UL, 0, 1};
  const std::vector<std::size_t> collisions_fwd{0, 1, -1UL};
  DISABLE_WARNING_POP
  assert(lefs.size() == nlefs); // NOLINT
  assert(collisions_fwd.size() == nlefs); // NOLINT
  assert(collisions_rev.size() == nlefs); // NOLINT
  // clang-format on

  apply_lef_bar_stalls_wrapper(c, lefs, barriers, collisions_rev, collisions_fwd);
}

TEST_CASE("Apply LEF-BAR stalls (w hard-stall) simple 002", "[simulation][short]") {
  modle::Config c;
  c.bin_size = 2;
  c.soft_stall_multiplier =
      0.5;  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  c.hard_stall_multiplier =
      2.0;  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  const std::size_t nlefs = 3;

  // clang-format off
  std::vector<Lef> lefs{
      Lef{{0}, {1}, 1},  // NOLINT
      Lef{{3}, {4}, 1},  // NOLINT
      Lef{{5}, {8}, 1}   // NOLINT
  };

  std::vector<ExtrusionBarrier> barriers{ // ExtrusionBarrier{pos, prob_of_block, motif_direction}
      ExtrusionBarrier{2, 0.99, '+'},  // NOLINT
      ExtrusionBarrier{4, 0.99, '+'},  // NOLINT
      ExtrusionBarrier{8, 0.99, '-'}   // NOLINT
  };

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_CONVERSION
  const std::vector<std::size_t> collisions_rev{-1UL, 0, 1};
  const std::vector<std::size_t> collisions_fwd{0, 1, 2};
  DISABLE_WARNING_POP
  (void)nlefs;
  assert(lefs.size() == nlefs); // NOLINT
  assert(collisions_fwd.size() == nlefs); // NOLINT
  assert(collisions_rev.size() == nlefs); // NOLINT
  // clang-format on

  apply_lef_bar_stalls_wrapper(c, lefs, barriers, collisions_rev, collisions_fwd);
}
*/

}  // namespace modle::test::libmodle
