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

inline void apply_lef_lef_stalls_wrapper(
    const modle::Config& c, std::vector<Lef>& lefs,
    const std::vector<Simulation::collision_t>& rev_collision_buff,
    const std::vector<Simulation::collision_t>& fwd_collision_buff,
    const std::vector<std::size_t>& rev_lef_rank_buff,
    const std::vector<std::size_t>& fwd_lef_rank_buff) {
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  seeder s{1583769804295373920};
  PRNG rand_eng{s};

  modle::Simulation(c, false).test_apply_lef_lef_stalls(
      absl::MakeSpan(lefs), rev_collision_buff, fwd_collision_buff,
      absl::MakeSpan(rev_lef_rank_buff), absl::MakeSpan(fwd_lef_rank_buff), rand_eng);
  for (auto i = 0UL; i < lefs.size(); ++i) {
    /*
    if (((lefs[i].rev_unit.lef_lef_stalls() > 0) != (rev_collision_buff[i] > 0))) {
      fmt::print(stderr, "i={}; rev_lef_lef_stalls={}; rev_collisions={};\n", i,
                 lefs[i].rev_unit.lef_lef_stalls(), rev_collision_buff[i]);
    }
    if (((lefs[i].fwd_unit.lef_lef_stalls() > 0) != (fwd_collision_buff[i] > 0))) {
      fmt::print(stderr, "i={}; fwd_lef_lef_stalls={}; fwd_collisions={};\n", i,
                 lefs[i].fwd_unit.lef_lef_stalls(), fwd_collision_buff[i]);
    }
     */

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

  modle::Simulation(c, false).test_check_lef_bar_collisions(
      lefs, rev_lef_rank_buff, fwd_lef_rank_buff, barriers,
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
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  seeder s{437181240247933207};
  PRNG rand_eng{s};

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
  c.bin_size = 5;
  const std::size_t nlefs = 4;

  std::vector<Bp> fwd_lef_rank_buff{0, 1, 2, 3};
  std::vector<Bp> rev_lef_rank_buff{0, 1, 2, 3};

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

  modle::Simulation(c, false).test_check_lef_lef_collisions(
      lefs, rev_lef_rank_buff, fwd_lef_rank_buff, absl::MakeSpan(rev_collision_counter_buff),
      absl::MakeSpan(fwd_collision_counter_buff));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_counter_buff[i] == expected_collisions_rev[i]);
    CHECK(fwd_collision_counter_buff[i] == expected_collisions_fwd[i]);
  }
}

TEST_CASE("Detect LEF-LEF collision simple 002", "[simulation][short]") {
  modle::Config c;
  c.bin_size = 5;  // NOLINT
  const std::size_t nlefs = 4;

  std::vector<Bp> fwd_lef_rank_buff{0, 1, 2, 3};
  std::vector<Bp> rev_lef_rank_buff{0, 1, 2, 3};

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

  modle::Simulation(c, false).test_check_lef_lef_collisions(
      lefs, rev_lef_rank_buff, fwd_lef_rank_buff, absl::MakeSpan(rev_collision_counter_buff),
      absl::MakeSpan(fwd_collision_counter_buff));

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

  std::vector<Bp> fwd_lef_rank_buff{0, 1, 2, 3};
  std::vector<Bp> rev_lef_rank_buff{0, 1, 2, 3};

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
  c.bin_size = 5;
  c.probability_of_extrusion_unit_bypass = 0.05;
  const std::size_t nlefs = 4;

  std::vector<Bp> fwd_lef_rank_buff{0, 1, 2, 3};
  std::vector<Bp> rev_lef_rank_buff{0, 1, 2, 3};

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
  c.bin_size = 2;
  const std::size_t nlefs = 3;

  std::vector<Bp> fwd_lef_rank_buff{0, 1, 2};
  std::vector<Bp> rev_lef_rank_buff{0, 1, 2};

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

  std::vector<Bp> fwd_lef_rank_buff{0, 1, 2};
  std::vector<Bp> rev_lef_rank_buff{0, 1, 2};

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
  c.soft_stall_multiplier = 0.5;
  c.hard_stall_multiplier = 2.0;
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

}  // namespace modle::test::libmodle
