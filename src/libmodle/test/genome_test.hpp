#pragma once

#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <catch2/catch.hpp>
#include <random>
#include <vector>

#include "modle/common.hpp"
#include "modle/config.hpp"
#include "modle/extrusion_factors.hpp"
#include "modle/genome.hpp"

namespace modle::test::libmodle {

inline void check_lef_lef_collisions_wrapper(
    const modle::config& c, const std::vector<Lef>& lefs,
    const std::vector<std::size_t>& expected_collisions_rev,
    const std::vector<std::size_t>& expected_collisions_fwd) {
  const auto& nlefs = lefs.size();
  std::vector<std::size_t> fwd_lef_rank_buff(nlefs);
  std::vector<std::size_t> rev_lef_rank_buff(nlefs);

  std::vector<std::size_t> fwd_collision_counter_buff(nlefs, 0);
  std::vector<std::size_t> rev_collision_counter_buff(nlefs, 0);

  std::transform(lefs.begin(), lefs.end(), rev_lef_rank_buff.begin(),
                 [](const auto& lef) { return lef.rev_unit.rank(); });
  std::transform(lefs.begin(), lefs.end(), fwd_lef_rank_buff.begin(),
                 [](const auto& lef) { return lef.fwd_unit.rank(); });

  modle::Genome g(c, false);
  g.test_check_lef_lef_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff,
                                  rev_collision_counter_buff, fwd_collision_counter_buff);

  // fmt::print(stderr, "rev_collisions={}\n", absl::StrJoin(rev_collision_counter_buff, ", "));
  // fmt::print(stderr, "fwd_collisions={}\n", absl::StrJoin(fwd_collision_counter_buff, ", "));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_counter_buff[i] == expected_collisions_rev[i]);
    CHECK(fwd_collision_counter_buff[i] == expected_collisions_fwd[i]);
  }
}

inline void apply_lef_lef_stalls_wrapper(const modle::config& c, std::vector<Lef>& lefs,
                                         const std::vector<std::size_t>& rev_collision_buff,
                                         const std::vector<std::size_t>& fwd_collision_buff) {
  const auto& nlefs = lefs.size();
  std::vector<std::size_t> fwd_lef_rank_buff(nlefs);
  std::vector<std::size_t> rev_lef_rank_buff(nlefs);

  std::transform(lefs.begin(), lefs.end(), rev_lef_rank_buff.begin(),
                 [](const auto& lef) { return lef.rev_unit.rank(); });
  std::transform(lefs.begin(), lefs.end(), fwd_lef_rank_buff.begin(),
                 [](const auto& lef) { return lef.fwd_unit.rank(); });

  modle::Genome g(c, false);
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  seeder s{1583769804295373920};
  PRNG rand_eng{s};

  g.test_apply_lef_lef_stalls(lefs, rev_collision_buff, fwd_collision_buff, rev_lef_rank_buff,
                              fwd_lef_rank_buff, rand_eng, c.probability_of_extrusion_unit_bypass);
  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(((lefs[i].rev_unit.lef_lef_stalls() > 0) == (rev_collision_buff[i] > 0)));
    CHECK(((lefs[i].fwd_unit.lef_lef_stalls() > 0) == (fwd_collision_buff[i] > 0)));
  }
}

inline void check_lef_bar_collisions_wrapper(
    const modle::config& c, const std::vector<Lef>& lefs,
    const std::vector<ExtrusionBarrier>& barriers,
    const std::vector<std::size_t>& expected_collisions_rev,
    const std::vector<std::size_t>& expected_collisions_fwd) {
  const auto& nlefs = lefs.size();
  std::vector<std::size_t> fwd_lef_rank_buff(nlefs);
  std::vector<std::size_t> rev_lef_rank_buff(nlefs);

  std::vector<std::size_t> fwd_collision_counter_buff(nlefs, 0);
  std::vector<std::size_t> rev_collision_counter_buff(nlefs, 0);

  std::transform(lefs.begin(), lefs.end(), rev_lef_rank_buff.begin(),
                 [](const auto& lef) { return lef.rev_unit.rank(); });
  std::transform(lefs.begin(), lefs.end(), fwd_lef_rank_buff.begin(),
                 [](const auto& lef) { return lef.fwd_unit.rank(); });

  modle::Genome g(c, false);
  g.test_check_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, barriers,
                                  rev_collision_counter_buff, fwd_collision_counter_buff);

  // fmt::print(stderr, "rev_collisions={}\n", absl::StrJoin(rev_collision_counter_buff, ", "));
  // fmt::print(stderr, "fwd_collisions={}\n", absl::StrJoin(fwd_collision_counter_buff, ", "));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(rev_collision_counter_buff[i] == expected_collisions_rev[i]);
    CHECK(fwd_collision_counter_buff[i] == expected_collisions_fwd[i]);
  }
}

inline void apply_lef_bar_stalls_wrapper(const modle::config& c, std::vector<Lef>& lefs,
                                         const std::vector<ExtrusionBarrier>& barriers,
                                         const std::vector<std::size_t>& rev_collision_buff,
                                         const std::vector<std::size_t>& fwd_collision_buff) {
  const auto& nlefs = lefs.size();
  std::vector<std::size_t> fwd_lef_rank_buff(nlefs);
  std::vector<std::size_t> rev_lef_rank_buff(nlefs);

  std::transform(lefs.begin(), lefs.end(), rev_lef_rank_buff.begin(),
                 [](const auto& lef) { return lef.rev_unit.rank(); });
  std::transform(lefs.begin(), lefs.end(), fwd_lef_rank_buff.begin(),
                 [](const auto& lef) { return lef.fwd_unit.rank(); });

  modle::Genome g(c, false);
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  seeder s{437181240247933207};
  PRNG rand_eng{s};

  g.test_apply_lef_bar_stalls(lefs, rev_collision_buff, fwd_collision_buff, barriers,
                              rev_lef_rank_buff, fwd_lef_rank_buff, rand_eng, 0.5, 2.0);

  // TODO Figure out a way to test hard stalls
  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(((lefs[i].rev_unit.lef_bar_stalls() > 0) == (rev_collision_buff[i] != -1UL)));
    CHECK(((lefs[i].fwd_unit.lef_bar_stalls() > 0) == (fwd_collision_buff[i] != -1UL)));
  }
}

TEST_CASE("Detect LEF-LEF collision simple 001", "[simulation][short]") {
  modle::config c;
  c.bin_size = 5;
  const std::size_t nlefs = 4;

  std::vector<Bp> fwd_lef_rank_buff(nlefs);
  std::vector<Bp> rev_lef_rank_buff(nlefs);

  // clang-format off
  std::vector<Lef> lefs{ // Lef{pos, lifetime, rank}
     Lef{{0, 0, 0}, {2, 0, 0}},
     Lef{{4, 0, 1}, {8, 0, 1}},
     Lef{{14, 0, 2}, {14, 0, 2}},
     Lef{{18, 0, 3}, {23, 0, 3}}
  };

  const std::vector<std::size_t> expected_collisions_rev{0, 1, 0, 1};
  const std::vector<std::size_t> expected_collisions_fwd{1, 0, 1, 0};
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  check_lef_lef_collisions_wrapper(c, lefs, expected_collisions_rev, expected_collisions_fwd);
}

TEST_CASE("Detect LEF-LEF collision simple 002", "[simulation][short]") {
  modle::config c;
  c.bin_size = 5;
  const std::size_t nlefs = 4;

  std::vector<Bp> fwd_lef_rank_buff(nlefs);
  std::vector<Bp> rev_lef_rank_buff(nlefs);

  // clang-format off
  std::vector<Lef> lefs{ // Lef{pos, lifetime, rank}
      Lef{{1, 0, 0}, {5, 0, 0}},
      Lef{{4, 0, 1}, {6, 0, 1}},
      Lef{{9, 0, 2}, {14, 0, 2}},
      Lef{{11, 0, 3}, {15, 0, 3}}
  };

  const std::vector<std::size_t> expected_collisions_rev{0, 0, 2, 1};
  const std::vector<std::size_t> expected_collisions_fwd{1, 2, 0, 0};
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  check_lef_lef_collisions_wrapper(c, lefs, expected_collisions_rev, expected_collisions_fwd);
}

TEST_CASE("Apply LEF-LEF stalls simple 001", "[simulation][short]") {
  modle::config c;
  c.bin_size = 5;
  c.probability_of_extrusion_unit_bypass = 0.05;
  const std::size_t nlefs = 4;

  std::vector<Bp> fwd_lef_rank_buff(nlefs);
  std::vector<Bp> rev_lef_rank_buff(nlefs);

  // clang-format off
  std::vector<Lef> lefs{ // Lef{pos, lifetime, rank}
      Lef{{0, 0, 0}, {2, 0, 0}},
      Lef{{4, 0, 1}, {8, 0, 1}},
      Lef{{14, 0, 2}, {14, 0, 2}},
      Lef{{18, 0, 3}, {23, 0, 3}}
  };

  const std::vector<std::size_t> expected_collisions_rev{0, 1, 0, 1};
  const std::vector<std::size_t> expected_collisions_fwd{1, 0, 1, 0};
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  apply_lef_lef_stalls_wrapper(c, lefs, expected_collisions_rev, expected_collisions_fwd);
}

TEST_CASE("Apply LEF-LEF stalls simple 002", "[simulation][short]") {
  modle::config c;
  c.bin_size = 5;
  c.probability_of_extrusion_unit_bypass = 0.05;
  const std::size_t nlefs = 4;

  std::vector<Bp> fwd_lef_rank_buff(nlefs);
  std::vector<Bp> rev_lef_rank_buff(nlefs);

  // clang-format off
  std::vector<Lef> lefs{ // Lef{pos, lifetime, rank}
      Lef{{1, 0, 0}, {5, 0, 0}},
      Lef{{4, 0, 1}, {6, 0, 1}},
      Lef{{9, 0, 2}, {14, 0, 2}},
      Lef{{11, 0, 3}, {15, 0, 3}}
  };

  const std::vector<std::size_t> expected_collisions_rev{0, 0, 2, 1};
  const std::vector<std::size_t> expected_collisions_fwd{1, 2, 0, 0};
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  apply_lef_lef_stalls_wrapper(c, lefs, expected_collisions_rev, expected_collisions_fwd);
}

TEST_CASE("Apply LEF-LEF stalls (reset old lef_lef stalls) simple 003", "[simulation][short]") {
  modle::config c;
  c.bin_size = 5;
  c.probability_of_extrusion_unit_bypass = 0.05;
  const std::size_t nlefs = 4;

  std::vector<Bp> fwd_lef_rank_buff(nlefs);
  std::vector<Bp> rev_lef_rank_buff(nlefs);

  // clang-format off
  std::vector<Lef> lefs{ // Lef{pos, lifetime, rank, nstalls_lef_lef, nstalls_lef_bar}
      Lef{{1, 0, 0, 10}, {5, 0, 0, 10}},
      Lef{{4, 0, 1, 10}, {6, 0, 1, 10}},
      Lef{{9, 0, 2, 10}, {14, 0, 2, 10}},
      Lef{{11, 0, 3, 10}, {15, 0, 3, 10}}
  };

  const std::vector<std::size_t> expected_collisions_rev{0, 0, 2, 1};
  const std::vector<std::size_t> expected_collisions_fwd{1, 2, 0, 0};
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  apply_lef_lef_stalls_wrapper(c, lefs, expected_collisions_rev, expected_collisions_fwd);
}

TEST_CASE("Detect LEF-BAR collisions simple 001", "[simulation][short]") {
  modle::config c;
  c.bin_size = 2;
  const std::size_t nlefs = 3;

  std::vector<Bp> fwd_lef_rank_buff(nlefs);
  std::vector<Bp> rev_lef_rank_buff(nlefs);

  // clang-format off
  std::vector<Lef> lefs{ // Lef{pos, lifetime, rank, nstalls_lef_lef, nstalls_lef_bar}
      Lef{{0, 0, 0}, {1, 0, 0}},
      Lef{{3, 0, 1}, {4, 0, 1}},
      Lef{{5, 0, 2}, {5, 0, 2}}
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
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  check_lef_bar_collisions_wrapper(c, lefs, barriers, expected_collisions_rev,
                                   expected_collisions_fwd);
}

TEST_CASE("Detect LEF-BAR collisions simple 002", "[simulation][short]") {
  modle::config c;
  c.bin_size = 2;
  const std::size_t nlefs = 3;

  std::vector<Bp> fwd_lef_rank_buff(nlefs);
  std::vector<Bp> rev_lef_rank_buff(nlefs);

  // clang-format off
  std::vector<Lef> lefs{ // Lef{pos, lifetime, rank, nstalls_lef_lef, nstalls_lef_bar}
      Lef{{0, 0, 0}, {1, 0, 0}},
      Lef{{3, 0, 1}, {4, 0, 1}},
      Lef{{5, 0, 2}, {8, 0, 2}}
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
  assert(lefs.size() == nlefs); // NOLINT
  assert(expected_collisions_fwd.size() == nlefs); // NOLINT
  assert(expected_collisions_rev.size() == nlefs); // NOLINT
  // clang-format on
  check_lef_bar_collisions_wrapper(c, lefs, barriers, expected_collisions_rev,
                                   expected_collisions_fwd);
}

TEST_CASE("Apply LEF-BAR stalls simple 001", "[simulation][short]") {
  modle::config c;
  c.bin_size = 2;
  const std::size_t nlefs = 3;

  std::vector<Bp> fwd_lef_rank_buff(nlefs);
  std::vector<Bp> rev_lef_rank_buff(nlefs);

  // clang-format off
  std::vector<Lef> lefs{ // Lef{pos, lifetime, rank, nstalls_lef_lef, nstalls_lef_bar}
      Lef{{0, 0, 0}, {1, 0, 0}},
      Lef{{3, 0, 1}, {4, 0, 1}},
      Lef{{5, 0, 2}, {5, 0, 2}}
  };

  std::vector<ExtrusionBarrier> barriers{ // ExtrusionBarrier{pos, prob_of_block, motif_direction}
      ExtrusionBarrier{2, 0.99, '+'},
      ExtrusionBarrier{4, 0.99, '+'},
      ExtrusionBarrier{8, 0.99, '+'}
  };

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_CONVERSION
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
  modle::config c;
  c.bin_size = 2;
  const std::size_t nlefs = 3;

  std::vector<Bp> fwd_lef_rank_buff(nlefs);
  std::vector<Bp> rev_lef_rank_buff(nlefs);

  // clang-format off
  std::vector<Lef> lefs{ // Lef{pos, lifetime, rank, nstalls_lef_lef, nstalls_lef_bar}
      Lef{{0, 0, 0}, {1, 0, 0}},
      Lef{{3, 0, 1}, {4, 0, 1}},
      Lef{{5, 0, 2}, {8, 0, 2}}
  };

  std::vector<ExtrusionBarrier> barriers{ // ExtrusionBarrier{pos, prob_of_block, motif_direction}
      ExtrusionBarrier{2, 0.99, '+'},
      ExtrusionBarrier{4, 0.99, '+'},
      ExtrusionBarrier{8, 0.99, '-'}
  };

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_CONVERSION
  const std::vector<std::size_t> collisions_rev{-1UL, 0, 1};
  const std::vector<std::size_t> collisions_fwd{0, 1, 2};
  DISABLE_WARNING_POP
  assert(lefs.size() == nlefs); // NOLINT
  assert(collisions_fwd.size() == nlefs); // NOLINT
  assert(collisions_rev.size() == nlefs); // NOLINT
  // clang-format on

  apply_lef_bar_stalls_wrapper(c, lefs, barriers, collisions_rev, collisions_fwd);
}

}  // namespace modle::test::libmodle
