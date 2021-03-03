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

[[nodiscard]] inline modle::config generate_mock_config() {
  modle::config c;
  c.number_of_lefs = 4;
  c.bin_size = 5;
  return c;
}

inline void check_lef_lef_collisions_wrapper(
    const modle::config& c, const std::vector<Lef>& lefs,
    const std::vector<std::size_t>& expected_collisions_rev,
    const std::vector<std::size_t>& expected_collisions_fwd) {
  const auto& nlefs = lefs.size();
  std::vector<Bp> fwd_rank_buff(nlefs);
  std::vector<Bp> rev_rank_buff(nlefs);

  std::vector<std::size_t> fwd_collision_counter_buff(nlefs, 0);
  std::vector<std::size_t> rev_collision_counter_buff(nlefs, 0);

  std::transform(lefs.begin(), lefs.end(), fwd_rank_buff.begin(),
                 [](const auto& lef) { return lef.fwd_unit.rank(); });
  std::transform(lefs.begin(), lefs.end(), rev_rank_buff.begin(),
                 [](const auto& lef) { return lef.rev_unit.rank(); });

  modle::Genome g(c, false);
  g.test_check_lef_lef_collisions(lefs, fwd_rank_buff, rev_rank_buff, fwd_collision_counter_buff,
                                  rev_collision_counter_buff);

  // fmt::print(stderr, "rev_collisions={}\n", absl::StrJoin(rev_collision_counter_buff, ", "));
  // fmt::print(stderr, "fwd_collisions={}\n", absl::StrJoin(fwd_collision_counter_buff, ", "));

  for (auto i = 0UL; i < lefs.size(); ++i) {
    CHECK(fwd_collision_counter_buff[i] == expected_collisions_fwd[i]);
    CHECK(rev_collision_counter_buff[i] == expected_collisions_rev[i]);
  }
}

TEST_CASE("Check LEF-LEF collision simple 001", "[simulation][short]") {
  modle::config c;
  c.number_of_lefs = 4;
  c.bin_size = 5;
  const auto& nlefs = c.number_of_lefs;

  std::vector<Bp> fwd_rank_buff(nlefs);
  std::vector<Bp> rev_rank_buff(nlefs);

  // clang-format off
  const std::vector<Lef> lefs{
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

TEST_CASE("Check LEF-LEF collision simple 002", "[simulation][short]") {
  modle::config c;
  c.number_of_lefs = 4;
  c.bin_size = 5;
  const auto& nlefs = c.number_of_lefs;

  std::vector<Bp> fwd_rank_buff(nlefs);
  std::vector<Bp> rev_rank_buff(nlefs);

  // clang-format off
  const std::vector<Lef> lefs{
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

}  // namespace modle::test::libmodle
