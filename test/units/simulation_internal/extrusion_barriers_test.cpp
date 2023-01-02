// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "modle/common/common.hpp"  // for dna::Direction

namespace modle::libmodle::test {

// See https://github.com/catchorg/Catch2/blob/v3.2.1/src/catch2/catch_approx.cpp#L27-L32
constexpr double DEFAULT_FP_TOLERANCE =
    static_cast<double>(std::numeric_limits<float>::epsilon() * 100);

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Extrusion barriers - comparisons", "[barriers][simulation][short]") {
  const auto b1 = ExtrusionBarrier{100, 0.0, 0.0, dna::FWD};
  const auto b2 = ExtrusionBarrier{500, 0.0, 0.0, dna::REV};

  CHECK(b1 == b1);
  CHECK(b1 != b2);

  CHECK(b1 < b2);
  CHECK(b2 > b1);

  CHECK(b1 <= b2);
  CHECK(b2 >= b1);
  CHECK(b1 <= b1);
  CHECK(b2 >= b2);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Extrusion barriers - occupancy", "[barriers][simulation][short]") {
  const auto occupancy1 = 0.85;
  const auto occupancy2 = 0.93;

  const auto stp_inactive1 = 0.7;
  const auto stp_inactive2 = 0.65;
  const auto stp_active1 =
      ExtrusionBarrier::compute_stp_active_from_occupancy(stp_inactive1, occupancy1);
  const auto stp_active2 =
      ExtrusionBarrier::compute_stp_active_from_occupancy(stp_inactive2, occupancy2);

  using State = ExtrusionBarriers::State;
  const auto barriers = ExtrusionBarriers{{100, 200},
                                          {dna::FWD, dna::REV},
                                          {stp_active1, stp_active2},
                                          {stp_inactive1, stp_inactive2},
                                          {State::ACTIVE, State::ACTIVE}};

  CHECK_THAT(barriers.occupancy(0), Catch::Matchers::WithinRel(occupancy1, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(barriers.occupancy(1), Catch::Matchers::WithinRel(occupancy2, DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Extrusion barriers - sort", "[barriers][simulation][short]") {
  using State = ExtrusionBarriers::State;
  auto barriers = ExtrusionBarriers{{100, 400, 200, 300},
                                    {dna::FWD, dna::FWD, dna::FWD, dna::FWD},
                                    {0.0, 0.0, 0.0, 0.0},
                                    {0.0, 0.0, 0.0, 0.0},
                                    {State::ACTIVE, State::ACTIVE, State::ACTIVE, State::ACTIVE}};

  barriers.sort();
  CHECK(std::is_sorted(barriers.pos().begin(), barriers.pos().end()));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Extrusion barriers - compute_occupancy", "[barriers][simulation][short]") {
  auto stp_active = 1.0;
  auto stp_inactive = 0.0;
  CHECK(ExtrusionBarrier::compute_occupancy_from_stp(stp_active, stp_inactive) == 1.0);

  stp_active = 0.0;
  stp_inactive = 1.0;
  CHECK(ExtrusionBarrier::compute_occupancy_from_stp(stp_active, stp_inactive) == 0.0);

  stp_active = 0.7;
  stp_inactive = 0.7;
  CHECK(ExtrusionBarrier::compute_occupancy_from_stp(stp_active, stp_inactive) == 0.5);

  stp_active = 0.7;
  stp_inactive = 0.5;
  CHECK_THAT(ExtrusionBarrier::compute_occupancy_from_stp(stp_active, stp_inactive),
             Catch::Matchers::WithinRel(0.625, DEFAULT_FP_TOLERANCE));

  stp_active = 1.0;
  stp_inactive = 0.5;
  CHECK(ExtrusionBarrier::compute_occupancy_from_stp(stp_active, stp_inactive) == 1.0);

  stp_active = 0.7;
  stp_inactive = 1.0;
  CHECK(ExtrusionBarrier::compute_occupancy_from_stp(stp_active, stp_inactive) == 0.0);
}

}  // namespace modle::libmodle::test
