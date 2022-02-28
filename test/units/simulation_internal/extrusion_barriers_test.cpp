// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier

#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr

#include "modle/common/common.hpp"  // for dna::Direction

namespace modle::test::libmodle {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Extrusion barriers - comparisons", "[barriers][simulation][short]") {
  const auto b1 = ExtrusionBarrier{100, 0.0, 0.0, dna::fwd};
  const auto b2 = ExtrusionBarrier{500, 0.0, 0.0, dna::rev};

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
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Extrusion barriers - pblock", "[barriers][simulation][short]") {
  const auto pb1 = 0.85;
  const auto pnn1 = 0.7;  // not-occ -> not-occ
  const auto poo1 =       // occ -> occ
      ExtrusionBarrier::compute_blocking_to_blocking_transition_probabilities_from_pblock(pb1,
                                                                                          pnn1);
  const auto pb2 = 0.93;
  const auto pnn2 = 0.65;  // not-occ -> not-occ
  const auto poo2 =        // occ -> occ
      ExtrusionBarrier::compute_blocking_to_blocking_transition_probabilities_from_pblock(pb2,
                                                                                          pnn2);

  const auto b1 = ExtrusionBarrier{100, poo1, pnn1, dna::fwd};
  const auto b2 = ExtrusionBarrier{200, poo2, pnn2, dna::rev};

  CHECK(b1.occupancy() == Approx(pb1));
  CHECK(b2.occupancy() == Approx(pb2));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Extrusion barriers - blocking direction", "[barriers][simulation][short]") {
  const auto b1 = ExtrusionBarrier{100, 0.0, 0.0, dna::fwd};
  const auto b2 = ExtrusionBarrier{200, 0.0, 0.0, dna::rev};
  const auto b3 = ExtrusionBarrier{300, 0.0, 0.0, '+'};
  const auto b4 = ExtrusionBarrier{400, 0.0, 0.0, '-'};

  CHECK(b1.blocking_direction_major() == dna::rev);
  CHECK(b1.blocking_direction_minor() == dna::fwd);

  CHECK(b2.blocking_direction_major() == dna::fwd);
  CHECK(b2.blocking_direction_minor() == dna::rev);

  CHECK(b3.blocking_direction_major() == dna::rev);
  CHECK(b3.blocking_direction_minor() == dna::fwd);

  CHECK(b4.blocking_direction_major() == dna::fwd);
  CHECK(b4.blocking_direction_minor() == dna::rev);
}

}  // namespace modle::test::libmodle
