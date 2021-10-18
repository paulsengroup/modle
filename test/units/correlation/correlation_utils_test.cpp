// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT
#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr, SourceLineInfo
#include <vector>            // for vector

#include "modle/common/common.hpp"  // for u32
#include "modle/correlation.hpp"

namespace modle::test::correlation {
using namespace modle::correlation;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test utils: sort vect. by idx", "[correlation][utils][short]") {
  const std::vector<u32> v{10, 5, 67, 3, 60, 45, 49, 1000};
  const std::vector<u32> expected{3, 1, 0, 5, 6, 4, 2, 7};
  const auto vi = modle::correlation::utils::sort_range_by_idx(v);
  REQUIRE(vi.size() == expected.size());
  for (usize i = 0; i < expected.size(); ++i) {
    CHECK(vi[i] == expected[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test utils: compute ranks wo ties", "[correlation][utils][short]") {
  const std::vector<u32> v{10, 5, 67, 3, 60, 45, 49, 1000};
  const std::vector<u32> expected{2, 1, 6, 0, 5, 3, 4, 7};
  const auto vi = modle::correlation::utils::compute_element_ranks(v);
  REQUIRE(vi.size() == expected.size());
  for (usize i = 0; i < expected.size(); ++i) {
    CHECK(vi[i] == expected[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test utils: compute ranks w ties", "[correlation][utils][short]") {
  const std::vector<u32> v{10, 5, 67, 3, 67, 45, 49, 1000};
  const std::vector<double> expected{2, 1, 5.5, 0, 5.5, 3, 4, 7};
  const auto vi = modle::correlation::utils::compute_element_ranks(v);
  REQUIRE(vi.size() == expected.size());
  for (usize i = 0; i < expected.size(); ++i) {
    CHECK(vi[i] == expected[i]);
  }
}
}  // namespace modle::test::correlation
