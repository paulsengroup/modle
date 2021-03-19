#pragma once

#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr, SourceLineInfo
#include <cstdint>           // for uint32_t
#include <vector>            // for vector

#include "../correlation_utils.hpp"  // for compute_element_ranks, sort_range_by_idx
#include "../correlation_utils.hpp"

namespace modle::test::correlation {
using namespace modle::correlation;

TEST_CASE("Corr. test utils: sort vect. by idx", "[correlation][utils][short]") {
  const std::vector<uint32_t> v{10, 5, 67, 3, 60, 45, 49, 1000};
  const std::vector<uint32_t> expected{3, 1, 0, 5, 6, 4, 2, 7};
  const auto vi = modle::correlation::utils::sort_range_by_idx(v);
  REQUIRE(vi.size() == expected.size());
  for (auto i = 0UL; i < expected.size(); ++i) {
    CHECK(vi[i] == expected[i]);
  }
}

TEST_CASE("Corr. test utils: compute ranks wo ties", "[correlation][utils][short]") {
  const std::vector<uint32_t> v{10, 5, 67, 3, 60, 45, 49, 1000};
  const std::vector<uint32_t> expected{2, 1, 6, 0, 5, 3, 4, 7};
  const auto vi = modle::correlation::utils::compute_element_ranks(v);
  REQUIRE(vi.size() == expected.size());
  for (auto i = 0UL; i < expected.size(); ++i) {
    CHECK(vi[i] == expected[i]);
  }
}

TEST_CASE("Corr. test utils: compute ranks w ties", "[correlation][utils][short]") {
  const std::vector<uint32_t> v{10, 5, 67, 3, 67, 45, 49, 1000};
  const std::vector<double> expected{2, 1, 5.5, 0, 5.5, 3, 4, 7};
  const auto vi = modle::correlation::utils::compute_element_ranks(v);
  REQUIRE(vi.size() == expected.size());
  for (auto i = 0UL; i < expected.size(); ++i) {
    CHECK(vi[i] == expected[i]);
  }
}
}  // namespace modle::test::correlation
