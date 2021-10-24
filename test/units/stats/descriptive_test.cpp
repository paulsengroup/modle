// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/stats/descriptive.hpp"

#include <algorithm>         // for transform
#include <catch2/catch.hpp>  // for Approx, operator==, AssertionHandler, operator""_catc...
#include <string_view>       // for string_view_literals
#include <vector>            // for vector

#include "modle/common/common.hpp"  // for u8, usize
#include "modle/common/utils.hpp"   // for identity::operator(), identity

namespace modle::test::stats {
using namespace modle::stats;
using namespace std::string_view_literals;

struct FP {
  float n;
};

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Mean", "[math][short]") {
  const std::vector<u8> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(5.0);

  CHECK(stats::mean(v1.begin(), v1.end()) == result);
  CHECK(stats::mean(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }) == result);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Moving average", "[math][short]") {
  const usize window_size = 3;
  const std::vector<u8> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const std::vector<double> results{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

  std::vector<double> output(results.size());
  REQUIRE(stats::moving_average(v1.begin(), v1.end(), output.begin(), window_size) ==
          v1.size() - window_size);

  for (usize i = 0; i < results.size(); ++i) {
    CHECK(results[i] == Approx(output[i]));
  }

  output.clear();
  output.resize(results.size());
  REQUIRE(stats::moving_average(v2.begin(), v2.end(), output.begin(), window_size,
                                [](const auto& fp) { return fp.n; }) == v1.size() - window_size);

  for (usize i = 0; i < results.size(); ++i) {
    CHECK(results[i] == Approx(output[i]));
  }

  output.clear();
  output.resize(1);
  REQUIRE(stats::moving_average(v1.begin(), v1.end(), output.begin(), v1.size() + 1) == 1);
  CHECK(stats::mean(v1.begin(), v1.end()) == Approx(output.front()));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Sum of squared deviations", "[math][short]") {
  const std::vector<u8> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(110.0);

  CHECK(stats::sum_of_squared_deviations(v1.begin(), v1.end()) == result);
  CHECK(stats::sum_of_squared_deviations(v2.begin(), v2.end(),
                                         [](const auto& fp) { return fp.n; }) == result);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Variance", "[math][short]") {
  const std::vector<u8> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(10.0);

  CHECK(stats::variance(v1.begin(), v1.end()) == result);
  CHECK(stats::variance(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }) == result);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Standard Deviation", "[math][short]") {
  const std::vector<u8> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(3.1622776601683795);

  CHECK(stats::standard_dev(v1.begin(), v1.end()) == result);
  CHECK(stats::standard_dev(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }) == result);
}
}  // namespace modle::test::stats
