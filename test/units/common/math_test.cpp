// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/common/math.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdint>
#include <vector>

namespace modle::test::math {
using namespace modle::math;

struct FP {
  float n;
};

TEST_CASE("Mean", "[utils][short]") {
  const std::vector<uint8_t> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(5.0);

  CHECK(math::mean(v1.begin(), v1.end()) == result);
  CHECK(math::mean(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }) == result);
}

TEST_CASE("Moving average", "[utils][short]") {
  const size_t window_size = 3;
  const std::vector<uint8_t> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const std::vector<double> results{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

  std::vector<double> output(results.size());
  REQUIRE(math::moving_average(v1.begin(), v1.end(), output.begin(), window_size) ==
          v1.size() - window_size);

  for (size_t i = 0; i < results.size(); ++i) {
    CHECK(results[i] == Approx(output[i]));
  }

  output.clear();
  output.resize(results.size());
  REQUIRE(math::moving_average(v2.begin(), v2.end(), output.begin(), window_size,
                               [](const auto& fp) { return fp.n; }) == v1.size() - window_size);

  for (size_t i = 0; i < results.size(); ++i) {
    CHECK(results[i] == Approx(output[i]));
  }

  output.clear();
  output.resize(1);
  REQUIRE(math::moving_average(v1.begin(), v1.end(), output.begin(), v1.size() + 1) == 1);
  CHECK(math::mean(v1.begin(), v1.end()) == Approx(output.front()));
}

TEST_CASE("Sum of squared deviations", "[utils][short]") {
  const std::vector<uint8_t> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(110.0);

  CHECK(math::sum_of_squared_deviations(v1.begin(), v1.end()) == result);
  CHECK(math::sum_of_squared_deviations(v2.begin(), v2.end(),
                                        [](const auto& fp) { return fp.n; }) == result);
}

TEST_CASE("Variance", "[utils][short]") {
  const std::vector<uint8_t> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(10.0);

  CHECK(math::variance(v1.begin(), v1.end()) == result);
  CHECK(math::variance(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }) == result);
}

TEST_CASE("Standard Deviation", "[utils][short]") {
  const std::vector<uint8_t> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(3.1622776601683795);

  CHECK(math::standard_dev(v1.begin(), v1.end()) == result);
  CHECK(math::standard_dev(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }) == result);
}
}  // namespace modle::test::math
