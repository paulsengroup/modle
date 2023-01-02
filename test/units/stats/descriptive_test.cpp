// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/stats/descriptive.hpp"

#include <algorithm>  // for transform
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <string_view>  // for string_view_literals
#include <vector>       // for vector

#include "modle/common/common.hpp"  // for u8, usize
#include "modle/common/utils.hpp"   // for identity::operator(), identity

namespace modle::stats::test {

// See https://github.com/catchorg/Catch2/blob/v3.2.1/src/catch2/catch_approx.cpp#L27-L32
constexpr double DEFAULT_FP_TOLERANCE =
    static_cast<double>(std::numeric_limits<float>::epsilon() * 100);

struct FP {
  float n;
};

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Mean", "[stats][short]") {
  const std::vector<u8> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  CHECK_THAT(stats::mean(v1.begin(), v1.end()),
             Catch::Matchers::WithinRel(5.0, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(stats::mean(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }),
             Catch::Matchers::WithinRel(5.0, DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Moving average", "[stats][short]") {
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
    CHECK_THAT(results[i], Catch::Matchers::WithinRel(output[i], DEFAULT_FP_TOLERANCE));
  }

  output.clear();
  output.resize(results.size());
  REQUIRE(stats::moving_average(v2.begin(), v2.end(), output.begin(), window_size,
                                [](const auto& fp) { return fp.n; }) == v1.size() - window_size);

  for (usize i = 0; i < results.size(); ++i) {
    CHECK_THAT(results[i], Catch::Matchers::WithinRel(output[i], DEFAULT_FP_TOLERANCE));
  }

  output.clear();
  output.resize(1);
  REQUIRE(stats::moving_average(v1.begin(), v1.end(), output.begin(), v1.size() + 1) == 1);
  CHECK_THAT(stats::mean(v1.begin(), v1.end()),
             Catch::Matchers::WithinRel(output.front(), DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Sum of squared deviations", "[stats][short]") {
  const std::vector<u8> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  CHECK_THAT(stats::sum_of_squared_deviations(v1.begin(), v1.end()),
             Catch::Matchers::WithinRel(110.0, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(
      stats::sum_of_squared_deviations(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }),
      Catch::Matchers::WithinRel(110.0, DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Variance", "[stats][short]") {
  const std::vector<u8> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  CHECK_THAT(stats::variance(v1.begin(), v1.end()),
             Catch::Matchers::WithinRel(10.0, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(stats::variance(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }),
             Catch::Matchers::WithinRel(10.0, DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Standard Deviation", "[stats][short]") {
  const std::vector<u8> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = 3.1622776601683795;

  CHECK_THAT(stats::standard_dev(v1.begin(), v1.end()),
             Catch::Matchers::WithinRel(result, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(stats::standard_dev(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }),
             Catch::Matchers::WithinRel(result, DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("SED", "[stats][short]") {
  const std::vector<double> v1{-21.65393487, 43.10503862, -28.29768021, 34.10045857, -43.84182432,
                               -44.43248959, -5.82815809, 17.21738318,  20.59312136, 12.29917826};
  const std::vector<double> v2{33.80650183,  -23.34010343, -43.32704456, -37.97380374, 20.67554989,
                               -45.59723581, 13.65381298,  -45.2267169,  -28.00073264, 26.58993515};

  const std::vector<double> weights{0.97302005, 0.05226173, 0.15995629, 0.31495018, 0.95241483,
                                    0.87420081, 0.21360278, 0.0822476,  0.26402032, 0.49666325};
  // computed with scipy.spatial.distance.euclidean
  const auto result = 154.65977964758991;
  const auto result_weighted = 99.94033647108283;

  CHECK_THAT(stats::sed(v1.begin(), v1.end(), v2.begin()),
             Catch::Matchers::WithinRel(result, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(stats::weighted_sed(v1.begin(), v1.end(), v2.begin(), weights.begin()),
             Catch::Matchers::WithinRel(result_weighted, DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("RMSE", "[stats][short]") {
  const std::vector<double> v1{-21.65393487, 43.10503862, -28.29768021, 34.10045857, -43.84182432,
                               -44.43248959, -5.82815809, 17.21738318,  20.59312136, 12.29917826};
  const std::vector<double> v2{33.80650183,  -23.34010343, -43.32704456, -37.97380374, 20.67554989,
                               -45.59723581, 13.65381298,  -45.2267169,  -28.00073264, 26.58993515};

  const std::vector<double> weights{0.97302005, 0.05226173, 0.15995629, 0.31495018, 0.95241483,
                                    0.87420081, 0.21360278, 0.0822476,  0.26402032, 0.49666325};

  // computed with mean_squared_error from sklearn.metrics
  const auto result = 48.90771661061378;
  const auto result_weighted = 47.73515476405454;

  CHECK_THAT(stats::rmse(v1.begin(), v1.end(), v2.begin()),
             Catch::Matchers::WithinRel(result, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(stats::weighted_rmse(v1.begin(), v1.end(), v2.begin(), weights.begin()),
             Catch::Matchers::WithinRel(result_weighted, DEFAULT_FP_TOLERANCE));
}
}  // namespace modle::stats::test
