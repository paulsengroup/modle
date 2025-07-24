// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/stats/correlation.hpp"

#include <cassert>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/common/numeric_utils.hpp"
#include "modle/common/string_utils.hpp"
#include "modle/compressed_io/compressed_io.hpp"

namespace modle::stats::test {

[[maybe_unused]] static const std::filesystem::path& data_dir() {
  static const std::filesystem::path data_dir{"test/data/unit_tests"};
  return data_dir;
}

// See https://github.com/catchorg/Catch2/blob/v3.2.1/src/catch2/catch_approx.cpp#L27-L32
constexpr double DEFAULT_FP_TOLERANCE =
    static_cast<double>(std::numeric_limits<float>::epsilon() * 100);

template <class N>
struct CorrelationBuffer {
  std::vector<N> v1{};
  std::vector<N> v2{};
  std::vector<double> weights{};

  double pcc{};
  double pcc_pval{};

  double rho{};
  double rho_pval{};

  double pccw{};
  // double pccw_pval{};

  double rhow{};
  // double rhow_pval{};

  explicit CorrelationBuffer(std::string_view buff) {
    // buff is expected to contain one line with the following syntax:
    // |v1|v2|weights|pcc|pcc_pval|rho|rho_pval|pccw|pccw_pval|rhow|rhow_pval|
    // v1, v2 and weights are comma-separated list of numbers (ints or floats)
    const std::vector<std::string_view> toks = str_split(buff, '|');
    REQUIRE(toks.size() == 3 + 8);

    for (const auto tok : str_split(toks[0], '\t')) {
      v1.push_back(utils::parse_numeric_or_throw<N>(tok));
    }

    v2.reserve(v1.size());
    for (const auto tok : str_split(toks[1], '\t')) {
      v2.push_back(utils::parse_numeric_or_throw<N>(tok));
    }

    weights.reserve(v1.size());
    for (const auto tok : str_split(toks[2], '\t')) {
      weights.push_back(utils::parse_numeric_or_throw<double>(tok));
    }

    REQUIRE(v1.size() == v2.size());
    REQUIRE(v1.size() == weights.size());

    pcc = utils::parse_numeric_or_throw<double>(toks[3]);
    pcc_pval = utils::parse_numeric_or_throw<double>(toks[4]);

    rho = utils::parse_numeric_or_throw<double>(toks[5]);
    rho_pval = utils::parse_numeric_or_throw<double>(toks[6]);

    pccw = utils::parse_numeric_or_throw<double>(toks[7]);
    // pccw_pval = utils::parse_numeric_or_throw<double>(toks[8]);

    rhow = utils::parse_numeric_or_throw<double>(toks[9]);
    // rhow_pval = utils::parse_numeric_or_throw<double>(toks[10]);
  }

  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
  void compute_and_compare_correlation(std::string_view method,
                                       double tolerance = DEFAULT_FP_TOLERANCE) const {
    const auto is_weighted = method.starts_with("weighted_");
    const auto is_pearson = method.ends_with("pearson");
    const auto is_spearman = method.ends_with("spearman");

    if (is_pearson && is_weighted) {
      const auto res = Pearson<>{}(v1, v2, weights);
      CHECK_THAT(res.pcc, Catch::Matchers::WithinRel(pccw, tolerance));
      CHECK(std::isnan(res.pvalue));
      return;
    }

    if (is_spearman && is_weighted) {
      const auto res = Spearman<>{}(v1, v2, weights);
      CHECK_THAT(res.rho, Catch::Matchers::WithinRel(rhow, tolerance));
      CHECK(std::isnan(res.pvalue));
      return;
    }

    if (is_pearson) {
      const auto res = Pearson<>{}(v1, v2);
      CHECK_THAT(res.pcc, Catch::Matchers::WithinRel(pcc, tolerance));
      CHECK_THAT(res.pvalue, Catch::Matchers::WithinRel(pcc_pval, tolerance));
      return;
    }

    assert(is_spearman);

    const auto res = Spearman<>{}(v1, v2);
    CHECK_THAT(res.rho, Catch::Matchers::WithinRel(rho, tolerance));
    CHECK_THAT(res.pvalue, Catch::Matchers::WithinRel(rho_pval, tolerance));
  }
};

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Pearson correlation wo/ ties", "[correlation][pearson][short]") {
  const std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};
  const std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto [pcc, pv] = Pearson<>{}(v1, v2);
  CHECK_THAT(pcc, Catch::Matchers::WithinRel(-0.033621194725622014, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(pv, Catch::Matchers::WithinRel(0.926536715854247, DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Weighted Pearson correlation wo/ ties", "[correlation][pearson][short]") {
  const std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};
  const std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const std::vector<double> w{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const auto [pcc, pv] = Pearson<>{}(v1, v2, w);
  CHECK_THAT(pcc, Catch::Matchers::WithinRel(0.1892337717235999250409, DEFAULT_FP_TOLERANCE));
  CHECK(std::isnan(pv));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Pearson correlation w/ ties", "[correlation][pearson][short]") {
  std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};
  std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto [pcc, pv] = Pearson<>{}(v1, v2);
  CHECK_THAT(pcc, Catch::Matchers::WithinRel(0.16426413174421572, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(pv, Catch::Matchers::WithinRel(0.6502118872600098, DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Weighted Pearson correlation w/ ties", "[correlation][pearson][short]") {
  const std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};
  const std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const std::vector<double> w{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const auto [pcc, pv] = Pearson<>{}(v1, v2, w);
  CHECK_THAT(pcc, Catch::Matchers::WithinRel(0.5009581087644285890548, DEFAULT_FP_TOLERANCE));
  CHECK(std::isnan(pv));
}

template <class N>
static void correlation_helper(const std::filesystem::path& path_to_reference_data,
                               std::string_view method) {
  static_assert(std::is_arithmetic_v<N>);

  compressed_io::Reader r(path_to_reference_data);
  for (std::string line; r.getline(line);) {
    CorrelationBuffer<N>(line).compute_and_compare_correlation(method);
  }
  REQUIRE(r.eof());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Pearson correlation", "[correlation][pearson][medium]") {
  const auto parent = data_dir() / "correlation";
  correlation_helper<u32>(parent / "correlation_uint.txt.xz", "pearson");
  correlation_helper<i32>(parent / "correlation_int.txt.xz", "pearson");
  correlation_helper<double>(parent / "correlation_float.txt.xz", "pearson");
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Pearson correlation (long vector)", "[correlation][pearson][medium]") {
  const auto parent = data_dir() / "correlation";
  correlation_helper<u32>(parent / "correlation_uint_long_vect.txt.xz", "pearson");
  correlation_helper<i32>(parent / "correlation_int_long_vect.txt.xz", "pearson");
  correlation_helper<double>(parent / "correlation_float_long_vect.txt.xz", "pearson");
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Weighted Pearson correlation", "[correlation][pearson][medium]") {
  const auto parent = data_dir() / "correlation";
  correlation_helper<u32>(parent / "correlation_uint.txt.xz", "weighted_pearson");
  correlation_helper<i32>(parent / "correlation_int.txt.xz", "weighted_pearson");
  correlation_helper<double>(parent / "correlation_float.txt.xz", "weighted_pearson");
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Weighted Pearson correlation (long vector)", "[correlation][pearson][long]") {
  const auto parent = data_dir() / "correlation";
  correlation_helper<u32>(parent / "correlation_uint_long_vect.txt.xz", "weighted_pearson");
  correlation_helper<i32>(parent / "correlation_int_long_vect.txt.xz", "weighted_pearson");
  correlation_helper<double>(parent / "correlation_float_long_vect.txt.xz", "weighted_pearson");
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Spearman correlation wo/ ties", "[correlation][spearman][short]") {
  const std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};
  const std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto [rho, pv] = Spearman<>{}(v1, v2);
  CHECK_THAT(rho, Catch::Matchers::WithinRel(-0.16363636363636364, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(pv, Catch::Matchers::WithinRel(0.6514773427962428, DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Spearman correlation w/ ties", "[correlation][spearman][short]") {
  const std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};
  const std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto [rho, pv] = Spearman<>{}(v1, v2);
  CHECK_THAT(rho, Catch::Matchers::WithinRel(0.024316221747202587, DEFAULT_FP_TOLERANCE));
  CHECK_THAT(pv, Catch::Matchers::WithinRel(0.9468397049085097, DEFAULT_FP_TOLERANCE));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Spearman correlation", "[correlation][spearman][long]") {
  const auto parent = data_dir() / "correlation";
  correlation_helper<u32>(parent / "correlation_uint.txt.xz", "spearman");
  correlation_helper<i32>(parent / "correlation_int.txt.xz", "spearman");
  correlation_helper<double>(parent / "correlation_float.txt.xz", "spearman");
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Spearman correlation (long vector)", "[correlation][spearman][long]") {
  const auto parent = data_dir() / "correlation";
  correlation_helper<u32>(parent / "correlation_uint_long_vect.txt.xz", "spearman");
  correlation_helper<i32>(parent / "correlation_int_long_vect.txt.xz", "spearman");
  correlation_helper<double>(parent / "correlation_float_long_vect.txt.xz", "spearman");
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Weighted Spearman correlation", "[correlation][spearman][medium]") {
  const auto parent = data_dir() / "correlation";
  correlation_helper<u32>(parent / "correlation_uint.txt.xz", "weighted_spearman");
  correlation_helper<i32>(parent / "correlation_int.txt.xz", "weighted_spearman");
  correlation_helper<double>(parent / "correlation_float.txt.xz", "weighted_spearman");
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Weighted Spearman correlation (long vector)", "[correlation][spearman][long]") {
  const auto parent = data_dir() / "correlation";
  correlation_helper<u32>(parent / "correlation_uint_long_vect.txt.xz", "weighted_spearman");
  correlation_helper<i32>(parent / "correlation_int_long_vect.txt.xz", "weighted_spearman");
  correlation_helper<double>(parent / "correlation_float_long_vect.txt.xz", "weighted_spearman");
}

}  // namespace modle::stats::test
