// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/stats/correlation.hpp"

#include <absl/container/flat_hash_set.h>  // for BitMask, operator!=
#include <absl/strings/match.h>            // for StartsWith
#include <fmt/format.h>                    // for format

#include <atomic>                         // for atomic
#include <boost/exception/exception.hpp>  // for clone_base
#include <cassert>                        // for assert
#include <catch2/catch.hpp>               // for Approx, operator==, AssertionHandler, operator...
#include <condition_variable>             // for condition_variable
#include <exception>                      // for exception
#include <mutex>                          // for mutex, unique_lock
#include <random>                         // for random_device
#include <stdexcept>                      // for overflow_error
#include <string>                         // for char_traits
#include <string_view>                    // for operator==, operator""sv, basic_string_view
#include <thread>                         // for thread
#include <utility>                        // for make_pair
#include <vector>                         // for vector

#include "./common.hpp"               // for generate_random_vect, external_corr, run_exter...
#include "modle/common/common.hpp"    // for u32, usize
#include "modle/common/random.hpp"    // for PRNG_t
#include "modle/common/smartdir.hpp"  // for SmartDir

namespace modle::test::stats {
using namespace modle::stats;

template <typename N>
static void test_correlation_w_random_vector(std::string_view method, usize vector_size,
                                             usize iterations, N min, N max) {
  static_assert(std::is_arithmetic_v<N>, "N should be an arithmetic type.");
  random::PRNG_t rand_eng{std::random_device{}()};

  std::mutex data_mutex;
  std::condition_variable input_data_cv;
  std::condition_variable output_data_cv;
  std::atomic<bool> input_data_ready{false};
  std::atomic<bool> output_data_ready{false};

  std::vector<N> v1(vector_size);
  std::vector<N> v2(vector_size);
  std::vector<double> weights(vector_size);

  std::atomic<double> cfx_reference{};
  std::atomic<double> pv_reference{};

  std::thread t([&]() {
    run_external_corr(method, v1, v2, weights, cfx_reference, pv_reference, data_mutex,
                      input_data_cv, output_data_cv, input_data_ready, output_data_ready);
  });

  auto pearson = Pearson<>{};
  auto spearman = Spearman<>{};
  for (usize i = 0; i < iterations; ++i) {
    assert(!input_data_ready);
    generate_random_vect(rand_eng, v1, min, max);
    generate_random_vect(rand_eng, v2, min, max);
    if (absl::StartsWith(method, "weighted_")) {
      generate_random_vect(rand_eng, weights, 0.0, 1.0);
      // fmt::print(stderr, FMT_STRING("{}\n"), fmt::join(weights, ","));
    }
    input_data_ready = true;
    input_data_cv.notify_one();
    {
      std::unique_lock<std::mutex> l(data_mutex);
      output_data_cv.wait(l, [&]() { return output_data_ready.load(); });
      output_data_ready = false;
    }

    const auto [cfx, pv] = [&]() {
      if (method == "pearson"sv) {
        const auto res = pearson(v1, v2);
        return std::make_pair(res.pcc, res.pvalue);
      }

      if (method == "weighted_pearson"sv) {
        const auto res = pearson(v1, v2, weights);
        return std::make_pair(res.pcc, res.pvalue);
      }

      if (method == "spearman"sv) {
        const auto res = spearman(v1, v2);
        return std::make_pair(res.rho, res.pvalue);
      }

      assert(method == "weighted_spearman"sv);
      const auto res = spearman(v1, v2, weights);
      return std::make_pair(res.rho, res.pvalue);
    }();

    CHECK(Approx(cfx) == cfx_reference);
    CHECK(Approx(pv) == pv_reference);
  }
  v1.clear();
  v2.clear();
  input_data_ready = true;
  input_data_cv.notify_one();
  t.join();
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Pearson wo ties", "[correlation][pearson][short]") {
  const std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};
  const std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto [pcc, pv] = Pearson<>{}(v1, v2);
  CHECK(Approx(pcc) == -0.033621194725622014);
  CHECK(Approx(pv) == 0.926536715854247);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Pearson wo ties", "[correlation][pearson][short]") {
  const std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};
  const std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const std::vector<double> w{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const auto [pcc, pv] = Pearson<>{}(v1, v2, w);
  CHECK(Approx(pcc) == 0.1892337717235999250409);
  CHECK(Approx(pv) == -1.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Pearson w ties", "[correlation][pearson][short]") {
  std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};
  std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto [pcc, pv] = Pearson<>{}(v1, v2);
  CHECK(Approx(pcc) == 0.16426413174421572);
  CHECK(Approx(pv) == 0.6502118872600098);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Pearson w ties", "[correlation][pearson][short]") {
  const std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};
  const std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const std::vector<double> w{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const auto [pcc, pv] = Pearson<>{}(v1, v2, w);
  CHECK(Approx(pcc) == 0.5009581087644285890548);
  CHECK(Approx(pv) == -1.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Pearson", "[correlation][pearson][short]") {
  random::PRNG_t rand_eng{2427588200550938527ULL};
  const auto v1 = generate_random_vect(rand_eng, 1'000, 0, 15'000);
  const auto v2 = generate_random_vect(rand_eng, 1'000, 0, 15'000);
  std::vector<double> dummy_weights{};
  const auto [pcc, pv] = Pearson<>{}(v1, v2);
  const auto [pcc_reference, pv_reference] = external_corr(v1, v2, dummy_weights, "pearson");
  CHECK(Approx(pcc) == pcc_reference);
  CHECK(Approx(pv) == pv_reference);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Pearson", "[correlation][pearson][short]") {
  random::PRNG_t rand_eng{17383284879759537016ULL};
  const auto v1 = generate_random_vect(rand_eng, 1'000, 0, 15'000);
  const auto v2 = generate_random_vect(rand_eng, 1'000, 0, 15'000);
  const auto w = generate_random_vect(rand_eng, 1'000, 0.0, 1.0);
  const auto [pcc, pv] = Pearson<>{}(v1, v2, w);
  const auto [pcc_reference, pv_reference] = external_corr(v1, v2, w, "weighted_pearson");
  CHECK(Approx(pcc) == pcc_reference);
  CHECK(Approx(pv) == pv_reference);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Pearson long", "[correlation][pearson][medium]") {
  test_correlation_w_random_vector("pearson", 1'000, 250, 0U, 15'000U);
  test_correlation_w_random_vector("pearson", 1'000, 250, -7'250, 7'250);
  test_correlation_w_random_vector("pearson", 1'000, 250, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Pearson long vect.", "[correlation][pearson][long]") {
  test_correlation_w_random_vector("pearson", 500'000, 2, 0U, 15'000U);
  test_correlation_w_random_vector("pearson", 500'000, 2, -7'250, 7'250);
  test_correlation_w_random_vector("pearson", 500'000, 2, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Pearson long", "[correlation][pearson][medium]") {
  test_correlation_w_random_vector("weighted_pearson", 1'000, 250, 0U, 15'000U);
  test_correlation_w_random_vector("weighted_pearson", 1'000, 250, -7'250, 7'250);
  test_correlation_w_random_vector("weighted_pearson", 1'000, 250, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Pearson long vect.", "[correlation][pearson][long]") {
  test_correlation_w_random_vector("weighted_pearson", 500'000, 2, 0U, 15'000U);
  test_correlation_w_random_vector("weighted_pearson", 500'000, 2, -7'250, 7'250);
  test_correlation_w_random_vector("weighted_pearson", 500'000, 2, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Spearman wo ties", "[correlation][spearman][short]") {
  const std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};
  const std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto [rho, pv] = Spearman<>{}(v1, v2);
  CHECK(Approx(rho) == -0.16363636363636364);
  CHECK(Approx(pv) == 0.6514773427962428);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Spearman w ties", "[correlation][spearman][short]") {
  const std::vector<u32> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};
  const std::vector<u32> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto [rho, pv] = Spearman<>{}(v1, v2);
  CHECK(Approx(rho) == 0.024316221747202587);
  CHECK(Approx(pv) == 0.9468397049085097);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Spearman", "[correlation][spearman][short]") {
  random::PRNG_t rand_eng{3860329809333667103ULL};
  const auto v1 = generate_random_vect(rand_eng, 1'000, 0, 15'000);
  const auto v2 = generate_random_vect(rand_eng, 1'000, 0, 15'000);
  std::vector<double> dummy_weights{};
  const auto [rho, pv] = Spearman<>{}(v1, v2);
  const auto [rho_reference, pv_reference] = external_corr(v1, v2, dummy_weights, "spearman");
  CHECK(Approx(rho) == rho_reference);
  CHECK(Approx(pv) == pv_reference);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Spearman", "[correlation][spearman][short]") {
  random::PRNG_t rand_eng{8469130800688738654ULL};
  const auto v1 = generate_random_vect(rand_eng, 1'000, 0, 15'000);
  const auto v2 = generate_random_vect(rand_eng, 1'000, 0, 15'000);
  const auto w = generate_random_vect(rand_eng, 1'000, 0.0, 1.0);
  const auto [rho, pv] = Spearman<>{}(v1, v2, w);
  const auto [rho_reference, pv_reference] = external_corr(v1, v2, w, "weighted_spearman");
  CHECK(Approx(rho) == rho_reference);
  CHECK(Approx(pv) == pv_reference);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Spearman long", "[correlation][spearman][long]") {
  test_correlation_w_random_vector("spearman", 1'000, 250, 0U, 15'000U);
  test_correlation_w_random_vector("spearman", 1'000, 250, -7'250, 7'250);
  test_correlation_w_random_vector("spearman", 1'000, 250, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Spearman long vect.", "[correlation][spearman][long]") {
  test_correlation_w_random_vector("spearman", 500'000, 2, 0U, 15'000U);
  test_correlation_w_random_vector("spearman", 500'000, 2, -7'250, 7'250);
  test_correlation_w_random_vector("spearman", 500'000, 2, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Spearman long", "[correlation][spearman][medium]") {
  test_correlation_w_random_vector("weighted_spearman", 1'000, 250, 0U, 15'000U);
  test_correlation_w_random_vector("weighted_spearman", 1'000, 250, -7'250, 7'250);
  test_correlation_w_random_vector("weighted_spearman", 1'000, 250, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Spearman long vect.", "[correlation][spearman][long]") {
  test_correlation_w_random_vector("weighted_spearman", 500'000, 2, 0U, 15'000U);
  test_correlation_w_random_vector("weighted_spearman", 500'000, 2, -7'250, 7'250);
  test_correlation_w_random_vector("weighted_spearman", 500'000, 2, -7'250.0, 7'250.0);
}

}  // namespace modle::test::stats
