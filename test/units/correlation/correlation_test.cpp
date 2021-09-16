// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/correlation.hpp"  // for compute_pearson, compute_pearson_significance, compute_...

#include <absl/container/flat_hash_set.h>  // for BitMask, operator!=

#include <atomic>              // for atomic
#include <cassert>             // for assert
#include <catch2/catch.hpp>    // for Approx, operator==, AssertionHandler, operator...
#include <condition_variable>  // for condition_variable
#include <cstdint>             // for uint32_t
#include <mutex>               // for mutex, unique_lock
#include <random>              // for random_device
#include <string>              // for char_traits
#include <string_view>         // for operator==, basic_string_view, string_view
#include <thread>              // for thread
#include <vector>              // for vector

#include "./common.hpp"               // for generate_random_vect, corr_scipy, run_scipy_corr
#include "modle/common/random.hpp"    // for PRNG_t
#include "modle/common/smartdir.hpp"  // for SmartDir

namespace modle::test {
const auto cleanup_on_exit{true};         // Useful for debugging
const SmartDir testdir{cleanup_on_exit};  // NOLINT Using auto here upsets GCC8
}  // namespace modle::test

namespace modle::test::correlation {
using namespace modle::correlation;

template <typename N>
static void test_correlation_w_random_vector(std::string_view method, size_t vector_size,
                                             size_t iterations, N min, N max) {
  static_assert(std::is_arithmetic_v<N>, "N should be an arithmetic type.");
  random::PRNG_t rand_eng{std::random_device{}()};

  std::mutex data_mutex;
  std::condition_variable input_data_cv;
  std::condition_variable output_data_cv;
  std::atomic<bool> input_data_ready{false};
  std::atomic<bool> output_data_ready{false};

  std::vector<N> v1(vector_size);
  std::vector<N> v2(vector_size);

  std::atomic<double> rho_py{};
  std::atomic<double> pv_py{};

  std::thread t([&]() {
    run_scipy_corr(method, v1, v2, rho_py, pv_py, data_mutex, input_data_cv, output_data_cv,
                   input_data_ready, output_data_ready);
  });

  for (size_t i = 0; i < iterations; ++i) {        // NOLINT
    assert(!input_data_ready);                     // NOLINT
    generate_random_vect(rand_eng, v1, min, max);  // NOLINT
    generate_random_vect(rand_eng, v2, min, max);  // NOLINT
    input_data_ready = true;
    input_data_cv.notify_one();
    {
      std::unique_lock<std::mutex> l(data_mutex);
      output_data_cv.wait(l, [&]() { return output_data_ready.load(); });
      output_data_ready = false;
    }

    const auto [rho, pv] = [&]() {
      double rho_;
      double pv_;
      if (method == "pearson") {
        rho_ = compute_pearson(v1, v2);
        pv_ = compute_pearson_significance(rho_, v1.size());
      } else if (method == "spearman") {
        rho_ = compute_spearman(v1, v2);
        pv_ = compute_spearman_significance(rho_, v1.size());
      }
      /* else if (method == "kendall") {
        rho = compute_kendall(v1, v2);
        pv_ = compute_kendall_significance(rho, v1.size());
      }
      */
      return std::make_pair(rho_, pv_);
    }();

    CHECK(Approx(rho) == rho_py);
    CHECK(Approx(pv) == pv_py);
  }
  v1.clear();
  v2.clear();
  input_data_ready = true;
  input_data_cv.notify_one();
  t.join();
}

TEST_CASE("Corr. test: Pearson wo ties", "[correlation][pearson][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  const auto pcc = compute_pearson(v1, v2);
  const auto pv = compute_pearson_significance(pcc, v1.size());
  CHECK(Approx(pcc) == -0.033621194725622014);  // NOLINT
  CHECK(Approx(pv) == 0.926536715854247);       // NOLINT
}

TEST_CASE("Corr. test: Pearson w ties", "[correlation][pearson][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  const auto pcc = compute_pearson(v1, v2);
  const auto pv = compute_pearson_significance(pcc, v1.size());
  CHECK(Approx(pcc) == 0.16426413174421572);  // NOLINT
  CHECK(Approx(pv) == 0.6502118872600098);    // NOLINT
}

TEST_CASE("Corr. test: Pearson Scipy", "[correlation][pearson][short]") {
  random::PRNG_t rand_eng{std::random_device{}()};
  auto v1 = generate_random_vect(rand_eng, 1'000, 0, 15'000);  // NOLINT
  auto v2 = generate_random_vect(rand_eng, 1'000, 0, 15'000);  // NOLINT
  const auto pcc = compute_pearson(v1, v2);
  const auto pv = compute_pearson_significance(pcc, v1.size());
  const auto [pcc_py, pv_py] = corr_scipy(v1, v2, "pearson");
  CHECK(Approx(pcc) == pcc_py);
  CHECK(Approx(pv) == pv_py);
}

TEST_CASE("Corr. test: Pearson Scipy long", "[correlation][pearson][medium]") {
  test_correlation_w_random_vector("pearson", 1'000, 250, 0U, 15'000U);
  test_correlation_w_random_vector("pearson", 1'000, 250, -7'250, 7'250);
  test_correlation_w_random_vector("pearson", 1'000, 250, -7'250.0, 7'250.0);
}

TEST_CASE("Corr. test: Pearson Scipy long vect.", "[correlation][pearson][long]") {
  test_correlation_w_random_vector("pearson", 500'000, 2, 0U, 15'000U);
  test_correlation_w_random_vector("pearson", 500'000, 2, -7'250, 7'250);
  test_correlation_w_random_vector("pearson", 500'000, 2, -7'250.0, 7'250.0);
}

TEST_CASE("Corr. test: Spearman wo ties", "[correlation][spearman][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  const auto rho = compute_spearman(v1, v2);
  const auto pv = compute_spearman_significance(rho, v1.size());
  CHECK(Approx(rho) == -0.16363636363636364);  // NOLINT
  CHECK(Approx(pv) == 0.6514773427962428);     // NOLINT
}

TEST_CASE("Corr. test: Spearman w ties", "[correlation][spearman][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  const auto rho = compute_spearman(v1, v2);
  const auto pv = compute_spearman_significance(rho, v1.size());
  CHECK(Approx(rho) == 0.024316221747202587);  // NOLINT
  CHECK(Approx(pv) == 0.9468397049085097);     // NOLINT
}

TEST_CASE("Corr. test: Spearman Scipy", "[correlation][spearman][short]") {
  random::PRNG_t rand_eng{std::random_device{}()};
  auto v1 = generate_random_vect(rand_eng, 1'000, 0, 15'000);  // NOLINT
  auto v2 = generate_random_vect(rand_eng, 1'000, 0, 15'000);  // NOLINT
  const auto rho = compute_spearman(v1, v2);
  const auto pv = compute_spearman_significance(rho, v1.size());
  const auto [rho_py, pv_py] = corr_scipy(v1, v2, "spearman");
  CHECK(Approx(rho) == rho_py);
  CHECK(Approx(pv) == pv_py);
}

TEST_CASE("Corr. test: Spearman Scipy long", "[correlation][spearman][long]") {
  test_correlation_w_random_vector("spearman", 1'000, 250, 0U, 15'000U);
  test_correlation_w_random_vector("spearman", 1'000, 250, -7'250, 7'250);
  test_correlation_w_random_vector("spearman", 1'000, 250, -7'250.0, 7'250.0);
}

TEST_CASE("Corr. test: Spearman Scipy long vect.", "[correlation][spearman][long]") {
  test_correlation_w_random_vector("spearman", 500'000, 2, 0U, 15'000U);
  test_correlation_w_random_vector("spearman", 500'000, 2, -7'250, 7'250);
  test_correlation_w_random_vector("spearman", 500'000, 2, -7'250.0, 7'250.0);
}

TEST_CASE("SED", "[correlation][sed][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  const auto sed = compute_sed(v1, v2);
  CHECK(Approx(sed) == 13125.999999999998);  // NOLINT
}

}  // namespace modle::test::correlation
