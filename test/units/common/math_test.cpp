// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/common/math.hpp"

#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>
#include <atomic>
#include <boost/process.hpp>
#include <catch2/catch.hpp>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <string_view>
#include <thread>
#include <vector>

#include "modle/common/random.hpp"
#include "modle/common/utils.hpp"

namespace modle::test::math {
using namespace modle::math;
using namespace std::string_view_literals;

struct FP {
  float n;
};

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Mean", "[math][short]") {
  const std::vector<uint8_t> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(5.0);

  CHECK(math::mean(v1.begin(), v1.end()) == result);
  CHECK(math::mean(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }) == result);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Moving average", "[math][short]") {
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Sum of squared deviations", "[math][short]") {
  const std::vector<uint8_t> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(110.0);

  CHECK(math::sum_of_squared_deviations(v1.begin(), v1.end()) == result);
  CHECK(math::sum_of_squared_deviations(v2.begin(), v2.end(),
                                        [](const auto& fp) { return fp.n; }) == result);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Variance", "[math][short]") {
  const std::vector<uint8_t> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(10.0);

  CHECK(math::variance(v1.begin(), v1.end()) == result);
  CHECK(math::variance(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }) == result);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Standard Deviation", "[math][short]") {
  const std::vector<uint8_t> v1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<FP> v2(v1.size());
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 [](const auto n) { return FP{static_cast<float>(n)}; });

  const auto result = Approx(3.1622776601683795);

  CHECK(math::standard_dev(v1.begin(), v1.end()) == result);
  CHECK(math::standard_dev(v2.begin(), v2.end(), [](const auto& fp) { return fp.n; }) == result);
}

// clang-format off
#define MAKE_BINOM_TEST_CASE(ALTERNATIVE) \
  std::make_pair(std::string_view{ALTERNATIVE},                           \
  "try:\n"                                                                \
  "  from scipy.stats import binomtest\n"                                 \
  "except ImportError:\n"                                                 \
  "  from scipy.stats import binom_test as binomtest\n"                   \
  "import fileinput\n"                                                    \
  "import sys\n"                                                          \
                                                                          \
  "for line in fileinput.input():\n"                                      \
  "  k, _, n = line.strip().partition(\"\\t\")\n"                         \
                                                                          \
  "  result = binomtest(int(k), int(n), alternative='" ALTERNATIVE "')\n" \
  "  if isinstance(result, float):\n"                                     \
  "    print(f\"{result:.32e}\", flush=True)\n"                           \
  "  else:\n"                                                             \
  "    print(f\"{result.pvalue:.32e}\", flush=True)\n"sv)


static constexpr utils::ConstMap<std::string_view, std::string_view, 3> python_cmds{
    MAKE_BINOM_TEST_CASE("two-sided"),
    MAKE_BINOM_TEST_CASE("less"),
    MAKE_BINOM_TEST_CASE("greater")
};  // clang-format on

template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
static void run_scipy(std::string_view method, N& n1, N& n2, std::atomic<double>& pv,
                      std::mutex& data_mutex, std::condition_variable& input_data_cv,
                      std::condition_variable& output_data_cv, std::atomic<bool>& input_data_ready,
                      std::atomic<bool>& output_data_ready) {
  boost::process::ipstream stdout_stream;
  boost::process::opstream stdin_stream;
  boost::process::child py(
      boost::process::search_path("python3").string(), "-c", std::string{python_cmds.at(method)},
      boost::process::std_in<stdin_stream, boost::process::std_out> stdout_stream);
  assert(py.running());  // NOLINT

  std::string sbuff;
  while (true) {
    {
      std::unique_lock<std::mutex> l(data_mutex);
      input_data_cv.wait(l, [&]() { return input_data_ready.load(); });
      input_data_ready = false;
    }

    if (n1 == N(0) && n2 == N(0)) {  // EOQ signal
      stdin_stream.pipe().close();
      py.wait();
      return;
    }

    sbuff = fmt::format(FMT_COMPILE("{}\t{}\n"), n1, n2);
    stdin_stream.write(sbuff.data(), static_cast<std::streamsize>(sbuff.size()));
    stdin_stream.flush();

    std::getline(stdout_stream, sbuff);
    pv = utils::parse_numeric_or_throw<double>(std::string_view(sbuff.data(), sbuff.size()));
    output_data_ready = true;
    output_data_cv.notify_one();
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Binom test - two-sided", "[math][short]") {
  const auto k = 100;
  const auto n1 = 143;
  const auto n2 = 220;

  const auto result_two_sided_n1 = Approx(2.1004568301653535e-06);
  const auto result_two_sided_n2 = Approx(0.20009222902827126);

  CHECK(math::binomial_test<TWO_SIDED>(k, n1) == result_two_sided_n1);
  CHECK(math::binomial_test<TWO_SIDED>(k, n2) == result_two_sided_n2);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Binom test - less", "[math][short]") {
  const auto k = 100;
  const auto n = 143;

  const auto result_less = Approx(0.9999995628468261);

  CHECK(math::binomial_test<LESS>(k, n) == result_less);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Binom test - greater", "[math][short]") {
  const auto k = 100;
  const auto n = 143;

  const auto result_greater = Approx(1.0502284150826767e-06);

  CHECK(math::binomial_test<GREATER>(k, n) == result_greater);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Binom test - two-sided randomized", "[math][long]") {
  // random::PRNG_t rand_eng{4888025265521095494};//std::random_device{}()};
  random::PRNG_t rand_eng{std::random_device{}()};
  const size_t iterations = 5'000;

  std::mutex data_mutex;
  std::condition_variable input_data_cv;
  std::condition_variable output_data_cv;
  std::atomic<bool> input_data_ready{false};
  std::atomic<bool> output_data_ready{false};

  int64_t k{-1};
  int64_t n{-1};

  std::atomic<double> pv_py{};

  std::thread t([&]() {
    run_scipy("two-sided"sv, k, n, pv_py, data_mutex, input_data_cv, output_data_cv,
              input_data_ready, output_data_ready);
  });

  for (size_t i = 0; i < iterations; ++i) {  // NOLINT
    assert(!input_data_ready);               // NOLINT
    do {
      k = random::uniform_int_distribution<int64_t>{0, 5000}(rand_eng);      // NOLINT
      n = k + random::uniform_int_distribution<int64_t>{0, 5000}(rand_eng);  // NOLINT
    } while (k + n == 0);

    input_data_ready = true;
    input_data_cv.notify_one();

    {
      std::unique_lock<std::mutex> l(data_mutex);
      output_data_cv.wait(l, [&]() { return output_data_ready.load(); });
      output_data_ready = false;
    }
    const auto pv = math::binomial_test(k, n);
    CHECK(Approx(pv).margin(1.0e-250) == pv_py);  // NOLINT
  }
  k = 0;
  n = 0;
  input_data_ready = true;
  input_data_cv.notify_one();
  t.join();
}

}  // namespace modle::test::math
