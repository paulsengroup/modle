// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/stats/tests.hpp"

#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>
#include <atomic>
#include <boost/process.hpp>
#include <catch2/catch.hpp>
#include <condition_variable>
#include <mutex>
#include <string_view>
#include <thread>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/common/random.hpp"
#include "modle/common/utils.hpp"

namespace modle::test::stats {
using namespace modle::stats;
using namespace std::string_view_literals;

struct FP {
  float n;
};

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

  CHECK(stats::binomial_test<TWO_SIDED>(k, n1) == result_two_sided_n1);
  CHECK(stats::binomial_test<TWO_SIDED>(k, n2) == result_two_sided_n2);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Binom test - less", "[math][short]") {
  const auto k = 100;
  const auto n = 143;

  const auto result_less = Approx(0.9999995628468261);

  CHECK(stats::binomial_test<LESS>(k, n) == result_less);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Binom test - greater", "[math][short]") {
  const auto k = 100;
  const auto n = 143;

  const auto result_greater = Approx(1.0502284150826767e-06);

  CHECK(stats::binomial_test<GREATER>(k, n) == result_greater);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Binom test - two-sided randomized", "[math][long]") {
  // random::PRNG_t rand_eng{4888025265521095494};//std::random_device{}()};
  random::PRNG_t rand_eng{std::random_device{}()};
  const usize iterations = 5'000;

  std::mutex data_mutex;
  std::condition_variable input_data_cv;
  std::condition_variable output_data_cv;
  std::atomic<bool> input_data_ready{false};
  std::atomic<bool> output_data_ready{false};

  i64 k{-1};
  i64 n{-1};

  std::atomic<double> pv_py{};

  std::thread t([&]() {
    run_scipy("two-sided"sv, k, n, pv_py, data_mutex, input_data_cv, output_data_cv,
              input_data_ready, output_data_ready);
  });

  for (usize i = 0; i < iterations; ++i) {  // NOLINT
    assert(!input_data_ready);              // NOLINT
    do {
      k = random::uniform_int_distribution<i64>{0, 5000}(rand_eng);      // NOLINT
      n = k + random::uniform_int_distribution<i64>{0, 5000}(rand_eng);  // NOLINT
    } while (k + n == 0);

    input_data_ready = true;
    input_data_cv.notify_one();

    {
      std::unique_lock<std::mutex> l(data_mutex);
      output_data_cv.wait(l, [&]() { return output_data_ready.load(); });
      output_data_ready = false;
    }
    const auto pv = stats::binomial_test<TWO_SIDED>(k, n);
    CHECK(Approx(pv).margin(1.0e-250) == pv_py);  // NOLINT
  }
  k = 0;
  n = 0;
  input_data_ready = true;
  input_data_cv.notify_one();
  t.join();
}

}  // namespace modle::test::stats
