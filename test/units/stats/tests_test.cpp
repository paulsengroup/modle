// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/stats/tests.hpp"

#include <fmt/compile.h>  // for format, FMT_COMPILE
#include <fmt/format.h>   // for format

#include <algorithm>                                              // for find_if
#include <atomic>                                                 // for atomic
#include <boost/exception/exception.hpp>                          // for clone_base
#include <boost/filesystem/path.hpp>                              // for path
#include <boost/fusion/sequence/intrinsic/at_key.hpp>             // for at_key
#include <boost/math/distributions/complement.hpp>                // for complement
#include <boost/math/distributions/detail/derived_accessors.hpp>  // for cdf, pdf
#include <boost/process/child.hpp>                                // for child::child
#include <boost/process/detail/child_decl.hpp>                    // for child
#include <boost/process/detail/posix/basic_pipe.hpp>              // for basic_pipe
#include <boost/process/io.hpp>                                   // for std_in, std_in_, std_out
#include <boost/process/pipe.hpp>                                 // for opstream, ipstream
#include <boost/process/search_path.hpp>                          // for search_path
#include <boost/random/uniform_int_distribution.hpp>              // for uniform_int_distribution
#include <cassert>                                                // for assert
#include <catch2/catch.hpp>                                       // for Approx, operator==, Ass...
#include <cmath>                                                  // for trunc
#include <condition_variable>                                     // for condition_variable
#include <memory>                                                 // for make_shared
#include <mutex>                                                  // for mutex, unique_lock
#include <ostream>                                                // for basic_ostream::flush
#include <random>                                                 // for random_device
#include <stdexcept>                                              // for overflow_error
#include <string>                                                 // for string, getline
#include <string_view>                                            // for string_view, basic_stri...
#include <thread>                                                 // for thread
#include <utility>                                                // for make_pair
#include <vector>                                                 // for vector

#include "modle/common/common.hpp"  // for i64, usize
#include "modle/common/random.hpp"  // for PRNG_t, uniform_int_dis...
#include "modle/common/utils.hpp"   // for parse_numeric_or_throw

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
TEST_CASE("Binom test - two-sided", "[stats][short]") {
  const auto k = 100;
  const auto n1 = 143;
  const auto n2 = 220;

  const auto result_two_sided_n1 = Approx(2.1004568301653535e-06);
  const auto result_two_sided_n2 = Approx(0.20009222902827126);

  CHECK(stats::binomial_test<TWO_SIDED>(k, n1) == result_two_sided_n1);
  CHECK(stats::binomial_test<TWO_SIDED>(k, n2) == result_two_sided_n2);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Binom test - less", "[stats][short]") {
  const auto k = 100;
  const auto n = 143;

  const auto result_less = Approx(0.9999995628468261);

  CHECK(stats::binomial_test<LESS>(k, n) == result_less);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Binom test - greater", "[stats][short]") {
  const auto k = 100;
  const auto n = 143;

  const auto result_greater = Approx(1.0502284150826767e-06);

  CHECK(stats::binomial_test<GREATER>(k, n) == result_greater);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Binom test - two-sided randomized", "[stats][long]") {
  random::PRNG_t rand_eng{4888025265521095494};
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
