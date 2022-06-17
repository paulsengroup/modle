// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/stats/correlation.hpp"

#include <absl/container/flat_hash_set.h>  // for BitMask, operator!=
#include <absl/strings/match.h>            // for StartsWith
#include <fmt/compile.h>                   // for FMT_COMPILE
#include <fmt/format.h>                    // for format

#include <boost/exception/exception.hpp>
#include <boost/process.hpp>
#include <cassert>           // for assert
#include <catch2/catch.hpp>  // for Approx, operator==, AssertionHandler, operator...
#include <exception>         // for exception
#include <ostream>           // for ostream
#include <random>            // for random_device
#include <string>            // for char_traits
#include <string_view>       // for operator==, operator""sv, basic_string_view
#include <utility>           // for make_pair
#include <vector>            // for vector

#include "modle/common/common.hpp"  // for u32, usize
#include "modle/common/const_map.hpp"
#include "modle/common/numeric_utils.hpp"
#include "modle/common/random.hpp"  // for PRNG_t

namespace modle::test::stats {
using namespace modle::stats;
using namespace std::string_view_literals;

// clang-format off
#define MAKE_CORR_TEST_CASE(METHOD)                           \
  std::make_pair(std::string_view{METHOD},                    \
     "#!/usr/bin/env python3\n"                               \
     "from scipy.stats import " METHOD "r\n"                  \
     "from numpy import fromstring\n"                         \
     "import fileinput\n"                                     \
                                                              \
     "for line in fileinput.input():\n"                       \
     "  v1, _, v2 = line.strip().partition(\"\\t\")\n"        \
     "  v1 = fromstring(v1, sep=\",\")\n"                     \
     "  v2 = fromstring(v2, sep=\",\")\n"                     \
                                                              \
     "  corr, pv = " METHOD "r(v1, v2)\n"                     \
     "  print(f\"{corr:.16e}\\t{pv:.16e}\", flush=True)\n"sv)

#define MAKE_WEIGHTED_CORR_TEST_CASE(METHOD)                       \
  std::make_pair(std::string_view{METHOD},                         \
     "#!/usr/bin/env Rscript\n"                                    \
     "library(\"wCorr\")\n"                                        \
     "f <- file(\"stdin\")\n"                                      \
     "open(f)\n"                                                   \
     "method <- sub(\"weighted_\", \"\", \"" METHOD "\")\n"        \
     "while(length(line <- readLines(f, n=1)) > 0) {\n"            \
     "   toks <- strsplit(line, '\\t')[[1]]\n"                     \
     "   v1 <- as.numeric(strsplit(toks[[1]], ',')[[1]])\n"        \
     "   v2 <- as.numeric(strsplit(toks[[2]], ',')[[1]])\n"        \
     "   w <- as.numeric(strsplit(toks[[3]], ',')[[1]])\n"         \
                                                                   \
     "   corr <- weightedCorr(v1, v2, weights=w, method=method)\n" \
     "   write(sprintf(\"%.22f\\t-1\", corr), stdout())\n"         \
     "   flush.console()\n"                                        \
     "}"sv)

static constexpr utils::ConstMap<std::string_view, std::string_view, 4> external_cmds{
    MAKE_CORR_TEST_CASE("pearson"),
    MAKE_CORR_TEST_CASE("spearman"),
    MAKE_WEIGHTED_CORR_TEST_CASE("weighted_pearson"),
    MAKE_WEIGHTED_CORR_TEST_CASE("weighted_spearman")
};  // clang-format on

template <class N>
void generate_random_vector(random::PRNG_t& rand_eng, std::vector<N>& buff, N min, N max,
                            bool allow_duplicates = true) {
  using DistrT =
      typename std::conditional<std::is_floating_point_v<N>, random::uniform_real_distribution<N>,
                                random::uniform_int_distribution<N>>::type;
  DistrT dist(min, max);
  if (allow_duplicates) {
    std::generate(buff.begin(), buff.end(), [&]() { return dist(rand_eng); });
    return;
  }

  absl::flat_hash_set<u32> s(buff.size());
  while (s.size() < buff.size()) {
    s.insert(dist(rand_eng));
  }
  std::copy(s.begin(), s.end(), buff.begin());
}

template <class N>
[[nodiscard]] std::vector<N> generate_random_vector(random::PRNG_t& rand_eng, usize size, N min,
                                                    N max, bool allow_duplicates = true) {
  std::vector<N> v(size);
  generate_random_vector(rand_eng, v, min, max, allow_duplicates);
  return v;
}

template <class N1, class N2>
class CorrelationBuff {
  static_assert(std::is_arithmetic_v<N1>);
  static_assert(std::is_arithmetic_v<N2>);
  std::vector<N1> _v1;
  std::vector<N1> _v2;
  std::vector<N2> _weights;

 public:
  CorrelationBuff() = delete;
  CorrelationBuff(usize size, bool with_weights)
      : _v1(size, N1(0)), _v2(size, N1(0)), _weights(with_weights ? size : 0, N2(0)) {}
  [[nodiscard]] const std::vector<N1>& v1() const noexcept { return this->_v1; }
  [[nodiscard]] const std::vector<N1>& v2() const noexcept { return this->_v2; }
  [[nodiscard]] const std::vector<N2>& weights() const noexcept { return this->_weights; }

  void generate(random::PRNG_t& rand_eng, N1 min, N1 max) {
    assert(this->_v1.size() == this->_v2.size());
    generate_random_vector(rand_eng, this->_v1, min, max);
    generate_random_vector(rand_eng, this->_v2, min, max);
    assert(this->_v1.size() == this->_v2.size());

    if (!this->_weights.empty()) {
      assert(this->_v1.size() == this->_weights.size());
      generate_random_vector(rand_eng, this->_weights, N2(0), N2(1));
      assert(this->_v1.size() == this->_weights.size());
    }
  }

  void serialize(std::string& buff) const {
    if (!this->_weights.empty()) {
      buff = fmt::format(FMT_COMPILE("{}\t{}\t{}\n"), fmt::join(this->_v1, ","),
                         fmt::join(this->_v2, ","), fmt::join(this->_weights, ","));
      return;
    }

    buff =
        fmt::format(FMT_COMPILE("{}\t{}\n"), fmt::join(this->_v1, ","), fmt::join(this->_v2, ","));
  }
  [[nodiscard]] std::string serialize() const {
    std::string buff;
    this->serialize(buff);
    return buff;
  }
};

template <class N1, class N2>
[[nodiscard]] static std::pair<double, double> compute_correlation(
    std::string_view method, const CorrelationBuff<N1, N2>& data) {
  const auto is_weighted = absl::StartsWith(method, "weighted_");
  const auto is_pearson = absl::EndsWith(method, "pearson");
  const auto is_spearman = absl::EndsWith(method, "spearman");

  if (is_pearson && is_weighted) {
    const auto res = Pearson<>{}(data.v1(), data.v2(), data.weights());
    return std::make_pair(res.pcc, res.pvalue);
  }

  if (is_spearman && is_weighted) {
    const auto res = Spearman<>{}(data.v1(), data.v2(), data.weights());
    return std::make_pair(res.rho, res.pvalue);
  }

  if (is_pearson) {
    const auto res = Pearson<>{}(data.v1(), data.v2());
    return std::make_pair(res.pcc, res.pvalue);
  }

  assert(is_spearman);

  const auto res = Spearman<>{}(data.v1(), data.v2());
  return std::make_pair(res.rho, res.pvalue);
}

class ExternalCorrelationRunner {
  boost::process::ipstream _stdout{};
  boost::process::ipstream _stderr{};
  boost::process::opstream _stdin{};

  boost::process::child _c{};

 public:
  ExternalCorrelationRunner() = delete;
  explicit ExternalCorrelationRunner(std::string_view method) {
    if (absl::StartsWith(method, "weighted_")) {
      this->_c = boost::process::child(
          boost::process::search_path("Rscript").string(), "--quiet", "-e",
          std::string{external_cmds.at(method)},
          boost::process::std_in<this->_stdin, boost::process::std_out> this->_stdout,
          boost::process::std_err > this->_stderr);
    } else {
      this->_c = boost::process::child(
          boost::process::search_path("python3").string(), "-c",
          std::string{external_cmds.at(method)},
          boost::process::std_in<this->_stdin, boost::process::std_out> this->_stdout,
          boost::process::std_err > this->_stderr);
    }
    assert(this->_c.running());
  }

  ~ExternalCorrelationRunner() noexcept {
    try {
      assert(!this->_c.running());
      if (!this->_c) {
        std::ostringstream buff;
        buff << this->_stderr.rdbuf();
        if (const auto s = buff.str(); !s.empty()) {
          throw std::runtime_error(fmt::format(FMT_STRING("{}\n"), s));
        }
      }
    } catch (const std::exception& e) {
      fmt::print(
          stderr,
          FMT_STRING("An exception was raised while calling ~ExternalCorrelationRunner():  {}"),
          e.what());
    } catch (...) {
      fmt::print(
          stderr,
          FMT_STRING("An unknown error occurred while calling ~ExternalCorrelationRunner()\n"));
    }
  }

  void wait(std::chrono::milliseconds duration = std::chrono::seconds(15)) {
    this->_c.wait_for(duration);
  }

  void write_to_stdin(const std::string& msg) {
    assert(!msg.empty() && msg.back() == '\n');
    this->_stdin.write(msg.data(), static_cast<std::streamsize>(msg.size()));
    this->_stdin.flush();
  }

  void read_from_stdout(std::string& buff) { std::getline(this->_stdout, buff); }
  [[nodiscard]] std::string read_from_stdout() {
    std::string buff;
    this->read_from_stdout(buff);
    return buff;
  }

  void read_from_stderr(std::string& buff) { std::getline(this->_stderr, buff); }
  [[nodiscard]] std::string read_from_stderr() {
    std::string buff;
    this->read_from_stderr(buff);
    return buff;
  }

  void signal_end_of_input() { this->_stdin.pipe().close(); }
};

template <class N1, class N2 = double, usize iterations = 1>
static void run_correlation_test(std::string_view method, usize vector_size, N1 min, N1 max,
                                 u64 seed = std::random_device{}()) {
  static_assert(std::is_arithmetic_v<N1>);
  static_assert(std::is_arithmetic_v<N2>);

  const auto weighted_corr = absl::StartsWith(method, "weighted_");
  CorrelationBuff<N1, N2> buff(vector_size, weighted_corr);
  std::string sbuff;
  std::vector<std::string_view> tok_buff;

  auto rand_eng = random::PRNG(seed);

  ExternalCorrelationRunner runner(method);
  try {
    for (usize i = 0; i < iterations; ++i) {
      buff.generate(rand_eng, min, max);
      buff.serialize(sbuff);

      runner.write_to_stdin(sbuff);
      const auto [corr, pv] = compute_correlation(method, buff);
      runner.read_from_stdout(sbuff);
      tok_buff = absl::StrSplit(sbuff, '\t');
      REQUIRE(tok_buff.size() == 2);

      const auto expected_corr = utils::parse_numeric_or_throw<double>(tok_buff.front());
      const auto expected_pv = utils::parse_numeric_or_throw<double>(tok_buff.back());

      CHECK(Approx(expected_corr) == corr);
      CHECK(Approx(expected_pv) == pv);
    }
    runner.signal_end_of_input();
    runner.wait();
  } catch (const std::exception& e) {
    const auto stderr_ = runner.read_from_stderr();
    if (!stderr_.empty()) {
      throw std::runtime_error(fmt::format(FMT_STRING("{}:\n{}"), e.what(), stderr_));
    }
    throw;
  }
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
TEST_CASE("Corr. test: Pearson (SciPy)", "[correlation][pearson][short]") {
  run_correlation_test("pearson", 1'000, 0, 15'000, 2427588200550938527ULL);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Pearson (wCorr)", "[correlation][pearson][short]") {
  run_correlation_test("weighted_pearson", 1'000, 0, 15'000, 17383284879759537016ULL);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Pearson long (SciPy)", "[correlation][pearson][medium]") {
  run_correlation_test<u32, double, 250>("pearson", 1'000, 0, 15'000);
  run_correlation_test<i32, double, 250>("pearson", 1'000, -7'250, 7'250);
  run_correlation_test<double, double, 250>("pearson", 1'000, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Pearson long vect. (SciPy)", "[correlation][pearson][long]") {
  run_correlation_test<u32, double, 5>("pearson", 100'000, 0U, 15'000U);
  run_correlation_test<i32, double, 5>("pearson", 100'000, -7'250, 7'250);
  run_correlation_test<double, double, 5>("pearson", 100'000, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Pearson long (wCorr)", "[correlation][pearson][medium]") {
  run_correlation_test<u32, double, 250>("weighted_pearson", 1'000, 0U, 15'000U);
  run_correlation_test<i32, double, 250>("weighted_pearson", 1'000, -7'250, 7'250);
  run_correlation_test<double, double, 250>("weighted_pearson", 1'000, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Pearson long vect. (wCorr)", "[correlation][pearson][long]") {
  run_correlation_test<u32, double, 5>("weighted_pearson", 100'000, 0U, 15'000U);
  run_correlation_test<i32, double, 5>("weighted_pearson", 100'000, -7'250, 7'250);
  run_correlation_test<double, double, 5>("weighted_pearson", 100'000, -7'250.0, 7'250.0);
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
TEST_CASE("Corr. test: Spearman (SciPy)", "[correlation][spearman][short]") {
  run_correlation_test("spearman", 1'000, 0, 15'000, 3860329809333667103ULL);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Spearman (wCorr)", "[correlation][spearman][short]") {
  run_correlation_test("weighted_spearman", 1'000, 0, 15'000, 8469130800688738654ULL);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Spearman long (SciPy)", "[correlation][spearman][long]") {
  run_correlation_test<u32, double, 250>("spearman", 1'000, 0, 15'000);
  run_correlation_test<i32, double, 250>("spearman", 1'000, -7'250, 7'250);
  run_correlation_test<double, double, 250>("spearman", 1'000, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Spearman long vect. (SciPy)", "[correlation][spearman][long]") {
  run_correlation_test<u32, double, 5>("spearman", 100'000, 0U, 15'000U);
  run_correlation_test<i32, double, 5>("spearman", 100'000, -7'250, 7'250);
  run_correlation_test<double, double, 5>("spearman", 100'000, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Spearman long (wCorr)", "[correlation][spearman][medium]") {
  run_correlation_test<u32, double, 250>("weighted_spearman", 1'000, 0U, 15'000U);
  run_correlation_test<i32, double, 250>("weighted_spearman", 1'000, -7'250, 7'250);
  run_correlation_test<double, double, 250>("weighted_spearman", 1'000, -7'250.0, 7'250.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Corr. test: Weighted Spearman long vect. (wCorr)", "[correlation][spearman][long]") {
  run_correlation_test<u32, double, 5>("weighted_spearman", 100'000, 0U, 15'000U);
  run_correlation_test<i32, double, 5>("weighted_spearman", 100'000, -7'250, 7'250);
  run_correlation_test<double, double, 5>("weighted_spearman", 100'000, -7'250.0, 7'250.0);
}

}  // namespace modle::test::stats
