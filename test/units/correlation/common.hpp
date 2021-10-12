// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once
#include <absl/container/flat_hash_set.h>  // for flat_hash_set, operator!=
#include <absl/hash/hash.h>                // for Hash
#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>  // for max, generate
#include <array>      // for array
#include <boost/container_hash/hash.hpp>
#include <boost/filesystem/path.hpp>  // for create_directories, exists, path, remove
#include <boost/process.hpp>
#include <cstdint>    // for uint32_t, int32_t, int64_t
#include <cstdio>     // for pclose, fgets, popen, FILE
#include <fstream>    // for basic_ofstream, operator<<, basic_ostream
#include <memory>     // for allocator, unique_ptr
#include <numeric>    // for iota
#include <stdexcept>  // for runtime_error
#include <string>     // for string, operator+, char_traits, stod
#include <string_view>
#include <type_traits>  // for enable_if, is_arithmetic
#include <utility>      // for pair
#include <vector>       // for vector

#include "modle/common/random.hpp"
#include "modle/common/utils.hpp"

namespace modle::test::correlation {
using namespace std::string_view_literals;

// clang-format off
#define MAKE_CORR_TEST_CASE(METHOD)                           \
  std::make_pair(std::string_view{METHOD},                    \
     "from scipy.stats import " METHOD "\n"                   \
     "from numpy import fromstring\n"                         \
     "import fileinput\n"                                     \
                                                              \
     "for line in fileinput.input():\n"                       \
     "  v1, _, v2 = line.strip().partition(\"\\t\")\n"        \
     "  v1 = fromstring(v1, sep=\",\")\n"                     \
     "  v2 = fromstring(v2, sep=\",\")\n"                     \
                                                              \
     "  corr, pv = " METHOD "(v1, v2)\n"                      \
     "  print(f\"{corr:.16e}\\t{pv:.16e}\", flush=True)\n"sv)

static constexpr utils::ConstMap<std::string_view, std::string_view, 3> python_cmds{
    MAKE_CORR_TEST_CASE("pearsonr"),
    MAKE_CORR_TEST_CASE("spearmanr"),
    MAKE_CORR_TEST_CASE("kendalltau")
};  // clang-format on

template <class N, class = std::enable_if_t<std::is_arithmetic_v<N>>>
inline void generate_random_vect(random::PRNG_t& rand_eng, std::vector<N>& buff, N min, N max,
                                 bool allow_duplicates = true) {
  using dist_t =
      typename std::conditional<std::is_floating_point_v<N>, random::uniform_real_distribution<N>,
                                random::uniform_int_distribution<N>>::type;
  dist_t dist(min, max);
  if (allow_duplicates) {
    std::generate(buff.begin(), buff.end(), [&]() { return dist(rand_eng); });
  } else {
    absl::flat_hash_set<uint32_t> s;
    while (s.size() < buff.size()) {
      s.insert(dist(rand_eng));
    }
    std::copy(s.begin(), s.end(), buff.begin());
  }
}

template <class N = uint32_t, class = std::enable_if_t<std::is_arithmetic_v<N>>>
[[nodiscard]] inline std::vector<N> generate_random_vect(random::PRNG_t& rand_eng, size_t size,
                                                         uint32_t min, uint32_t max,
                                                         bool allow_duplicates = true) {
  std::vector<N> v(size);
  generate_random_vect(rand_eng, v, min, max, allow_duplicates);
  return v;
}

inline std::pair<std::vector<uint32_t>, std::vector<uint32_t>> generate_correlated_vects(
    random::PRNG_t& rand_eng, uint32_t size) {
  random::uniform_int_distribution<int32_t> dist(static_cast<int32_t>(size) / -50,  // NOLINT
                                                 static_cast<int32_t>(size / 50));  // NOLINT
  std::vector<uint32_t> v1(size);
  std::vector<uint32_t> v2(size);
  std::iota(v1.begin(), v1.end(), 0);
  std::iota(v2.begin(), v2.end(), 0);
  for (size_t i = 0; i < size; ++i) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_USELESS_CAST
    auto n = static_cast<int64_t>(v1[i]) + dist(rand_eng);
    v1[i] = static_cast<uint32_t>(std::max(int64_t(0), n));
    n = static_cast<int64_t>(v2[i]) + dist(rand_eng);
    v2[i] = static_cast<uint32_t>(std::max(int64_t(0), n));
    DISABLE_WARNING_POP
  }
  return {v1, v2};
}

template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
inline void run_scipy_corr(std::string_view method, std::vector<N>& v1, std::vector<N>& v2,
                           std::atomic<double>& rho, std::atomic<double>& pv,
                           std::mutex& data_mutex, std::condition_variable& input_data_cv,
                           std::condition_variable& output_data_cv,
                           std::atomic<bool>& input_data_ready,
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

    if (v1.empty()) {      // EOQ signal
      assert(v2.empty());  // NOLINT
      stdin_stream.pipe().close();
      py.wait();
      return;
    }

    sbuff = fmt::format(FMT_COMPILE("{}\t{}\n"), fmt::join(v1, ","), fmt::join(v2, ","));
    stdin_stream.write(sbuff.data(), static_cast<std::streamsize>(sbuff.size()));
    stdin_stream.flush();

    std::getline(stdout_stream, sbuff);
    const auto sep_idx = sbuff.find('\t');
    assert(sep_idx < sbuff.size());  // NOLINT
    rho = utils::parse_numeric_or_throw<double>(std::string_view(sbuff.data(), sep_idx));
    pv = utils::parse_numeric_or_throw<double>(
        std::string_view(sbuff.data() + sep_idx + 1, sbuff.size() - sep_idx));
    output_data_ready = true;
    output_data_cv.notify_one();
  }
}

template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
[[nodiscard]] inline std::pair<double, double> corr_scipy(const std::vector<N>& v1,
                                                          const std::vector<N>& v2,
                                                          std::string_view method) {
  boost::process::ipstream stdout_stream;
  boost::process::opstream stdin_stream;
  boost::process::spawn(
      boost::process::search_path("python3").string(), "-c", std::string{python_cmds.at(method)},
      boost::process::std_in<stdin_stream, boost::process::std_out> stdout_stream);

  auto buff = fmt::format(FMT_COMPILE("{}\t{}\n"), fmt::join(v1, ","), fmt::join(v2, ","));
  stdin_stream.write(buff.data(), static_cast<std::streamsize>(buff.size()));
  stdin_stream.close();

  std::getline(stdout_stream, buff);
  const auto sep_idx = buff.find('\t');
  assert(sep_idx < buff.size());  // NOLINT
  const auto rho = utils::parse_numeric_or_throw<double>(std::string_view(buff.data(), sep_idx));
  const auto pv = utils::parse_numeric_or_throw<double>(
      std::string_view(buff.data() + sep_idx + 1, buff.size() - sep_idx));

  return {rho, pv};
}
}  // namespace modle::test::correlation
