// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
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

#include "modle/common/common.hpp"  // for u32, i32, i64
#include "modle/common/numeric_utils.hpp"
#include "modle/common/random.hpp"

namespace modle::test::stats {
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
     "method=sub(\"weighted_\", \"\", \"" METHOD "\")\n"           \
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
    absl::flat_hash_set<u32> s;
    while (s.size() < buff.size()) {
      s.insert(dist(rand_eng));
    }
    std::copy(s.begin(), s.end(), buff.begin());
  }
}

template <class N = u32, class = std::enable_if_t<std::is_arithmetic_v<N>>>
[[nodiscard]] inline std::vector<N> generate_random_vect(random::PRNG_t& rand_eng, usize size,
                                                         N min, N max,
                                                         bool allow_duplicates = true) {
  std::vector<N> v(size);
  generate_random_vect(rand_eng, v, min, max, allow_duplicates);
  return v;
}

inline std::pair<std::vector<u32>, std::vector<u32>> generate_correlated_vects(
    random::PRNG_t& rand_eng, u32 size) {
  random::uniform_int_distribution<i32> dist(static_cast<i32>(size) / -50,
                                             static_cast<i32>(size / 50));
  std::vector<u32> v1(size);
  std::vector<u32> v2(size);
  std::iota(v1.begin(), v1.end(), 0);
  std::iota(v2.begin(), v2.end(), 0);
  for (usize i = 0; i < size; ++i) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_USELESS_CAST
    auto n = static_cast<i64>(v1[i]) + dist(rand_eng);
    v1[i] = static_cast<u32>(std::max(i64(0), n));
    n = static_cast<i64>(v2[i]) + dist(rand_eng);
    v2[i] = static_cast<u32>(std::max(i64(0), n));
    DISABLE_WARNING_POP
  }
  return {v1, v2};
}

template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
inline void run_external_corr(std::string_view method, std::vector<N>& v1, std::vector<N>& v2,
                              std::vector<double>& weights, std::atomic<double>& rho,
                              std::atomic<double>& pv, std::mutex& data_mutex,
                              std::condition_variable& input_data_cv,
                              std::condition_variable& output_data_cv,
                              std::atomic<bool>& input_data_ready,
                              std::atomic<bool>& output_data_ready) {
  boost::process::ipstream stdout_stream;
  boost::process::ipstream stderr_stream;
  boost::process::opstream stdin_stream;

  try {
    auto c = [&]() {
      if (absl::StartsWith(method, "weighted_"sv)) {
        return boost::process::child(
            boost::process::search_path("Rscript").string(), "--quiet", "-e",
            std::string{external_cmds.at(method)},
            boost::process::std_in<stdin_stream, boost::process::std_out> stdout_stream,
            boost::process::std_err > stderr_stream);
      }

      return boost::process::child(
          boost::process::search_path("python3").string(), "-c",
          std::string{external_cmds.at(method)},
          boost::process::std_in<stdin_stream, boost::process::std_out> stdout_stream,
          boost::process::std_err > stderr_stream);
    }();
    assert(c.running());

    std::string sbuff;
    while (true) {
      {
        std::unique_lock<std::mutex> l(data_mutex);
        input_data_cv.wait(l, [&]() { return input_data_ready.load(); });
        input_data_ready = false;
      }

      if (v1.empty()) {  // EOQ signal
        assert(v2.empty());
        stdin_stream.pipe().close();
        c.wait();
        return;
      }

      sbuff = [&]() {
        if (absl::StartsWith(method, "weighted_")) {
          return fmt::format(FMT_COMPILE("{}\t{}\t{}\n"), fmt::join(v1, ","), fmt::join(v2, ","),
                             fmt::join(weights, ","));
        }
        return fmt::format(FMT_COMPILE("{}\t{}\n"), fmt::join(v1, ","), fmt::join(v2, ","));
      }();
      stdin_stream.write(sbuff.data(), static_cast<std::streamsize>(sbuff.size()));
      stdin_stream.flush();

      std::getline(stdout_stream, sbuff);
      const auto sep_idx = sbuff.find('\t');
      assert(sep_idx < sbuff.size());
      rho = utils::parse_numeric_or_throw<double>(std::string_view(sbuff.data(), sep_idx));
      pv = utils::parse_numeric_or_throw<double>(
          std::string_view(sbuff.data() + sep_idx + 1, sbuff.size() - sep_idx));
      output_data_ready = true;
      output_data_cv.notify_one();
    }
  } catch (const std::exception& e) {
    std::string buff1;
    std::string buff2;
    while (std::getline(stderr_stream, buff1)) {
      buff2.append(buff1);
      buff2.append("\n");
    }

    if (!buff2.empty()) {
      throw std::runtime_error(fmt::format(FMT_STRING("{}:\n{}"), e.what(), buff2));
    }
    throw;
  }
}

template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
[[nodiscard]] inline std::pair<double, double> external_corr(const std::vector<N>& v1,
                                                             const std::vector<N>& v2,
                                                             const std::vector<double>& w,
                                                             std::string_view method) {
  boost::process::ipstream stdout_stream;
  boost::process::opstream stdin_stream;
  if (absl::StartsWith(method, "weighted_"sv)) {
    boost::process::spawn(
        boost::process::search_path("Rscript").string(), "--quiet", "-e",
        std::string{external_cmds.at(method)},
        boost::process::std_in<stdin_stream, boost::process::std_out> stdout_stream);
  } else {
    boost::process::spawn(
        boost::process::search_path("python3").string(), "-c",
        std::string{external_cmds.at(method)},
        boost::process::std_in<stdin_stream, boost::process::std_out> stdout_stream);
  }

  std::string buff;
  if (absl::StartsWith(method, "weighted_"sv)) {
    buff = fmt::format(FMT_COMPILE("{}\t{}\t{}\n"), fmt::join(v1, ","), fmt::join(v2, ","),
                       fmt::join(w, ","));
  } else {
    buff = fmt::format(FMT_COMPILE("{}\t{}\n"), fmt::join(v1, ","), fmt::join(v2, ","));
  }
  stdin_stream.write(buff.data(), static_cast<std::streamsize>(buff.size()));
  stdin_stream.close();

  std::getline(stdout_stream, buff);
  const auto sep_idx = buff.find('\t');
  assert(sep_idx < buff.size());
  const auto rho = utils::parse_numeric_or_throw<double>(std::string_view(buff.data(), sep_idx));
  const auto pv = utils::parse_numeric_or_throw<double>(
      std::string_view(buff.data() + sep_idx + 1, buff.size() - sep_idx));

  return {rho, pv};
}
}  // namespace modle::test::stats
