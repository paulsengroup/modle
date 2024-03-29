// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "modle/bigwig.hpp"

#include <absl/types/span.h>  // for Span, MakeSpan
#include <fmt/format.h>       // for FMT_STRING, format

#include <cassert>      // for assert
#include <limits>       // for numeric_limits
#include <stdexcept>    // for runtime_error
#include <string>       // for string
#include <string_view>  // for string_view
#include <type_traits>  // for is_arithmetic, is_signed
#include <vector>       // for vector

#include "libbigwig/bigWig.h"                           // for bwAddIntervalSpanSteps
#include "modle/common/common.hpp"                      // for u32, u64, i32
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...

namespace modle::io::bigwig {

constexpr Reader::operator bool() const noexcept { return this->_fp != nullptr; }

constexpr const std::filesystem::path& Reader::path() const noexcept { return this->_fname; }
constexpr auto Reader::chromosomes() const noexcept -> const Chromosomes& { return this->_chroms; }

template <Reader::StatsType stat>
inline double Reader::stats(const std::string& chrom, bp_t start, bp_t end) {
  static_assert(stat != bwStatsType::doesNotExist);
  assert(!!*this);
  this->validate_query(chrom, start, end);

  const std::unique_ptr<double, decltype(&std::free)> stats{
      bwStatsFromFull(this->_fp, chrom.c_str(), utils::conditional_static_cast<u32>(start),
                      utils::conditional_static_cast<u32>(end), 1, stat),
      &std::free};
  if (!stats) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return *stats;
}

template <Reader::StatsType stat>
inline void Reader::stats(const std::string& chrom, bp_t start, bp_t end, bp_t window_size,
                          std::vector<double>& buff) {
  static_assert(stat != bwStatsType::doesNotExist);
  assert(window_size != 0);
  assert(!!*this);
  this->validate_query(chrom, start, end);

  const auto num_bins =
      utils::conditional_static_cast<u32>((end - start + window_size - 1) / window_size);

  buff.clear();
  const std::unique_ptr<double, decltype(&std::free)> stats{
      bwStatsFromFull(this->_fp, chrom.c_str(), utils::conditional_static_cast<u32>(start),
                      utils::conditional_static_cast<u32>(end), num_bins, stat),
      &std::free};
  if (!stats) {
    return;
  }

  buff.resize(num_bins);
  std::copy_n(stats.get(), buff.size(), buff.begin());
}

template <Reader::StatsType stat>
inline std::vector<double> Reader::stats(const std::string& chrom, bp_t start, bp_t end,
                                         bp_t window_size) {
  std::vector<double> buff{};
  this->stats<stat>(chrom, start, end, window_size, buff);
  return buff;
}

constexpr Writer::operator bool() const noexcept { return this->_fp != nullptr; }

template <class Chromosomes>
inline void Writer::write_chromosomes(Chromosomes& chroms) {
  const auto num_chroms = utils::conditional_static_cast<usize>(chroms.size());

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SHADOW
  auto chrom_names = [&]() {
    std::vector<char*> chrom_names(num_chroms);
    std::transform(chroms.begin(), chroms.end(), chrom_names.begin(),
                   [](auto& chrom) { return chrom.first.data(); });
    return chrom_names;
  }();

  auto chrom_sizes = [&]() {
    std::vector<u32> chrom_sizes(num_chroms);

    using chrom_size_t = decltype(chroms.begin()->second);
    std::transform(chroms.begin(), chroms.end(), chrom_sizes.begin(), [](auto& chrom) {
      if constexpr (const auto max_val = (std::numeric_limits<u32>::max)();
                    (std::numeric_limits<chrom_size_t>::max)() > max_val) {
        if (chrom.second > max_val) {
          throw std::runtime_error(fmt::format(
              FMT_STRING(
                  "We currently don't support writing chromosomes longer than ~4.29 Gbp (2^32 "
                  "bp) to bigWig files. \"{}\" has a length of {:.4f} Gbp"),
              chrom.first, static_cast<double>(chrom.second) / 1e9));
        }
      }
      if constexpr (std::is_signed<chrom_size_t>()) {
        if (chrom.second < 0) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("\"{}\" appears to have a negative length of {} bp. This is not allowed"),
              chrom.first, chrom.second));
        }
      }
      return static_cast<u32>(chrom.second);
    });

    return chrom_sizes;
  }();
  DISABLE_WARNING_POP

  this->write_chromosomes(chrom_names.data(), chrom_sizes.data(), num_chroms);
}

template <class N, class>
inline void Writer::write_range(std::string_view chrom_name, const absl::Span<N> values, u64 span,
                                u64 step, u64 offset) {
  assert(this->_initialized);
  assert(this->_fp);
  std::vector<float> fvalues;
  auto fvalues_span = [&]() {
    if constexpr (std::is_same_v<N, float>) {
      return values;
    } else {
      fvalues.resize(values.size());
      std::transform(values.begin(), values.end(), fvalues.begin(),
                     [](const auto n) { return static_cast<float>(n); });

      return absl::MakeSpan(fvalues);
    }
  }();

  auto chrom_name_tmp = std::string{chrom_name};
  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  if (bwAddIntervalSpanSteps(this->_fp, chrom_name_tmp.data(), static_cast<u32>(offset),
                             static_cast<u32>(span), static_cast<u32>(step), fvalues_span.data(),
                             static_cast<u32>(fvalues_span.size()))) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to write data for chrom \"{}\""), chrom_name));
  }
}

constexpr const std::filesystem::path& Writer::path() const noexcept { return this->_fname; }

}  // namespace modle::io::bigwig
