// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "modle/bigwig.hpp"

#include <absl/container/flat_hash_map.h>  // for flat_hash_map
#include <absl/strings/match.h>            // for StartsWith
#include <fmt/format.h>                    // for FMT_STRING, format
#include <fmt/ostream.h>

#include <iterator>     // for pair
#include <limits>       // for numeric_limits
#include <stdexcept>    // for runtime_error
#include <string>       // for string
#include <string_view>  // for string_view
#include <type_traits>  // for is_arithmetic, is_signed
#include <vector>       // for vector

#include "libBigWig/bigWig.h"       // for bwAddIntervalSpanSteps
#include "modle/common/common.hpp"  // for u32, u64, i32

namespace modle::io::bigwig {

template <class Chromosomes>
void Writer::write_chromosomes(const Chromosomes& chroms) {
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  const auto num_chroms = static_cast<usize>(chroms.size());
  DISABLE_WARNING_POP

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SHADOW
  const auto chrom_names = [&]() {
    std::vector<const char*> chrom_names(num_chroms);
    std::transform(chroms.begin(), chroms.end(), chrom_names.begin(),
                   [](const auto& chrom) { return chrom.first.c_str(); });
    return chrom_names;
  }();

  const auto chrom_sizes = [&]() {
    std::vector<u32> chrom_sizes(num_chroms);

    using chrom_size_t = decltype(chroms.begin()->second);
    std::transform(chroms.begin(), chroms.end(), chrom_sizes.begin(), [](const auto& chrom) {
      if constexpr (const auto max_val = (std::numeric_limits<u32>::max)();
                    (std::numeric_limits<chrom_size_t>::max)() > max_val) {
        if (chrom.second > max_val) {
          throw std::runtime_error(fmt::format(
              FMT_STRING(
                  "We currently don't support writing chromosomes longer than ~4.29 Gbp (2^32 "
                  "bp) to bigWig files. '{}' has a length of {:.4f} Gbp"),
              chrom.first, static_cast<double>(chrom.second) / 1e9));  // NOLINT
        }
      }
      if constexpr (std::is_signed<chrom_size_t>()) {
        if (chrom.second < 0) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("'{}' appears to have a negative length of {} bp. This is not allowed"),
              chrom.first, chrom.second));
        }
      }
      return static_cast<u32>(chrom.second);
    });

    return chrom_sizes;
  }();
  DISABLE_WARNING_POP

  this->write_chromosomes(chrom_names.data(), chrom_sizes.data(), chroms.size());
}

template <class Str>
void Writer::write_chromosomes(const absl::Span<Str> chrom_names,
                               const absl::Span<const u32> chrom_sizes) {
  assert(chrom_names.size() == chrom_sizes.size());  // NOLINT
  if constexpr (std::is_same_v<std::remove_cv_t<Str>, char**>) {
    write_chromosomes(chrom_names.data(), chrom_sizes.data(), chrom_names.size());
  }

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SHADOW
  const auto chrom_names_c = [&chrom_names]() {
    std::vector<const char*> chrom_names_c(chrom_names.size());
    std::transform(chrom_names.begin(), chrom_names.end(), chrom_names_c.begin(),
                   [](const auto& name) { return name.c_str(); });
    return chrom_names;
  }();
  DISABLE_WARNING_POP

  this->write_chromosomes(chrom_names_c.data(), chrom_sizes.data(), chrom_names_c.size());
}

template <class N, class>
void Writer::write_range(std::string_view chrom_name, const absl::Span<N> values, u64 span,
                         u64 step) {
  assert(this->_initialized);  // NOLINT
  assert(this->_fp);           // NOLINT
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
  if (bwAddIntervalSpanSteps(this->_fp, chrom_name_tmp.data(), this->_offset,
                             static_cast<u32>(span), static_cast<u32>(step), fvalues_span.data(),
                             static_cast<u32>(fvalues_span.size()))) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to write data for chrom \"{}\""), chrom_name));
  }
}

}  // namespace modle::io::bigwig
