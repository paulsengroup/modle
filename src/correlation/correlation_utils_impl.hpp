// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private

#include <absl/types/span.h>  // for MakeConstSpan, Span

#include <algorithm>    // for sort
#include <numeric>      // for iota
#include <type_traits>  // for is_arithmetic
#include <vector>       // for vector

#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNING_PUSH

namespace modle::correlation::utils {

// Return the indices corresponding to the sorted vector
template <typename N>
std::vector<usize> sort_range_by_idx(absl::Span<const N> v) {
  static_assert(std::is_arithmetic<N>::value,
                "v should be convertible to a span of numeric types.");
  std::vector<usize> vi(v.size());
  std::iota(vi.begin(), vi.end(), 0);

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_CONVERSION
  std::sort(vi.begin(), vi.end(), [&](usize i1, usize i2) { return v[i1] < v[i2]; });
  DISABLE_WARNING_POP
  return vi;
}

template <typename N>
std::vector<usize> sort_range_by_idx(const std::vector<N>& v) {
  return sort_range_by_idx(absl::MakeConstSpan(v));
}

// This returns a vector of integers corresponding to the rank of the values of
// vector v
template <typename N>
std::vector<double> compute_element_ranks(absl::Span<const N> v) {
  static_assert(std::is_arithmetic<N>::value,
                "v should be convertible to a span of numeric types.");
  const auto vi = sort_range_by_idx(v);
  std::vector<double> vr(vi.size());

  vr[vi[0]] = 0;
  for (auto i = 1U; i < vr.size(); ++i) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_SIGN_CONVERSION
    auto& prev = v[vi[i - 1]];
    auto& current = v[vi[i]];
    if (prev != current) {
      vr[vi[i]] = i;
    } else {
      const auto min = i - 1;
      auto max = i;
      for (auto j = max + 1; j < vr.size() && current == v[vi[j]]; ++j) {
        max = j;
      }
      for (auto j = min; j <= max; ++j) {
        // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
        vr[vi[j]] = (max + min) / 2.0;
      }
      DISABLE_WARNING_POP
      i = max;
    }
  }
  return vr;
}

template <typename N>
std::vector<double> compute_element_ranks(const std::vector<N>& v) {
  return compute_element_ranks(absl::MakeConstSpan(v));
}

}  // namespace modle::correlation::utils
