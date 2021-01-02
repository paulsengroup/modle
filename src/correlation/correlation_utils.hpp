#pragma once

#include <absl/strings/str_join.h>
#include <absl/types/span.h>

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <type_traits>
#include <vector>

#include "modle/suppress_compiler_warnings.hpp"

namespace modle::correlation::utils {

// Return the indices corresponding to the sorted vector
template <typename N>  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
inline std::vector<std::size_t> sort_range_by_idx(absl::Span<const N> v) {
  static_assert(std::is_arithmetic<N>::value,
                "v should be convertible to a span of numeric types.");
  std::vector<std::size_t> vi(v.size());
  std::iota(vi.begin(), vi.end(), 0);

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_CONVERSION
  std::sort(vi.begin(), vi.end(), [&](std::size_t i1, std::size_t i2) { return v[i1] < v[i2]; });
  DISABLE_WARNING_POP
  return vi;
}

template <typename N>  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
inline std::vector<std::size_t> sort_range_by_idx(const std::vector<N>& v) {
  return sort_range_by_idx(absl::MakeConstSpan(v));
}

// This returns a vector of integers corresponding to the rank of the values of
// vector v
template <typename N>
inline std::vector<double> compute_element_ranks(absl::Span<const N> v) {
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
inline std::vector<double> compute_element_ranks(const std::vector<N>& v) {
  return compute_element_ranks(absl::MakeConstSpan(v));
}

}  // namespace modle::correlation::utils