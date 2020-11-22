#pragma once

#include <algorithm>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <cstdint>
#include <execution>
#include <numeric>
#include <type_traits>
#include <vector>

#include "absl/strings/str_join.h"
#include "fmt/printf.h"
#include "range/v3/numeric.hpp"
#include "range/v3/view.hpp"

namespace modle::correlation::utils {
// Return the indices corresponding to the sorted vector
template <typename Rng>
std::vector<std::size_t> sort_range_by_idx(Rng range) {
  static_assert(ranges::random_access_range<Rng>, "r should be a random access range.");
  std::vector<std::size_t> vi(range.size());
  ranges::iota(vi, 0);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
  std::sort(vi.begin(), vi.end(),
            [&](std::size_t i1, std::size_t i2) { return range[i1] < range[i2]; });
#pragma GCC diagnostic pop
  return vi;
}

// This returns a vector of integers corresponding to the rank of the values of
// vector v
template <typename Rng>
std::vector<double> compute_element_ranks(Rng range) {
  static_assert(ranges::random_access_range<Rng>, "r should be a random access range.");
  const auto vi = sort_range_by_idx(range);
  std::vector<double> vr(vi.size());

  vr[vi[0]] = 0;
  for (auto i = 1U; i < vr.size(); ++i) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
    auto& prev = range[vi[i - 1]];
    auto& current = range[vi[i]];
    if (prev != current) {
      vr[vi[i]] = i;
    } else {
      const auto min = i - 1;
      auto max = i;
      for (auto j = max + 1; j < vr.size() && current == range[vi[j]]; ++j) {
        max = j;
      }
      for (auto j = min; j <= max; ++j) {
        vr[vi[j]] = (max + min) / 2.0;
      }
#pragma GCC diagnostic pop
      i = max;
    }
  }
  return vr;
}
}  // namespace modle::correlation::utils