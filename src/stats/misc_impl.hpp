// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <cassert>
#include <cmath>
#include <vector>

#include "modle/common/common.hpp"

namespace modle::stats {

// https://github.com/scipy/scipy/blob/v1.9.3/scipy/ndimage/_filters.py#L256-L260
template <class FP, class>
void compute_gauss_kernel2d(std::vector<FP>& buff, const FP sigma, const FP truncate) {
  assert(truncate > 0.0);
  assert(sigma >= 0.0);
  const auto radius = static_cast<usize>(truncate * sigma + FP(0.5));
  compute_gauss_kernel2d(radius, buff, sigma);
}

template <class FP, class>
std::vector<FP> compute_gauss_kernel2d(const FP sigma, const FP truncate) {
  assert(truncate > 0.0);
  std::vector<FP> v;
  compute_gauss_kernel2d(v, sigma, truncate);
  return v;
}

template <class FP, class>
std::vector<FP> compute_gauss_kernel2d(const usize radius, const FP sigma) {
  const auto size = 2 * radius + 1;
  std::vector<FP> v(size * size);
  compute_gauss_kernel2d(radius, v, sigma);
  return v;
}

namespace internal {

template <class FP>
[[nodiscard]] inline std::vector<FP> compute_gauss_kernel1d(const usize size, const FP sigma) {
  // Inspired by https://codereview.stackexchange.com/a/169675
  std::vector<FP> buff(size);

  const auto spread = FP(1) / (2 * sigma * sigma);
  const auto mid = FP(size - 1) / FP(2);

  for (usize i = 0; i < size; ++i) {
    const auto n = FP(i) - mid;
    buff[i] = std::exp(-n * n * spread);
  }

  return buff;
}
}  // namespace internal

template <class FP, class>
void compute_gauss_kernel2d(const usize radius, std::vector<FP>& buff, const FP sigma) {
  // Inspired by https://codereview.stackexchange.com/a/169675
  const auto size = 2 * radius + 1;
  const auto gauss_kernel1d = internal::compute_gauss_kernel1d<FP>(size, sigma);

  buff.resize(size * size);

  FP sum = 0;
  for (usize i = 0; i < size; ++i) {
    for (usize j = 0; j < size; ++j) {
      const auto k = (i * size) + j;
      buff[k] = gauss_kernel1d[i] * gauss_kernel1d[j];
      sum += buff[k];
    }
  }

  std::transform(buff.begin(), buff.end(), buff.begin(), [&](const FP n) { return n / sum; });
}

template <class InputIt1, class InputIt2, class N, class>
constexpr N cross_correlation(InputIt1 kernel_first, InputIt1 kernel_last, InputIt2 buff_first) {
  N tot = 0;
  for (; kernel_first != kernel_last; ++kernel_first, ++buff_first) {
    tot += utils::conditional_static_cast<N>(*kernel_first) *
           utils::conditional_static_cast<N>(*buff_first);
  }

  return tot;
}

template <class Rng1, class Rng2, class N, class>
constexpr N cross_correlation(const Rng1& kernel, const Rng2& buff) {
  assert(kernel.size() == buff.size());
  return cross_correlation(kernel.begin(), kernel.end(), buff.begin());
}

}  // namespace modle::stats
