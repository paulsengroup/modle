// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>                           // for transform
#include <boost/math/constants/constants.hpp>  // for pi
#include <cassert>                             // for assert
#include <cmath>                               // for exp, sqrt
#include <vector>                              // for vector

#include "modle/common/common.hpp"  // for isize, usize

namespace modle::stats {

// https://patrickfuller.github.io/gaussian-blur-image-processing-for-scientists-and-engineers-part-4/
template <class FP, class>
void compute_gauss_kernel(std::vector<FP>& buff, const FP sigma, const FP cutoff) {
  const auto size = static_cast<usize>(1 + 2 * std::sqrt(-2 * sigma * sigma * std::log(cutoff)));
  compute_gauss_kernel(size + (size % 2 == 0), buff, sigma);
}

template <class FP, class>
std::vector<FP> compute_gauss_kernel(const FP sigma, const FP cutoff) {
  std::vector<FP> v;
  compute_gauss_kernel(v, sigma, cutoff);
  return v;
}

template <class FP, class>
std::vector<FP> compute_gauss_kernel(const usize size, const FP sigma) {
  std::vector<FP> v(size * size);
  compute_gauss_kernel(size, v, sigma);
  return v;
}

template <class FP, class>
void compute_gauss_kernel(const usize size, std::vector<FP>& buff, const FP sigma) {
  assert(size % 2 != 0);
  buff.resize(size * size);
  const auto q = FP(2) * sigma * sigma;
  const auto denom = q * boost::math::constants::pi<FP>();

  const auto i1 = static_cast<isize>(size / 2);
  const auto i0 = i1 * isize(-1);

  FP sum = 0;
  for (auto i = i0; i <= i1; ++i) {
    for (auto j = i0; j <= i1; ++j) {
      const auto p = std::sqrt(static_cast<FP>(i * i + j * j));
      const auto ii = static_cast<usize>(i + i1);
      const auto jj = static_cast<usize>(j + i1);
      const auto idx = (ii * size) + jj;
      buff[idx] = (std::exp(-(p * p) / q)) / denom;
      sum += buff[idx];
    }
  }

  std::transform(buff.begin(), buff.end(), buff.begin(), [&](const FP n) { return n / sum; });
}

}  // namespace modle::stats
