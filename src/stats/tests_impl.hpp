// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>                                // for min, max
#include <boost/math/distributions/binomial.hpp>    // for cdf, pdf, binomial_distribution
#include <boost/math/distributions/complement.hpp>  // for complement
#include <cassert>                                  // for assert
#include <cmath>                                    // for ceil, floor

#include "modle/common/random.hpp"                      // for binomial_distribution
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...

namespace modle::stats {
// https://github.com/scipy/scipy/blob/47bb6febaa10658c72962b9615d5d5aa2513fa3a/scipy/stats/_binomtest.py#L204
template <binomial_test_alternative alternative, class FP, class I, class>
FP binomial_test(const I k, const I n, const FP p) {
  assert(k >= 0);
  assert(n > 0);
  assert(n >= k);
  assert(p >= FP(0) && p <= FP(1));

  using binom_distr = boost::math::binomial_distribution<FP>;
  auto distr = binom_distr(FP(n), p);

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SHADOW
  auto pmf = [](auto&& distr, const auto k) { return boost::math::pdf(distr, k); };
  auto cdf = [](auto&& distr, const auto k) { return boost::math::cdf(distr, k); };
  auto sf = [](auto&& distr, const auto k) {
    return boost::math::cdf(boost::math::complement(distr, k));
  };
  DISABLE_WARNING_POP

  using ALTERNATIVE = binomial_test_alternative;
  if constexpr (alternative == ALTERNATIVE::LESS) {
    return std::min(FP(1), cdf(distr, k));
  }

  if constexpr (alternative == ALTERNATIVE::GREATER) {
    return std::min(FP(1), sf(distr, k - 1));
  }

  assert(alternative == ALTERNATIVE::TWO_SIDED);
  if (FP(k) == p * FP(n)) {
    return FP(1);
  }

  // https://github.com/scipy/scipy/blob/47bb6febaa10658c72962b9615d5d5aa2513fa3a/scipy/stats/_binomtest.py#L337
  auto binary_search_over_binom_pmf = [&pmf, &distr](auto d, auto x0, auto x1, auto cfx) {
    assert(cfx == 1 || cfx == -1);
    auto fx = [&](auto x) { return cfx * pmf(distr, x); };
    while (x0 < x1) {
      const auto x = x0 + ((x1 - x0) / 2);
      const auto y = fx(x);
      if (y < d) {
        x0 = x + 1;
      } else if (y > d) {
        x1 = x - 1;
      } else {
        return x;
      }
    }

    if (fx(x0) <= d) {
      return x0;
    }
    return x0 - 1;
  };

  const auto d = pmf(distr, k);
  const auto rerr = FP(1) + FP(1e-7);
  const auto m = p * FP(n);

  // https://github.com/scipy/scipy/blob/47bb6febaa10658c72962b9615d5d5aa2513fa3a/scipy/stats/_binomtest.py#L310
  if (FP(k) < m) {
    const auto x0 = static_cast<I>(std::ceil(m));  // lower bound
    const auto x1 = n;                             // upper bound
    const auto ix = binary_search_over_binom_pmf(-d * rerr, x0, x1, -1);
    const auto y = n - ix + I(d * rerr == pmf(distr, ix));
    return std::min(FP(1), cdf(distr, k) + sf(distr, n - y));
  }

  // https://github.com/scipy/scipy/blob/47bb6febaa10658c72962b9615d5d5aa2513fa3a/scipy/stats/_binomtest.py#L320
  const I x0 = 0;                                 // lower bound
  const auto x1 = static_cast<I>(std::floor(m));  // upper bound
  const auto ix = binary_search_over_binom_pmf(d * rerr, x0, x1, 1);
  const auto y = ix + 1;
  return std::min(FP(1), cdf(distr, y - 1) + sf(distr, k - 1));
}

}  // namespace modle::stats

// IWYU pragma: private, include "modle/tests.hpp"
