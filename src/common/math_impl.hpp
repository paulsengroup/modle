// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <boost/math/distributions/binomial.hpp>  // for binomial
#include <cmath>                                  // for sqrt
#include <cstddef>                                // for size_t, ptrdiff_t
#include <iterator>                               // for distance
#include <numeric>                                // for accumulate
#include <type_traits>
#include <vector>  // for vector

namespace modle::math {
template <class T, class>
constexpr T abs_diff(const T a, const T b) noexcept {
  return a > b ? a - b : b - a;
}

template <class InputIt, class UnaryOperation, class>
double mean(InputIt begin, InputIt end, UnaryOperation op) {
  if (begin == end) {
    return 0.0;
  }
  return std::accumulate(
             begin, end, 0.0,
             [&](auto accumulator, const auto& val) { return accumulator + double(op(val)); }) /
         static_cast<double>(std::distance(begin, end));
}

template <class InputIt, class OutputIt, class I, class UnaryOperation, class>
size_t moving_average(InputIt input_begin, InputIt input_end, OutputIt output_begin,
                      const I window_size, UnaryOperation op) {
  using OutputItValueT = typename OutputIt::value_type;
  if (static_cast<std::ptrdiff_t>(window_size) >= std::distance(input_begin, input_end)) {
    *output_begin = OutputItValueT(math::mean(input_begin, input_end, op));
    return 1;
  }

  for (auto it1 = input_begin, it2 = input_begin + static_cast<std::ptrdiff_t>(window_size);
       it2 != input_end; ++it1, ++it2, ++output_begin) {
    *output_begin = OutputItValueT(mean(it1, it2, op));
  }
  return size_t(std::distance(input_begin, input_end)) - size_t(window_size);
}

template <class InputIt, class FP, class UnaryOperation, class>
double sum_of_squared_deviations(InputIt input_begin, InputIt input_end, const FP sample_mean,
                                 UnaryOperation op) {
  if (input_begin == input_end) {
    return 0.0;
  }
  return std::accumulate(input_begin, input_end, 0.0,
                         [&, sm = double(sample_mean)](const auto accumulator, const auto& val) {
                           const auto n = double(op(val)) - sm;
                           return accumulator + (n * n);
                         });
}

template <class InputIt, class UnaryOperation, class>
double sum_of_squared_deviations(InputIt input_begin, InputIt input_end, UnaryOperation op) {
  return math::sum_of_squared_deviations(input_begin, input_end,
                                         math::mean(input_begin, input_end, op), op);
}

template <class InputIt, class FP, class UnaryOperation, class>
double variance(InputIt input_begin, InputIt input_end, const FP sample_mean, UnaryOperation op) {
  if (input_begin == input_end) {
    return 0.0;
  }
  return math::sum_of_squared_deviations(input_begin, input_end, sample_mean, op) /
         static_cast<double>(std::distance(input_begin, input_end));
}

template <class InputIt, class UnaryOperation, class>
double variance(InputIt input_begin, InputIt input_end, UnaryOperation op) {
  return math::variance(input_begin, input_end, math::mean(input_begin, input_end, op), op);
}

template <class InputIt, class FP, class UnaryOperation, class>
double standard_dev(InputIt input_begin, InputIt input_end, const FP sample_mean,
                    UnaryOperation op) {
  if (input_begin == input_end) {
    return 0.0;
  }
  return std::sqrt(math::variance(input_begin, input_end, sample_mean, op));
}

template <class InputIt, class UnaryOperation, class>
double standard_dev(InputIt input_begin, InputIt input_end, UnaryOperation op) {
  return math::standard_dev(input_begin, input_end, math::mean(input_begin, input_end, op), op);
}

// https://github.com/scipy/scipy/blob/47bb6febaa10658c72962b9615d5d5aa2513fa3a/scipy/stats/_binomtest.py#L204
template <binomial_test_alternative alternative, class FP, class I, class>
FP binomial_test(const I k, const I n, const FP p) {
  assert(k >= 0);                    // NOLINT
  assert(n > 0);                     // NOLINT
  assert(n >= k);                    // NOLINT
  assert(p >= FP(0) && p <= FP(1));  // NOLINT

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

  assert(alternative == ALTERNATIVE::TWO_SIDED);  // NOLINT
  if (FP(k) == p * FP(n)) {
    return FP(1);
  }

  // https://github.com/scipy/scipy/blob/47bb6febaa10658c72962b9615d5d5aa2513fa3a/scipy/stats/_binomtest.py#L337
  auto binary_search_over_binom_pmf = [&pmf, &distr](auto d, auto x0, auto x1, auto cfx) {
    assert(cfx == FP(1) || cfx == FP(-1));  // NOLINT
    auto fx = [&](auto x) { return cfx * pmf(distr, x); };
    while (x0 < x1) {
      const auto x = x0 + std::trunc((x1 - x0) / 2);
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

  auto d = pmf(distr, k);
  const auto rerr = FP(1) + FP(1e-7);
  // https://github.com/scipy/scipy/blob/47bb6febaa10658c72962b9615d5d5aa2513fa3a/scipy/stats/_binomtest.py#L310
  if (FP(k) < p * FP(n)) {
    d *= -rerr;                      // lower bound
    auto x0 = std::ceil(p * FP(n));  // upper bound
    auto x1 = FP(n);
    const auto ix = binary_search_over_binom_pmf(d, x0, x1, FP(-1));
    const auto y = FP(n) - ix + FP(d * rerr == pmf(distr, ix));
    return std::min(FP(1), cdf(distr, k) + sf(distr, FP(n) - y));
  }

  // https://github.com/scipy/scipy/blob/47bb6febaa10658c72962b9615d5d5aa2513fa3a/scipy/stats/_binomtest.py#L320
  d *= rerr;        // lower bound
  auto x0 = FP(0);  // upper bound
  auto x1 = std::floor(p * FP(n));
  const auto ix = binary_search_over_binom_pmf(d, x0, x1, FP(1));
  const auto y = ix + 1;
  return std::min(FP(1), cdf(distr, y - 1) + sf(distr, k - 1));
}
}  // namespace modle::math

// IWYU pragma: private, include "modle/math.hpp"
