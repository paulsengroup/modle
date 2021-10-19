#include <algorithm>                              // for min
#include <boost/math/distributions/binomial.hpp>  // for binomial_distribution
#include <cassert>

namespace modle::stats {
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
}  // namespace modle::stats

// IWYU pragma: private, include "modle/tests.hpp"
