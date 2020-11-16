#pragma once
#include <algorithm>
#include <atomic>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <execution>
#include <numeric>
#include <vector>

#include "./correlation_utils.hpp"
#include "modle/correlation.hpp"

namespace modle::correlation {
template <typename N>
CorrelationTest<N>::CorrelationTest(const std::vector<N> &v1, const std::vector<N> &v2)
    : _v1(v1), _v2(v2) {
  assert(v1.size() == v2.size());
}

template <typename N>
double CorrelationTest<N>::compute_spearman_significance(double rho, std::size_t n) {
  assert(n > 2);
  const auto dof = static_cast<const double>(n - 2);
  const double tscore = rho * std::sqrt(dof / ((1.0 + rho) * (1.0 - rho)));
  boost::math::students_t_distribution<double> dist(dof);
  // See line 4229 at
  // https://github.com/scipy/scipy/blob/6703631bcd15750e86f4098b0421efabcac0f7c2/scipy/stats/stats.py
  // 2 * survival function
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  return 2.0 * boost::math::cdf<double>(boost::math::complement(dist, std::fabs<double>(tscore)));
}

template <typename N>
template <typename Iterator, typename>
std::pair<double, double> CorrelationTest<N>::compute_spearman(const Iterator v1_b,
                                                               const Iterator v1_e,
                                                               const Iterator v2_b) {
  const uint32_t size = std::distance(v1_b, v1_e);
  if (std::accumulate(v1_b, v1_e, 0.0) == 0 || std::accumulate(v2_b, v2_b + size, 0.0) == 0) {
    return {0, 0};
  }
  auto v1r = utils::compute_element_ranks(v1_b, v1_b + size);
  auto v2r = utils::compute_element_ranks(v2_b, v2_b + size);

  auto v1r_avg = std::reduce(v1r.begin(), v1r.end(), 0.0) / size;
  auto v2r_avg = std::reduce(v2r.begin(), v2r.end(), 0.0) / size;
  double n = 0;  // numerator: sum (x - xm) * (y - ym)
  // denominator: sum (x - xm)^2 * sum (y - ym)^2
  double d1 = 0;
  double d2 = 0;

  for (auto i = 0U; i < size; ++i) {
    const auto &r1 = v1r[i];
    const auto &r2 = v2r[i];
    n += (r1 - v1r_avg) * (r2 - v2r_avg);
    d1 += std::pow(r1 - v1r_avg, 2);
    d2 += std::pow(r2 - v2r_avg, 2);
  }
  auto rho = n / std::sqrt(d1 * d2);
  if (std::isnan(rho)) {
    throw std::logic_error("compute_spearman: rho cannot be nan!");
  }

  return {rho, compute_spearman_significance(rho, size)};
};

template <typename N>
double CorrelationTest<N>::compute_kendall_b_significance(uint32_t nc, uint32_t nd, uint32_t size,
                                                          uint32_t tie1, uint32_t tie2) {
  const uint32_t tot = size * (size - 1) / 2;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  if ((tie1 == 0 && tie2 == 0) && (size <= 33 || std::min(nd, tot - nd) <= 1)) {
    // Exact p-value, see p. 68 of Maurice G. Kendall, "Rank Correlation Methods" (4th Edition),
    // Charles Griffin & Co., 1970.
    auto c = std::min(nd, tot - nd);
    assert(2 * c < size * (size - 1));
    if (size == 1 || size == 2 || (2 * c) == tot) {
      return 1.0;
    }
    if (c == 0 || c == 1) {
      // std::tgamma is used to approximate factorial
      // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
      return 2.0 / std::tgamma(size - (c == 1));  // NOLINT(readability-implicit-bool-conversion)
    }
    std::vector<double> v_new(c + 1, 0.0);
    v_new[0] = 1;
    v_new[1] = 1;
    for (auto i = 3U; i <= size; ++i) {
      auto v_old = v_new;
      for (auto j = 1U; j < std::min(i, c + 1); ++j) {
        v_new[j] += v_new[j - 1];
      }
      for (auto j = i; j <= c; ++j) {
        v_new[j] += v_new[j - 1] - v_old[j - i];
      }
    }
    // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
    return 0.2 * std::reduce(v_new.begin(), v_new.end(), 0.0) /
           std::tgamma(size);  // looks like pv is multiplied by 10 for some reason
  }
  const uint64_t x0 = tie1 * (tie1 - 1) * (tie1 - 2);
  const uint64_t y0 = tie2 * (tie2 - 1) * (tie2 - 2);
  const uint64_t x1 = tie1 * (tie1 - 1) * (2 * tie1 + 5);
  const uint64_t y1 = tie2 * (tie2 - 1) * (2 * tie2 + 5);
  static constexpr auto sqrt2 = boost::math::constants::root_two<double>();
  const double var = (size * (size - 1) * (2.0 * size + 5) - x1 - y1) / 18.0 +
                     (2.0 * tie1 * tie2) / (size * (size - 1)) +
                     x0 * y0 / (9.0 * size * (size - 1) * (size - 2));
  return boost::math::erfc<double>(std::abs<double>(static_cast<int64_t>(nc) - nd) /
                                   std::sqrt(var) / sqrt2);
}
template <typename N>
template <typename Iterator, typename>
std::pair<double, double> CorrelationTest<N>::compute_kendall_b_n2(const Iterator v1_b,
                                                                   const Iterator v1_e,
                                                                   const Iterator v2_b) {
  const std::size_t size = std::distance(v1_b, v1_e);
  const uint32_t n_pairs = size * (size - 1) / 2;
  double tau = 0.0;

  std::atomic<uint32_t> tie1 = 0;  // n of ties in v1
  std::atomic<uint32_t> tie2 = 0;  // n of ties in v2
  std::atomic<int64_t> nc = 0;     // n of concordant ranks
  std::atomic<int64_t> nd = 0;     // n of discordant ranks

  std::vector<uint32_t> idx(size);
  std::iota(idx.begin(), idx.end(), 0);

  auto ktb_fx = [&](uint32_t i) {
    uint32_t m1_local = 0;   // n of ties in v1
    uint32_t m2_local = 0;   // n of ties in v2
    int64_t nc_local = 0;    // n of concordant ranks
    int64_t nd_local = 0;    // n of discordant ranks
    auto v1i = *(v1_b + i);  // Equivalent to v1[i]
    auto v2i = *(v2_b + i);  // Equivalent to v2[i]
    for (auto j = i + 1; j < size; ++j) {
      auto v1j = *(v1_b + j);  // Equivalent to v1[j]
      auto v2j = *(v2_b + j);  // Equivalent to v2[j]
      if (v2i > v2j) {
        nc_local += v1i > v1j;
        nd_local += v1i < v1j;
        m1_local += v1i == v1j;
      } else if (v2i < v2j) {
        nd_local += v1i > v1j;
        nc_local += v1i < v1j;
        m1_local += v1i == v1j;
      } else {
        ++m2_local;
        m1_local += v1i == v1j;
      }
    }
    nc.fetch_add(nc_local, std::memory_order_relaxed);
    nd.fetch_add(nd_local, std::memory_order_relaxed);
    tie1.fetch_add(m1_local, std::memory_order_relaxed);
    tie2.fetch_add(m2_local, std::memory_order_relaxed);
  };

  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  if (size < 1000) {  // Don't bother running for_each in parallel for small arrays
    std::for_each(std::execution::seq, idx.begin(), idx.end(), ktb_fx);
  } else {
    std::for_each(std::execution::par, idx.begin(), idx.end(), ktb_fx);
  }

  if (tie1 < n_pairs && tie2 < n_pairs) {
    tau = static_cast<double>(nc - nd) / (std::sqrt(n_pairs - tie1) * std::sqrt(n_pairs - tie2));
  }

  const auto pv = compute_kendall_b_significance(nc, nd, size, tie1, tie2);
  return {tau, pv};
}

template <typename N>
template <typename Iterator, typename>
std::pair<double, double> CorrelationTest<N>::compute_kendall(const Iterator v1_b,
                                                              const Iterator v1_e,
                                                              const Iterator v2_b) {
  return compute_kendall_b_n2(v1_b, v1_e, v2_b);
}

template <typename N>
std::pair<double, double> CorrelationTest<N>::compute_kendall() const {
  return this->compute_kendall(this->_v1.begin(), this->_v1.end(), this->_v2.begin());
}

template <typename N>
std::pair<double, double> CorrelationTest<N>::compute_spearman() const {
  return this->compute_spearman(this->_v1.begin(), this->_v1.end(), this->_v2.begin());
}

template <typename N>
std::vector<std::pair<double, double>> CorrelationTest<N>::compute_spearman(
    uint64_t window_size, uint64_t window_overlap) const {
  std::vector<uint32_t> windows_start((this->_v1.size() - window_size) / window_overlap);
  for (auto i = 0UL; i < windows_start.size(); ++i) {
    windows_start[i] = i * window_overlap;
  }
  assert(windows_start.size() == (this->_v1.size() - window_size) / window_overlap);

  std::vector<std::pair<double, double>> correlations(windows_start.size());
  std::transform(std::execution::par_unseq, windows_start.begin(), windows_start.end(),
                 correlations.begin(), [&](uint32_t offset) {
                   auto v1s = this->_v1.begin() + offset;
                   auto v1e = v1s + window_size;
                   auto v2s = this->_v2.begin() + offset;
                   return compute_spearman(v1s, v1e, v2s);
                 });

  return correlations;
}

template <typename N>
std::vector<std::pair<double, double>> CorrelationTest<N>::compute_kendall(
    uint64_t window_size, uint64_t window_overlap) const {
  std::vector<uint32_t> windows_start((this->_v1.size() - window_size) / window_overlap);
  for (auto i = 0UL; i < windows_start.size(); ++i) {
    windows_start[i] = i * window_overlap;
  }
  assert(windows_start.size() == (this->_v1.size() - window_size) / window_overlap);
  std::vector<std::pair<double, double>> correlations(windows_start.size());
  std::transform(std::execution::par_unseq, windows_start.begin(), windows_start.end(),
                 correlations.begin(), [&](uint32_t offset) {
                   auto v1_s = this->_v1.begin() + offset;
                   auto v1_e = v1_s + window_size;
                   auto v2_s = this->_v2.begin() + offset;
                   return compute_kendall(v1_s, v1_e, v2_s);
                 });
  return correlations;
}

}  // namespace modle::correlation
