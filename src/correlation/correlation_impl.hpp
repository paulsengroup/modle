#pragma once

#include <algorithm>
#include <atomic>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

#include "./correlation_utils.hpp"
#include "modle/correlation.hpp"
#include "range/v3/algorithm.hpp"
#include "range/v3/range.hpp"

namespace modle::correlation {

double compute_pearson_significance(double pcc, std::size_t n) {
  assert(n > 2);
  const auto ab = static_cast<const double>(n) / 2.0 - 1;
  boost::math::beta_distribution<double> dist(ab, ab);
  // https://github.com/scipy/scipy/blob/6703631bcd15750e86f4098b0421efabcac0f7c2/scipy/stats/stats.py#L3885
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  return 2.0 * boost::math::cdf<double>(dist, 0.5 * (1 - std::abs(pcc)));
}

double compute_spearman_significance(double rho, std::size_t n) {
  assert(n > 2);
  const auto dof = static_cast<const double>(n - 2);
  const double tscore = rho * std::sqrt(dof / ((1.0 + rho) * (1.0 - rho)));
  boost::math::students_t_distribution<double> dist(dof);
  // https://github.com/scipy/scipy/blob/6703631bcd15750e86f4098b0421efabcac0f7c2/scipy/stats/stats.py#L4229
  // 2 * survival function
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  return 2.0 * boost::math::cdf<double>(boost::math::complement(dist, std::fabs<double>(tscore)));
}

/*
 * J. Bennett, R. Grout, P. Pebay, D. Roe and D. Thompson, "Numerically stable, single-pass,
 * parallel statistics algorithms," 2009 IEEE International Conference on Cluster Computing and
 * Workshops, New Orleans, LA, 2009, pp. 1-8, doi: 10.1109/CLUSTR.2009.5289161.
 */
template <typename Rng>
double compute_pearson(Rng r1, Rng r2) {
  static_assert(ranges::random_access_range<Rng>, "r1 and r2 should be a random access range");
  double cov = 0;
  double r1_avg = r1[0];
  double r2_avg = r2[0];
  double d1 = 0;
  double d2 = 0;
  for (std::size_t i = 1; i < r1.size(); ++i) {
    auto r1_tmp = r1[i] - r1_avg;
    auto r2_tmp = r2[i] - r2_avg;
    d1 += (i * r1_tmp * r1_tmp) / (i + 1);
    d2 += (i * r2_tmp * r2_tmp) / (i + 1);
    cov += i * r1_tmp * r2_tmp / (i + 1);
    r1_avg += r1_tmp / (i + 1);
    r2_avg += r2_tmp / (i + 1);
  }

  // Both datasets are constant (i.e. perfectly correlated)
  if (d1 == 0 && d2 == 0) {
    return 1.0;
  }
  // Only one of the two datasets is constant (i.e. there's no correlation)
  if (d1 == 0 || d2 == 0) {
    return 0.0;
  }

  const auto pcc = std::clamp(cov / std::sqrt(d1 * d2), -1.0, 1.0);
  if (std::isnan(pcc)) {
    throw std::logic_error("compute_pearson: pcc cannot be nan!");
  }

  return pcc;
};

template <typename Rng>
double compute_spearman(Rng r1, Rng r2) {
  static_assert(ranges::random_access_range<Rng>, "r1 and r2 should be a random access range");
  // Shortcut to avoid computing the correlation when any of the vector is all zeros
  if (ranges::all_of(r1, [](auto n) { return n == 0; }) ||
      ranges::all_of(r2, [](auto n) { return n == 0; })) {
    return 1.0;
  }
  return compute_pearson(utils::compute_element_ranks(r1), utils::compute_element_ranks(r2));
}

template <typename N>
std::pair<std::vector<double>, std::vector<double>> compute_pearson(const std::vector<N>& v1,
                                                                    const std::vector<N>& v2,
                                                                    std::size_t window_span,
                                                                    std::size_t window_overlap) {
  // TODO: compute_pearson and compute_spearman (implemented below this function), basically do the
  // same thing. Figure out a way to remove redundant code (i.e. everything excepr compute_* and
  // compute_*_significance
  static_assert(std::is_arithmetic<N>::value,
                "compute_pearson requires a numeric type as template argument.");
  if (v1.size() != v2.size()) {
    throw std::runtime_error(
        fmt::format("compute_pearson expects a pair of vectors of the same size, got {} "
                    "and {} respectively",
                    v1.size(), v2.size()));
  }
  assert(window_span > window_overlap);
  std::vector<double> pcc_vals(v1.size() / (window_span - window_overlap));
  std::vector<double> p_vals(pcc_vals.size());

  for (std::size_t i = 0; i < pcc_vals.size(); ++i) {
    const auto window_start = i * (window_span - window_overlap);
    const auto window_end = window_start + window_span;
    const auto slice1 = v1 | ranges::views::slice(window_start, window_end);
    const auto slice2 = v2 | ranges::views::slice(window_start, window_end);
    pcc_vals[i] = compute_pearson(slice1, slice2);
    p_vals[i] = compute_pearson_significance(pcc_vals[i], slice1.size());
  }

  return std::make_pair(pcc_vals, p_vals);
}

template <typename N>
std::pair<std::vector<double>, std::vector<double>> compute_spearman(const std::vector<N>& v1,
                                                                     const std::vector<N>& v2,
                                                                     std::size_t window_span,
                                                                     std::size_t window_overlap) {
  static_assert(std::is_arithmetic<N>::value,
                "compute_spearman requires a numeric type as template argument.");
  if (v1.size() != v2.size()) {
    throw std::runtime_error(
        fmt::format("compute_spearman expects a pair of vectors of the same size, got {} "
                    "and {} respectively",
                    v1.size(), v2.size()));
  }
  assert(window_span > window_overlap);
  std::vector<double> rho_vals(v1.size() / (window_span - window_overlap));
  std::vector<double> p_vals(rho_vals.size());

  for (std::size_t i = 0; i < rho_vals.size(); ++i) {
    const auto window_start = i * (window_span - window_overlap);
    const auto window_end = window_start + window_span;
    const auto slice1 = v1 | ranges::views::slice(window_start, window_end);
    const auto slice2 = v2 | ranges::views::slice(window_start, window_end);
    rho_vals[i] = compute_spearman(slice1, slice2);
    p_vals[i] = compute_spearman_significance(rho_vals[i], slice1.size());
  }

  return std::make_pair(rho_vals, p_vals);
}
/*
template <typename N>
std::pair<std::vector<double>, std::vector<double>> compute_pearson_L(const std::vector<N>& v1,
                                                                      const std::vector<N>& v2,
                                                                      std::size_t nrows,
                                                                      std::size_t ncols) {
  // TODO: compute_pearson and compute_spearman (implemented below this function), basically do the
  // same thing. Figure out a way to remove redundant code (i.e. everything excepr compute_* and
  // compute_*_significance
  static_assert(std::is_arithmetic<N>::value,
                "compute_pearson requires a numeric type as template argument.");
  if (v1.size() != v2.size()) {
    throw std::runtime_error(
        fmt::format("compute_pearson expects a pair of vectors of the same size, got {} "
                    "and {} respectively",
                    v1.size(), v2.size()));
  }
  assert(ncols > nrows);
  std::vector<double> pcc_vals(ncols);
  std::vector<double> p_vals(pcc_vals.size());

  std::vector<N> vl1((2 * nrows) - 1);
  std::vector<N> vl2(vl1.size());
  for (std::size_t i = 0; i < pcc_vals.size(); ++i) {
    std::size_t idx = i * nrows;
    for (std::size_t j = 0; j < nrows; ++j) {
      if (idx >= v1.size()) {
        std::fill(vl1.begin() + j, vl1.begin() + nrows, 0);
        std::fill(vl2.begin() + j, vl2.begin() + nrows, 0);
        break;
      }
      vl1[j] = v1[idx];
      vl2[j] = v2[idx];
      idx += nrows + 1;
    }
    idx = i * nrows - 1;
    for (std::size_t j = nrows; j < 2 * nrows; ++j) {
      if (idx > (i * nrows) + i) {
        std::fill(vl1.begin() + j, vl1.end(), 0);
        std::fill(vl2.begin() + j, vl2.end(), 0);
        break;
      }
      vl1[j] = v1[idx];
      vl2[j] = v2[idx];
      ++idx;
    }

    // fmt::print(stderr, "{} - {}\n{} - {}\n", i, absl::StrJoin(vl1, ", "), i, absl::StrJoin(vl2,
    // ", ")); usleep(250000);

    // assert(idx < v1.size());
    pcc_vals[i] = compute_pearson(vl1, vl2);
    p_vals[i] = compute_pearson_significance(pcc_vals[i], vl1.size());
  }

  return std::make_pair(pcc_vals, p_vals);
}
*/

/*
template <typename N>
double compute_kendall_b_significance(uint32_t nc, uint32_t nd, uint32_t size,
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
  const auto x0 = static_cast<double>(tie1 * (tie1 - 1) * (tie1 - 2));
  const auto y0 = static_cast<double>(tie2 * (tie2 - 1) * (tie2 - 2));
  const auto x1 = static_cast<double>(tie1 * (tie1 - 1) * (2 * tie1 + 5));
  const auto y1 = static_cast<double>(tie2 * (tie2 - 1) * (2 * tie2 + 5));
  static constexpr auto sqrt2 = boost::math::constants::root_two<double>();
  const double var = (size * (size - 1) * (2.0 * size + 5) - x1 - y1) / 18.0 +
                     (2.0 * tie1 * tie2) / (size * (size - 1)) +
                     x0 * y0 / (9.0 * size * (size - 1) * (size - 2));
  return boost::math::erfc<double>(std::abs<double>(static_cast<double>(nc) - nd) / std::sqrt(var) /
                                   sqrt2);
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
 */

}  // namespace modle::correlation
