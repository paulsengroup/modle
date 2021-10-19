// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cpp-sort/sorters/pdq_sorter.h>  // for pdq_sorter
#include <cpp-sort/sorters/ska_sorter.h>  // for ska_sort

#include <algorithm>                                    // for min, all_of, max, clamp
#include <boost/math/distributions/beta.hpp>            // for cdf, beta_distribution
#include <boost/math/distributions/complement.hpp>      // for complement
#include <boost/math/distributions/students_t.hpp>      // for cdf, students_t_distribution
#include <boost/math/special_functions/fpclassify.hpp>  // for isinf, isnan
#include <cassert>                                      // for assert
#include <cmath>                                        // for sqrt, isnan
#include <cstdlib>                                      // for abs
#include <numeric>
#include <stdexcept>    // for overflow_error, logic_error, runt...
#include <type_traits>  // for is_arithmetic
#include <utility>      // for make_pair, pair
#include <vector>       // for vector

#include "modle/common/common.hpp"                      // for u64
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_CONVERSION, DISAB...
#include "modle/common/utils.hpp"                       // for

namespace modle::correlation {

template <class FP>
template <class It1, class It2>
typename Pearson<FP>::Result Pearson<FP>::operator()(It1 begin1, It1 end1, It2 begin2) const {
  Result result;
  result.pcc = Pearson<FP>::compute_pcc(begin1, end1, begin2);
  result.pvalue = Pearson<FP>::compute_significance(result.pcc, std::distance(begin1, end1));

  return result;
}

template <class FP>
template <class Range1, class Range2>
typename Pearson<FP>::Result Pearson<FP>::operator()(const Range1& r1, const Range2& r2) const {
  assert(std::distance(std::cbegin(r1), std::cend(r1)) ==
         std::distance(std::cbegin(r2), std::cend(r2)));  // NOLINT
  return this->operator()(std::cbegin(r1), std::cend(r1), std::cbegin(r2));
}

/*
 * J. Bennett, R. Grout, P. Pebay, D. Roe and D. Thompson, "Numerically stable, single-pass,
 * parallel statistics algorithms," 2009 IEEE International Conference on Cluster Computing and
 * Workshops, New Orleans, LA, 2009, pp. 1-8, doi: 10.1109/CLUSTR.2009.5289161.
 */

template <class FP>
template <class It1, class It2>
FP Pearson<FP>::compute_pcc(It1 begin1, It1 end1, It2 begin2) {
  const auto size = std::distance(begin1, end1);
  assert(size >= 0);  // NOLINT

  FP cov = 0;
  FP r1_avg = size == 0 ? 0 : *begin1;
  FP r2_avg = size == 0 ? 0 : *begin2;
  FP d1 = 0;
  FP d2 = 0;

  usize i = 1;
  for (auto it1 = begin1 + 1, it2 = begin2 + 1; it1 != end1; ++it1, ++it2, ++i) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_USELESS_CAST
    const auto ii = static_cast<FP>(i);
    const auto r1 = static_cast<FP>(*it1) - r1_avg;
    const auto r2 = static_cast<FP>(*it2) - r2_avg;
    DISABLE_WARNING_POP

    d1 += (ii * r1 * r1) / (ii + 1);
    d2 += (ii * r2 * r2) / (ii + 1);
    cov += ii * r1 * r2 / (ii + 1);
    r1_avg += r1 / (ii + 1);
    r2_avg += r2 / (ii + 1);
  }

  // Both datasets are constant (i.e. perfectly correlated)
  if (d1 == 0 && d2 == 0) {
    return FP(1);
  }
  // Only one of the two datasets is constant (i.e. there's no correlation)
  if (d1 == 0 || d2 == 0) {
    return FP(0);
  }

  const auto pcc = std::clamp(cov / std::sqrt(d1 * d2), FP(-1), FP(1));
  assert(!std::isnan(pcc));  // NOLINT

  return pcc;
}

template <class FP>
template <class I, class>
FP Pearson<FP>::compute_significance(const FP pcc, const I n) {
  assert(n > I(2));  // NOLINT
  const auto ab = (static_cast<FP>(n) / FP(2)) - 1;
  boost::math::beta_distribution<FP> dist(ab, ab);
  // https://github.com/scipy/scipy/blob/6703631bcd15750e86f4098b0421efabcac0f7c2/scipy/stats/stats.py#L3885
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  return FP(2) * boost::math::cdf<FP>(dist, FP(0.5) * (FP(1) - std::abs(pcc)));
}

template <class FP>
Spearman<FP>::Spearman(const usize n) : _rank_buff1(n), _rank_buff2(n), _idx_buff(n) {}

template <class FP>
template <class It1, class It2>
FP Spearman<FP>::compute_rho(It1 begin1, It1 end1, It2 begin2) {
  assert(std::distance(begin1, end1) >= 0);  // NOLINT

  const auto size = static_cast<usize>(std::distance(begin1, end1));
  auto end2 = begin2 + static_cast<isize>(size);
  if (BOOST_UNLIKELY(std::all_of(begin1, end1, [](const auto n) { return n == FP(0); }) ||
                     std::all_of(begin2, end2, [](const auto n) { return n == FP(0); }))) {
    return FP(1);
  }

  this->_rank_buff1.resize(size);
  this->_rank_buff2.resize(size);
  Spearman<FP>::compute_element_ranks(begin1, end1, this->_rank_buff1, this->_idx_buff);
  Spearman<FP>::compute_element_ranks(begin2, end2, this->_rank_buff2, this->_idx_buff);
  return Pearson<FP>::compute_pcc(this->_rank_buff1.cbegin(), this->_rank_buff1.cend(),
                                  this->_rank_buff2.cbegin());
}

template <class FP>
template <class I, class>
FP Spearman<FP>::compute_significance(const FP rho, const I n) {
  assert(n > I(2));  // NOLINT
  const auto dof = static_cast<FP>(n - 2);
  const double tscore = rho * std::sqrt(dof / ((FP(1) + rho) * (FP(1) - rho)));
  boost::math::students_t_distribution<FP> dist(dof);
  // https://github.com/scipy/scipy/blob/6703631bcd15750e86f4098b0421efabcac0f7c2/scipy/stats/stats.py#L4229
  // 2 * survival function
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  return FP(2) * boost::math::cdf<FP>(boost::math::complement(dist, std::fabs(tscore)));
}

template <class FP>
template <class It>
void Spearman<FP>::sort_by_index(It it_begin, It it_end, std::vector<usize>& idx_buff) {
  assert(std::distance(it_begin, it_end) >= 0);  // NOLINT
  idx_buff.resize(static_cast<usize>(std::distance(it_begin, it_end)));
  std::iota(idx_buff.begin(), idx_buff.end(), 0);

  if constexpr (std::is_arithmetic_v<decltype(*it_begin)>) {
    cppsort::ska_sort(it_begin, it_end);
  } else {
    cppsort::pdq_sort(idx_buff.begin(), idx_buff.end(), [&](auto i1, auto i2) {
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      return *(it_begin + i1) < *(it_begin + i2);
      DISABLE_WARNING_POP
    });
  }
}

// This returns a vector of integers corresponding to the rank of the values of
// vector v
template <class FP>
template <class It>
void Spearman<FP>::compute_element_ranks(It it_begin, It it_end, std::vector<FP>& rank_buff,
                                         std::vector<usize>& idx_buff) {
  assert(std::distance(it_begin, it_end) >= 0);  // NOLINT
  const auto size = static_cast<usize>(std::distance(it_begin, it_end));
  rank_buff.resize(size);
  sort_by_index(it_begin, it_end, idx_buff);

  rank_buff[idx_buff[0]] = FP(0);
  for (usize i = 1; i < size; ++i) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_SIGN_CONVERSION
    DISABLE_WARNING_USELESS_CAST
    auto& prev = *(it_begin + idx_buff[i - 1]);
    auto& current = *(it_begin + idx_buff[i]);
    if (prev != current) {
      rank_buff[idx_buff[i]] = static_cast<FP>(i);
    } else {
      const auto min = i - 1;
      auto max = i;
      for (auto j = max + 1; j < size && current == *(it_begin + idx_buff[j]); ++j) {
        max = j;
      }
      for (auto j = min; j <= max; ++j) {
        // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
        rank_buff[idx_buff[j]] = static_cast<FP>(max + min) / FP(2);
      }
      i = max;
    }
    DISABLE_WARNING_POP
  }
}

template <class FP>
template <class It1, class It2>
typename Spearman<FP>::Result Spearman<FP>::operator()(It1 begin1, It1 end1, It2 begin2) {
  Result result;
  result.rho = this->compute_rho(begin1, end1, begin2);
  result.pvalue = Spearman<FP>::compute_significance(result.rho, std::distance(begin1, end1));

  return result;
}

template <class FP>
template <class Range1, class Range2>
typename Spearman<FP>::Result Spearman<FP>::operator()(const Range1& r1, const Range2& r2) {
  assert(std::distance(std::cbegin(r1), std::cend(r1)) ==
         std::distance(std::cbegin(r2), std::cend(r2)));  // NOLINT
  return this->operator()(std::cbegin(r1), std::cend(r1), std::cbegin(r2));
}

template <class FP>
template <class It1, class It2>
inline FP SED<FP>::operator()(It1 begin1, It1 end1, It2 begin2) const {
  const isize size = std::distance(begin1, end1);
  assert(size >= 0);  // NOLINT

  if constexpr (std::is_integral_v<decltype(*begin1)>) {
    u64 sed = 0;
    for (isize i = 0; i < size; ++i) {
      const auto n1 = *(begin1 + i);
      const auto n2 = *(begin2 + i);
      const auto n = n1 > n2 ? n1 - n2 : n2 - n1;
      sed += n * n;
    }
    return static_cast<FP>(sed);
  } else {
    FP sed = 0;
    for (isize i = 0; i < size; ++i) {
      const auto n1 = *(begin1 + i);
      const auto n2 = *(begin2 + i);
      const auto n = n1 > n2 ? n1 - n2 : n2 - n1;
      sed += n * n;
    }
    return sed;
  }
}

template <class FP>
template <class Range1, class Range2>
[[nodiscard]] inline FP SED<FP>::operator()(const Range1& r1, const Range2& r2) const {
  assert(std::distance(std::cbegin(r1), std::cend(r1)) ==
         std::distance(std::cbegin(r2), std::cend(r2)));  // NOLINT
  return this->operator()(std::cbegin(r1), std::cend(r1), std::cbegin(r2));
}

/*
template <typename N1, typename N2>
std::pair<std::vector<double>, std::vector<double>> compute_corr(absl::Span<const N1> v1,
                                                                 absl::Span<const N2> v2,
                                                                 Algorithm type, usize window_span,
                                                                 usize window_overlap) {
  if constexpr (utils::ndebug_not_defined()) {
    if (v1.size() != v2.size()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("compute_pearson expects a pair of vectors of the same size, got {} "
                     "and {} respectively"),
          v1.size(), v2.size()));
    }
  }
  assert(window_span > window_overlap);
  std::vector<double> r_vals(v1.size() / (window_span - window_overlap));
  std::vector<double> p_vals(r_vals.size());

  double (*compute_corr_fp)(absl::Span<const N1>, absl::Span<const N2>) = nullptr;
  double (*compute_corr_significance_fp)(double, usize) = nullptr;
  if (type == Algorithm::pearson) {
    compute_corr_fp = &compute_pearson;
    compute_corr_significance_fp = (&compute_pearson_significance);
  } else if (type == Algorithm::spearman) {
    compute_corr_fp = &compute_spearman;
    compute_corr_significance_fp = (&compute_spearman_significance);
  }
  assert(compute_corr_fp);
  assert(compute_corr_significance_fp);

  for (usize i = 0; i < r_vals.size(); ++i) {
    const auto window_start = i * (window_span - window_overlap);
    const auto window_end = window_start + window_span;
    auto slice1 = absl::MakeConstSpan(v1).subspan(window_start, window_end);
    auto slice2 = absl::MakeConstSpan(v2).subspan(window_start, window_end);
    r_vals[i] = compute_corr_fp(slice1, slice2);
    p_vals[i] = compute_corr_significance_fp(r_vals[i], slice1.size());
  }

  return std::make_pair(r_vals, p_vals);
}

template <typename N1, typename N2>
double compute_sed(absl::Span<const N1> v1, absl::Span<const N2> v2) {
  static_assert(std::is_arithmetic<N1>::value,
                "v1 should be convertible to a Span of numeric type");
  static_assert(std::is_arithmetic<N2>::value,
                "v2 should be convertible to a Span of numeric type");

  u64 sed{0};
  for (usize i = 0; i < v1.size(); ++i) {
    const auto n = v1[i] > v2[i] ? v1[i] - v2[i] : v2[i] - v1[i];
    sed += n * n;
  }
  return static_cast<double>(sed);
}

template <typename N1, typename N2>
double compute_sed(const std::vector<N1>& v1, const std::vector<N2>& v2) {
  return compute_sed(absl::MakeConstSpan(v1), absl::MakeConstSpan(v2));
}
 */
}  // namespace modle::correlation

// IWYU pragma: private, include "modle/correlation.hpp"
// IWYU pragma: no_include <boost/exception/exception.hpp>
// IWYU pragma: no_include <boost/math/policies/policy.hpp>
