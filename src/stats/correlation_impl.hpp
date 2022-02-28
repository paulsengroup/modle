// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cpp-sort/sorters/pdq_sorter.h>  // for pdq_sorter

#include <algorithm>                                    // for min, all_of, max, clamp
#include <boost/math/distributions/beta.hpp>            // for cdf, beta_distribution
#include <boost/math/distributions/complement.hpp>      // for complement
#include <boost/math/distributions/students_t.hpp>      // for cdf, students_t_distribution
#include <boost/math/special_functions/fpclassify.hpp>  // for isinf, isnan
#include <cassert>                                      // for assert
#include <cmath>                                        // for sqrt, isnan
#include <iterator>                                     // for begin, end, distance
#include <limits>                                       // for numeric_limits
#include <numeric>
#include <type_traits>  // for is_arithmetic
#include <vector>       // for vector

#include "modle/common/common.hpp"                      // for u64
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_CONVERSION, DISAB...
#include "modle/common/utils.hpp"                       // for RepeatIterator
#include "modle/stats/descriptive.hpp"                  // for mean

namespace modle::stats {

template <class FP>
template <class It1, class It2, class It3>
typename Pearson<FP>::Result Pearson<FP>::operator()(It1 first1, It1 last1, It2 first2,
                                                     It3 weight_first) const {
  Result result;
  result.pcc = Pearson<FP>::compute_weighted_pcc(first1, last1, first2, weight_first);
  using weight_t = std::decay_t<decltype(*weight_first)>;
  if constexpr (std::is_same_v<utils::RepeatIterator<weight_t>, It3>) {
    result.pvalue = Pearson<FP>::compute_significance(result.pcc, std::distance(first1, last1));
  } else {
    result.pvalue = -1;
  }

  return result;
}

template <class FP>
template <class Range1, class Range2>
typename Pearson<FP>::Result Pearson<FP>::operator()(const Range1& r1, const Range2& r2) const {
  assert(std::distance(std::cbegin(r1), std::cend(r1)) ==
         std::distance(std::cbegin(r2), std::cend(r2)));
  return (*this)(std::cbegin(r1), std::cend(r1), std::cbegin(r2), utils::RepeatIterator<FP>(1));
}

template <class FP>
template <class Range1, class Range2, class Range3>
typename Pearson<FP>::Result Pearson<FP>::operator()(const Range1& r1, const Range2& r2,
                                                     const Range3& weights) const {
  assert(std::distance(std::cbegin(r1), std::cend(r1)) ==
         std::distance(std::cbegin(r2), std::cend(r2)));
  assert(std::distance(std::cbegin(r1), std::cend(r1)) ==
         std::distance(std::cbegin(weights), std::cend(weights)));

  return (*this)(std::cbegin(r1), std::cend(r1), std::cbegin(r2), std::cbegin(weights));
}

/*
 * J. Bennett, R. Grout, P. Pebay, D. Roe and D. Thompson, "Numerically stable, single-pass,
 * parallel statistics algorithms," 2009 IEEE International Conference on Cluster Computing and
 * Workshops, New Orleans, LA, 2009, pp. 1-8, doi: 10.1109/CLUSTR.2009.5289161.
 */

template <class FP>
template <class It1, class It2, class It3>
FP Pearson<FP>::compute_weighted_pcc(It1 first1, It1 last1, It2 first2, It3 weight_first) {
  const auto size = std::distance(first1, last1);

  const auto wmean1 = stats::weighted_mean(first1, last1, weight_first);
  const auto wmean2 = stats::weighted_mean(first2, first2 + size, weight_first);

  const auto cov1 = stats::weighted_covariance(first1, last1, first1, weight_first, wmean1, wmean1);
  const auto cov2 =
      stats::weighted_covariance(first2, first2 + size, first2, weight_first, wmean2, wmean2);
  const auto cov3 = stats::weighted_covariance(first1, last1, first2, weight_first, wmean1, wmean2);

  const auto pcc = std::clamp(cov3 / std::sqrt(cov1 * cov2), -1.0, 1.0);

  return pcc;
}

template <class FP>
template <class It1, class It2>
FP Pearson<FP>::compute_pcc(It1 first1, It1 last1, It2 first2) {
  return compute_weighted_pcc(first1, last1, first2, utils::RepeatIterator<FP>(1));
}

template <class FP>
template <class I, class>
FP Pearson<FP>::compute_significance(const FP pcc, const I n) {
  if (MODLE_UNLIKELY(std::isnan(pcc))) {
    return std::numeric_limits<FP>::quiet_NaN();
  }
  assert(n > I(2));
  const auto ab = (static_cast<FP>(n) / FP(2)) - 1;
  boost::math::beta_distribution<FP> dist(ab, ab);
  // https://github.com/scipy/scipy/blob/6703631bcd15750e86f4098b0421efabcac0f7c2/scipy/stats/stats.py#L3885
  return FP(2) * boost::math::cdf<FP>(dist, FP(0.5) * (FP(1) - std::abs(pcc)));
}

template <class FP>
Spearman<FP>::Spearman(const usize n) : _rank_buff1(n), _rank_buff2(n), _idx_buff(n) {}

template <class FP>
template <class It1, class It2>
FP Spearman<FP>::compute_rho(It1 first1, It1 last1, It2 first2) {
  assert(std::distance(first1, last1) >= 0);

  const auto size = static_cast<usize>(std::distance(first1, last1));
  auto last2 = first2 + static_cast<isize>(size);
  if (MODLE_UNLIKELY(std::all_of(first1, last1, [](const auto n) { return n == FP(0); }) ||
                     std::all_of(first2, last2, [](const auto n) { return n == FP(0); }))) {
    return FP(1);
  }

  this->_rank_buff1.resize(size);
  this->_rank_buff2.resize(size);
  internal::compute_element_ranks(first1, last1, this->_rank_buff1, this->_idx_buff);
  internal::compute_element_ranks(first2, last2, this->_rank_buff2, this->_idx_buff);
  return Pearson<FP>::compute_pcc(this->_rank_buff1.begin(), this->_rank_buff1.end(),
                                  this->_rank_buff2.begin());
}

template <class FP>
template <class It1, class It2, class It3>
FP Spearman<FP>::compute_weighted_rho(It1 first1, It1 last1, It2 first2, It3 weight_first) {
  assert(std::distance(first1, last1) >= 0);

  const auto size = static_cast<usize>(std::distance(first1, last1));
  auto last2 = first2 + static_cast<isize>(size);
  if (MODLE_UNLIKELY(std::all_of(first1, last1, [](const auto n) { return n == FP(0); }) ||
                     std::all_of(first2, last2, [](const auto n) { return n == FP(0); }))) {
    return FP(1);
  }

  this->_rank_buff1.resize(size);
  this->_rank_buff2.resize(size);
  internal::compute_weighted_element_ranks(first1, last1, weight_first, this->_rank_buff1,
                                           this->_idx_buff);
  internal::compute_weighted_element_ranks(first2, last2, weight_first, this->_rank_buff2,
                                           this->_idx_buff);
  return Pearson<FP>::compute_weighted_pcc(this->_rank_buff1.begin(), this->_rank_buff1.end(),
                                           this->_rank_buff2.begin(), weight_first);
}

template <class FP>
template <class I, class>
FP Spearman<FP>::compute_significance(const FP rho, const I n) {
  assert(n > I(2));
  if (MODLE_UNLIKELY(std::isnan(rho))) {
    return std::numeric_limits<FP>::quiet_NaN();
  }
  const auto dof = static_cast<FP>(n - 2);
  const double tscore = rho * std::sqrt(dof / ((FP(1) + rho) * (FP(1) - rho)));
  boost::math::students_t_distribution<FP> dist(dof);
  // https://github.com/scipy/scipy/blob/6703631bcd15750e86f4098b0421efabcac0f7c2/scipy/stats/stats.py#L4229
  // 2 * survival function
  return FP(2) * boost::math::cdf<FP>(boost::math::complement(dist, std::fabs(tscore)));
}

template <class FP>
template <class It1, class It2, class It3>
typename Spearman<FP>::Result Spearman<FP>::operator()(It1 first1, It1 last1, It2 first2,
                                                       It3 weight_first) {
  Result result;
  using weight_t = std::decay_t<decltype(*weight_first)>;
  if constexpr (std::is_same_v<utils::RepeatIterator<weight_t>, It3>) {
    result.rho = this->compute_rho(first1, last1, first2);
    result.pvalue = Spearman<FP>::compute_significance(result.rho, std::distance(first1, last1));
  } else {
    result.rho = this->compute_weighted_rho(first1, last1, first2, weight_first);
    result.pvalue = -1;
  }
  return result;
}

template <class FP>
template <class Range1, class Range2>
typename Spearman<FP>::Result Spearman<FP>::operator()(const Range1& r1, const Range2& r2) {
  assert(std::distance(std::cbegin(r1), std::cend(r1)) ==
         std::distance(std::cbegin(r2), std::cend(r2)));
  return (*this)(std::cbegin(r1), std::cend(r1), std::cbegin(r2), utils::RepeatIterator<FP>(1));
}

template <class FP>
template <class Range1, class Range2, class Range3>
typename Spearman<FP>::Result Spearman<FP>::operator()(const Range1& r1, const Range2& r2,
                                                       const Range3& weights) {
  assert(std::distance(std::cbegin(r1), std::cend(r1)) ==
         std::distance(std::cbegin(r2), std::cend(r2)));
  assert(std::distance(std::cbegin(r1), std::cend(r1)) ==
         std::distance(std::cbegin(weights), std::cend(weights)));
  return (*this)(std::cbegin(r1), std::cend(r1), std::cbegin(r2), std::cbegin(weights));
}

namespace internal {

template <class It>
void sort_by_index(It first, It last, std::vector<usize>& idx_buff) {
  assert(std::distance(first, last) >= 0);
  idx_buff.resize(static_cast<usize>(std::distance(first, last)));
  std::iota(idx_buff.begin(), idx_buff.end(), 0);

  cppsort::pdq_sort(idx_buff.begin(), idx_buff.end(), [&](auto i1, auto i2) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_SIGN_CONVERSION
    return *(first + i1) < *(first + i2);
    DISABLE_WARNING_POP
  });
}

template <class It, class FP, class>
void compute_element_ranks(It first, It last, std::vector<FP>& rank_buff,
                           std::vector<usize>& idx_buff) {
  assert(std::distance(first, last) >= 0);
  const auto size = static_cast<usize>(std::distance(first, last));
  rank_buff.resize(size);
  idx_buff.resize(size);
  sort_by_index(first, last, idx_buff);

  rank_buff[idx_buff.front()] = FP(0);
  for (usize i = 1; i < size; ++i) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_SIGN_CONVERSION
    const auto& prev = *(first + idx_buff[i - 1]);
    const auto& current = *(first + idx_buff[i]);
    DISABLE_WARNING_POP

    if (MODLE_UNLIKELY(prev == current)) {   // Process ties
      const usize lower_bound_idx = i - 1;   // first tied-element
      const usize upper_bound_idx = [&]() {  // last tied-element
        for (auto j = lower_bound_idx + 1; j < size; ++j) {
          const auto& next = *(first + static_cast<isize>(idx_buff[j]));
          if (current != next) {
            return j;
          }
        }
        return size;
      }();
      const auto avg_rank = static_cast<FP>(lower_bound_idx + upper_bound_idx - 1) / FP(2);
      // Assign the average rank to all tied elements
      for (i = lower_bound_idx; i < upper_bound_idx; ++i) {
        rank_buff[idx_buff[i]] = avg_rank;
      }
      // Deal with the remote possibility that the last element in the collection is tied
      if (MODLE_UNLIKELY(i == size)) {
        return;
      }
    }
    rank_buff[idx_buff[i]] = static_cast<FP>(i);
  }
}

// https://rdrr.io/cran/wCorr/f/inst/doc/wCorrFormulas.pdf
template <class It1, class It2, class FP, class>
void compute_weighted_element_ranks(It1 first, It1 last, It2 weight_first,
                                    std::vector<FP>& rank_buff, std::vector<usize>& idx_buff) {
  assert(std::distance(first, last) >= 0);
  const auto size = static_cast<usize>(std::distance(first, last));
  rank_buff.resize(size);
  idx_buff.resize(size);
  sort_by_index(first, last, idx_buff);

  double weight_sum = *weight_first;
  rank_buff[idx_buff.front()] = FP(weight_sum);
  for (usize i = 1; i < size; ++i) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_SIGN_CONVERSION
    const auto& prev = *(first + idx_buff[i - 1]);
    const auto& current = *(first + idx_buff[i]);
    const auto& prev_weight = *(weight_first + idx_buff[i - 1]);
    DISABLE_WARNING_POP

    if (MODLE_UNLIKELY(prev == current)) {   // Process ties
      const usize lower_bound_idx = i - 1;   // first tied-element
      const usize upper_bound_idx = [&]() {  // last tied-element
        for (auto j = lower_bound_idx + 1; j < size; ++j) {
          const auto& next = *(first + static_cast<isize>(idx_buff[j]));
          if (current != next) {
            return j;
          }
        }
        return size;
      }();

      // Subtract the previous weight now that we know that the i-1 element belongs to a tie
      weight_sum -= prev_weight;
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      const auto weight_sum_ties = std::accumulate(
          idx_buff.begin() + lower_bound_idx, idx_buff.begin() + upper_bound_idx, FP(0),
          [&](auto accumulator, auto ii) { return accumulator + FP(*(weight_first + ii)); });
      DISABLE_WARNING_POP

      const auto num_ties = upper_bound_idx - lower_bound_idx;
      const auto avg_tie_weight = weight_sum_ties / static_cast<FP>(num_ties);

      // weight_sum == aj, the rest corresponds to the bj term in the eq from wCorrFormulas.pdf
      const auto avg_rank = weight_sum + ((static_cast<FP>(num_ties + 1) / FP(2)) * avg_tie_weight);
      // Assign the average rank to all tied elements
      for (i = lower_bound_idx; i < upper_bound_idx; ++i) {
        rank_buff[idx_buff[i]] = avg_rank;
      }

      // This is equivalent to adding n * avg_tie_weight
      weight_sum += weight_sum_ties;
      // Deal with the remote possibility that the last element in the collection is tied
      if (MODLE_UNLIKELY(i == size)) {
        return;
      }
    }
    // Here aj = weight_sum + current_weight; bj = 0;
    const auto& current_weight = *(weight_first + static_cast<isize>(idx_buff[i]));
    rank_buff[idx_buff[i]] = (weight_sum += current_weight);
  }
}

}  // namespace internal

}  // namespace modle::stats

// IWYU pragma: private, include "modle/correlation.hpp"
// IWYU pragma: no_include <boost/exception/exception.hpp>
// IWYU pragma: no_include <boost/math/policies/policy.hpp>
