// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <iterator>  // for distance
#include <numeric>   // for accumulate
#include <vector>    // for vector

#include "modle/common/common.hpp"                      // for usize, isize
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...
#include "modle/common/utils.hpp"                       // for RepeatIterator::operator*, Repeat...

namespace modle::stats {
template <class T, class>
constexpr T abs_diff(const T a, const T b) noexcept {
  return a > b ? a - b : b - a;
}

template <class InputIt, class UnaryOperation, class>
double mean(InputIt first, InputIt last, UnaryOperation op) {
  if (first == last) {
    return 0.0;
  }
  return std::accumulate(
             first, last, 0.0,
             [&](auto accumulator, const auto& val) { return accumulator + double(op(val)); }) /
         static_cast<double>(std::distance(first, last));
}

template <class InputIt1, class InputIt2, class UnaryOperation, class>
double weighted_mean(InputIt1 first, InputIt1 last, InputIt2 weight_first, UnaryOperation op) {
  double total = 0;
  double weight_total = 0;

  for (; first != last; ++first, ++weight_first) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_USELESS_CAST
    total += static_cast<double>(op(*first)) * static_cast<double>(*weight_first);
    weight_total += static_cast<double>(*weight_first);
    DISABLE_WARNING_POP
  }

  return total / weight_total;
}

template <class InputIt, class OutputIt, class I, class UnaryOperation, class>
usize moving_average(InputIt first, InputIt last, OutputIt output_begin, const I window_size,
                     UnaryOperation op) {
  using OutputItValueT = typename OutputIt::value_type;
  if (static_cast<isize>(window_size) >= std::distance(first, last)) {
    *output_begin = OutputItValueT(stats::mean(first, last, op));
    return 1;
  }

  for (auto it1 = first, it2 = first + static_cast<isize>(window_size); it2 != last;
       ++it1, ++it2, ++output_begin) {
    *output_begin = OutputItValueT(mean(it1, it2, op));
  }
  return usize(std::distance(first, last)) - usize(window_size);
}

template <class InputIt, class FP, class UnaryOperation, class>
double sum_of_squared_deviations(InputIt first, InputIt last, const FP sample_mean,
                                 UnaryOperation op) {
  if (first == last) {
    return 0.0;
  }
  return std::accumulate(first, last, 0.0,
                         [&, sm = double(sample_mean)](const auto accumulator, const auto& val) {
                           const auto n = double(op(val)) - sm;
                           return accumulator + (n * n);
                         });
}

template <class InputIt, class UnaryOperation, class>
double sum_of_squared_deviations(InputIt first, InputIt last, UnaryOperation op) {
  return stats::sum_of_squared_deviations(first, last, stats::mean(first, last, op), op);
}

template <class InputIt, class FP, class UnaryOperation, class>
double variance(InputIt first, InputIt last, const FP sample_mean, UnaryOperation op) {
  if (first == last) {
    return 0.0;
  }
  return stats::sum_of_squared_deviations(first, last, sample_mean, op) /
         static_cast<double>(std::distance(first, last));
}

template <class InputIt, class UnaryOperation, class>
double variance(InputIt first, InputIt last, UnaryOperation op) {
  return stats::variance(first, last, stats::mean(first, last, op), op);
}

template <class InputIt, class FP, class UnaryOperation, class>
double standard_dev(InputIt first, InputIt last, const FP sample_mean, UnaryOperation op) {
  if (first == last) {
    return 0.0;
  }
  return std::sqrt(stats::variance(first, last, sample_mean, op));
}

template <class InputIt, class UnaryOperation, class>
double standard_dev(InputIt first, InputIt last, UnaryOperation op) {
  return stats::standard_dev(first, last, stats::mean(first, last, op), op);
}

template <class InputIt, class FP, class UnaryOperation, class>
double covariance(InputIt first1, InputIt last1, InputIt first2, FP mean1, FP mean2,
                  UnaryOperation op) {
  double sum_of_squared_devs = 0;
  usize i = 0;
  for (; first1 != last1; ++first1, ++first2, ++i) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_USELESS_CAST
    const auto n1 = static_cast<double>(op(*first1)) - static_cast<double>(mean1);
    const auto n2 = static_cast<double>(op(*first2)) - static_cast<double>(mean2);
    DISABLE_WARNING_POP
    sum_of_squared_devs += n1 * n2;
  }

  return sum_of_squared_devs / static_cast<double>(i);
}

template <class InputIt, class UnaryOperation, class>
double covariance(InputIt first1, InputIt last1, InputIt first2, UnaryOperation op) {
  return covariance(first1, last1, first2, mean(first1, last1),
                    mean(first2, first2 + std::distance(first1, last1)), op);
}

template <class InputIt1, class InputIt2, class FP, class UnaryOperation, class>
double weighted_covariance(InputIt1 first1, InputIt1 last1, InputIt1 first2, InputIt2 weight_first,
                           FP mean1, FP mean2, UnaryOperation op) {
  double sum_of_squared_devs = 0;
  double weight_sum = 0;
  for (; first1 != last1; ++first1, ++first2, ++weight_first) {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_USELESS_CAST
    const auto n1 = static_cast<double>(op(*first1)) - static_cast<double>(mean1);
    const auto n2 = static_cast<double>(op(*first2)) - static_cast<double>(mean2);
    sum_of_squared_devs += static_cast<double>(*weight_first) * n1 * n2;
    weight_sum += static_cast<double>(*weight_first);
    DISABLE_WARNING_POP
  }

  return sum_of_squared_devs / weight_sum;
}

template <class InputIt1, class InputIt2, class UnaryOperation, class>
double weighted_covariance(InputIt1 first1, InputIt1 last1, InputIt1 first2, InputIt2 weight_first,
                           UnaryOperation op) {
  return weighted_covariance(first1, last1, first2, weight_first, mean(first1, last1),
                             mean(first2, first2 + std::distance(first1, last1)), op);
}

}  // namespace modle::stats

// IWYU pragma: private, include "modle/math.hpp"
