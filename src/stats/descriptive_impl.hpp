// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>
#include <iterator>
#include <numeric>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"
#include "modle/common/utils.hpp"

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
    total += utils::conditional_static_cast<double>(op(*first)) *
             utils::conditional_static_cast<double>(*weight_first);
    weight_total += utils::conditional_static_cast<double>(*weight_first);
  }

  return total / weight_total;
}

template <class InputIt, class OutputIt, class I, class UnaryOperation, class>
usize moving_average(InputIt first, InputIt last, OutputIt output_first, const I window_size,
                     UnaryOperation op) {
  using OutputItValueT = typename OutputIt::value_type;
  if (static_cast<isize>(window_size) >= std::distance(first, last)) {
    *output_first = OutputItValueT(stats::mean(first, last, op));
    return 1;
  }

  for (auto it1 = first, it2 = first + static_cast<isize>(window_size); it2 != last;
       ++it1, ++it2, ++output_first) {
    *output_first = OutputItValueT(mean(it1, it2, op));
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
    const auto n1 = utils::conditional_static_cast<double>(op(*first1)) -
                    utils::conditional_static_cast<double>(mean1);
    const auto n2 = utils::conditional_static_cast<double>(op(*first2)) -
                    utils::conditional_static_cast<double>(mean2);
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
    const auto n1 = utils::conditional_static_cast<double>(op(*first1)) -
                    utils::conditional_static_cast<double>(mean1);
    const auto n2 = utils::conditional_static_cast<double>(op(*first2)) -
                    utils::conditional_static_cast<double>(mean2);
    sum_of_squared_devs += utils::conditional_static_cast<double>(*weight_first) * n1 * n2;
    weight_sum += utils::conditional_static_cast<double>(*weight_first);
  }

  return sum_of_squared_devs / weight_sum;
}

template <class InputIt1, class InputIt2, class UnaryOperation, class>
double weighted_covariance(InputIt1 first1, InputIt1 last1, InputIt1 first2, InputIt2 weight_first,
                           UnaryOperation op) {
  return weighted_covariance(first1, last1, first2, weight_first, mean(first1, last1),
                             mean(first2, first2 + std::distance(first1, last1)), op);
}

template <class InputIt, class UnaryOperation, class>
double sed(InputIt first1, InputIt last1, InputIt first2, UnaryOperation op) noexcept {
  return weighted_sed(first1, last1, first2, utils::RepeatIterator<double>(1), op);
}
template <class InputIt1, class InputIt2, class UnaryOperation, class>
double weighted_sed(InputIt1 first1, InputIt1 last1, InputIt1 first2, InputIt2 weight_first,
                    UnaryOperation op) noexcept {
  const isize size = std::distance(first1, last1);
  assert(size >= 0);

  double sed = 0;
  for (isize i = 0; i < size; ++i) {
    const auto n1 = op(*(first1 + i));
    const auto n2 = op(*(first2 + i));
    const auto w = *(weight_first + i);
    const auto n = [&]() {
      if constexpr (std::is_unsigned_v<decltype(op(*first1))>) {
        return n1 > n2 ? n1 - n2 : n2 - n1;
      } else {
        return n2 - n1;
      }
    }();
    sed +=
        utils::conditional_static_cast<double>(w) * utils::conditional_static_cast<double>(n * n);
  }

  return std::sqrt(sed);
}

template <class InputIt, class UnaryOperation, class>
double rmse(InputIt first1, InputIt last1, InputIt first2, UnaryOperation op) {
  return weighted_rmse(first1, last1, first2, utils::RepeatIterator<double>(1), op);
}

template <class InputIt1, class InputIt2, class UnaryOperation, class>
double weighted_rmse(InputIt1 first1, InputIt1 last1, InputIt1 first2, InputIt2 weight_first,
                     UnaryOperation op) {
  isize n = 0;
  double tot = 0;
  double tot_w = 0;
  auto it1 = first1;
  auto it2 = first2;
  auto it3 = weight_first;
  for (; it1 != last1; ++n, ++it1, ++it2, ++it3) {
    tot += std::pow(utils::conditional_static_cast<double>(op(*it1) - op(*it2)), 2.0) *
           utils::conditional_static_cast<double>(*it3);
    tot_w += utils::conditional_static_cast<double>(*it3);
  }

  return std::sqrt(tot / tot_w);
}

}  // namespace modle::stats

// IWYU pragma: private, include "modle/math.hpp"
