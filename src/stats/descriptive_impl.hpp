// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>     // for sqrt
#include <iterator>  // for distance
#include <numeric>   // for accumulate
#include <type_traits>
#include <vector>  // for vector

#include "modle/common/common.hpp"  // for isize

namespace modle::stats {
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
usize moving_average(InputIt input_begin, InputIt input_end, OutputIt output_begin,
                     const I window_size, UnaryOperation op) {
  using OutputItValueT = typename OutputIt::value_type;
  if (static_cast<isize>(window_size) >= std::distance(input_begin, input_end)) {
    *output_begin = OutputItValueT(stats::mean(input_begin, input_end, op));
    return 1;
  }

  for (auto it1 = input_begin, it2 = input_begin + static_cast<isize>(window_size);
       it2 != input_end; ++it1, ++it2, ++output_begin) {
    *output_begin = OutputItValueT(mean(it1, it2, op));
  }
  return usize(std::distance(input_begin, input_end)) - usize(window_size);
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
  return stats::sum_of_squared_deviations(input_begin, input_end,
                                          stats::mean(input_begin, input_end, op), op);
}

template <class InputIt, class FP, class UnaryOperation, class>
double variance(InputIt input_begin, InputIt input_end, const FP sample_mean, UnaryOperation op) {
  if (input_begin == input_end) {
    return 0.0;
  }
  return stats::sum_of_squared_deviations(input_begin, input_end, sample_mean, op) /
         static_cast<double>(std::distance(input_begin, input_end));
}

template <class InputIt, class UnaryOperation, class>
double variance(InputIt input_begin, InputIt input_end, UnaryOperation op) {
  return stats::variance(input_begin, input_end, stats::mean(input_begin, input_end, op), op);
}

template <class InputIt, class FP, class UnaryOperation, class>
double standard_dev(InputIt input_begin, InputIt input_end, const FP sample_mean,
                    UnaryOperation op) {
  if (input_begin == input_end) {
    return 0.0;
  }
  return std::sqrt(stats::variance(input_begin, input_end, sample_mean, op));
}

template <class InputIt, class UnaryOperation, class>
double standard_dev(InputIt input_begin, InputIt input_end, UnaryOperation op) {
  return stats::standard_dev(input_begin, input_end, stats::mean(input_begin, input_end, op), op);
}
}  // namespace modle::stats

// IWYU pragma: private, include "modle/math.hpp"
