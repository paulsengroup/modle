// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>  // for declval, enable_if_t

#include "modle/common/common.hpp"  // for usize
#include "modle/common/utils.hpp"   // for identity

namespace modle::stats {

template <class T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr T abs_diff(T a, T b) noexcept;

//  Mean
//
// Compute the arithmetic mean of the items between the first and last iterators
template <class InputIt, typename UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double mean(InputIt first, InputIt last,
                                 UnaryOperation op = utils::identity());

template <class InputIt1, class InputIt2, typename UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt1>())>>>
[[nodiscard]] double weighted_mean(InputIt1 first, InputIt1 last, InputIt2 weight_first,
                                   UnaryOperation op = utils::identity());
/////////

//  Moving average
//
// Compute the moving average of the items between the first and last iterators given \p
// window_size If \p window_size >= std::distance(begin, end), return math::mean(begin, end)
template <class InputIt, class OutputIt, class I, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_integral_v<I> &&
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
inline usize moving_average(InputIt first, InputIt last, OutputIt output_first, I window_size,
                            UnaryOperation op = utils::identity());
/////////

//  Sum of squared deviations
//
template <class InputIt, class FP, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_floating_point_v<FP> &&
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double sum_of_squared_deviations(InputIt first, InputIt last, FP sample_mean,
                                                      UnaryOperation op = utils::identity());
template <class InputIt, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double sum_of_squared_deviations(InputIt first, InputIt last,
                                                      UnaryOperation op = utils::identity());
/////////

//  Variance
//
template <class InputIt, class FP, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_floating_point_v<FP> &&
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double variance(InputIt first, InputIt last, FP sample_mean,
                                     UnaryOperation op = utils::identity());
template <class InputIt, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double variance(InputIt first, InputIt last,
                                     UnaryOperation op = utils::identity());
/////////

//  Standard deviation
//
template <class InputIt, class FP, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_floating_point_v<FP> &&
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double standard_dev(InputIt first, InputIt last, FP sample_mean,
                                         UnaryOperation op = utils::identity());
template <class InputIt, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double standard_dev(InputIt first, InputIt last,
                                         UnaryOperation op = utils::identity());
/////////

//  Covariance
//
template <class InputIt, class FP, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_floating_point_v<FP> &&
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double covariance(InputIt first1, InputIt last1, InputIt first2, FP mean1,
                                       FP mean2, UnaryOperation op = utils::identity());
template <class InputIt, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double covariance(InputIt first1, InputIt last1, InputIt first2,
                                       UnaryOperation op = utils::identity());

template <class InputIt1, class InputIt2, class FP, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_floating_point_v<FP> &&
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt1>())>>>
[[nodiscard]] inline double weighted_covariance(InputIt1 first1, InputIt1 last1, InputIt1 first2,
                                                InputIt2 weight_first, FP mean1, FP mean2,
                                                UnaryOperation op = utils::identity());
template <class InputIt1, class InputIt2, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt1>())>>>
[[nodiscard]] inline double weighted_covariance(InputIt1 first1, InputIt1 last1, InputIt1 first2,
                                                InputIt2 weight_first,
                                                UnaryOperation op = utils::identity());

/////////

// SED
//

template <class InputIt, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double sed(InputIt first1, InputIt last1, InputIt first2,
                                UnaryOperation op = utils::identity()) noexcept;
template <class InputIt1, class InputIt2, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt1>())>>>
[[nodiscard]] inline double weighted_sed(InputIt1 first1, InputIt1 last1, InputIt1 first2,
                                         InputIt2 weight_first,
                                         UnaryOperation op = utils::identity()) noexcept;
/////////

// RMSE
//
template <class InputIt, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt>())>>>
[[nodiscard]] inline double rmse(InputIt first1, InputIt last1, InputIt first2,
                                 UnaryOperation op = utils::identity());
template <class InputIt1, class InputIt2, class UnaryOperation = utils::identity,
          class = std::enable_if_t<
              std::is_invocable_v<UnaryOperation, decltype(*std::declval<InputIt1>())>>>
[[nodiscard]] inline double weighted_rmse(InputIt1 first1, InputIt1 last1, InputIt1 first2,
                                          InputIt2 weight_first,
                                          UnaryOperation op = utils::identity());
/////////

}  // namespace modle::stats
#include "../../../descriptive_impl.hpp"  // IWYU pragma: export
