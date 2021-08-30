#pragma once

#include <absl/types/span.h>

#include <type_traits>

#include "modle/common/utils.hpp"

namespace modle::math {

template <class T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr T abs_diff(T a, T b) noexcept;

//  Mean
//
// Compute the arithmetic mean of the items between the begin and end iterators
template <
    class InputIt, class UnaryOperation = utils::identity,
    class = std::enable_if_t<std::is_invocable_v<UnaryOperation, typename InputIt::value_type>>>
[[nodiscard]] inline double mean(InputIt begin, InputIt end, UnaryOperation op = utils::identity());
/////////

//  Moving average
//
// Compute the moving average of the items between the begin and end iterators given \p window_size
// If \p window_size >= std::distance(begin, end), return math::mean(begin, end)
template <
    class InputIt, class OutputIt, class I, class UnaryOperation = utils::identity,
    class = std::enable_if_t<std::is_integral_v<I> &&
                             std::is_invocable_v<UnaryOperation, typename InputIt::value_type>>>
inline size_t moving_average(InputIt input_begin, InputIt input_end, OutputIt output_begin,
                             I window_size, UnaryOperation op = utils::identity());
/////////

//  Sum of squared deviations
//
template <
    class InputIt, class FP, class UnaryOperation = utils::identity,
    class = std::enable_if_t<std::is_floating_point_v<FP> &&
                             std::is_invocable_v<UnaryOperation, typename InputIt::value_type>>>
[[nodiscard]] inline double sum_of_squared_deviations(InputIt input_begin, InputIt input_end,
                                                      FP sample_mean,
                                                      UnaryOperation op = utils::identity());
template <
    class InputIt, class UnaryOperation = utils::identity,
    class = std::enable_if_t<std::is_invocable_v<UnaryOperation, typename InputIt::value_type>>>
[[nodiscard]] inline double sum_of_squared_deviations(InputIt input_begin, InputIt input_end,
                                                      UnaryOperation op = utils::identity());
/////////

//  Variance
//
template <
    class InputIt, class FP, class UnaryOperation = utils::identity,
    class = std::enable_if_t<std::is_floating_point_v<FP> &&
                             std::is_invocable_v<UnaryOperation, typename InputIt::value_type>>>
[[nodiscard]] inline double variance(InputIt input_begin, InputIt input_end, FP sample_mean,
                                     UnaryOperation op = utils::identity());
template <
    class InputIt, class UnaryOperation = utils::identity,
    class = std::enable_if_t<std::is_invocable_v<UnaryOperation, typename InputIt::value_type>>>
[[nodiscard]] inline double variance(InputIt input_begin, InputIt input_end,
                                     UnaryOperation op = utils::identity());
/////////

//  Standard deviation
//
template <
    class InputIt, class FP, class UnaryOperation = utils::identity,
    class = std::enable_if_t<std::is_floating_point_v<FP> &&
                             std::is_invocable_v<UnaryOperation, typename InputIt::value_type>>>
[[nodiscard]] inline double standard_dev(InputIt input_begin, InputIt input_end, FP sample_mean,
                                         UnaryOperation op = utils::identity());
template <
    class InputIt, class UnaryOperation = utils::identity,
    class = std::enable_if_t<std::is_invocable_v<UnaryOperation, typename InputIt::value_type>>>
[[nodiscard]] inline double standard_dev(InputIt input_begin, InputIt input_end,
                                         UnaryOperation op = utils::identity());
/////////

}  // namespace modle::math
#include "../../../math_impl.hpp"  // IWYU pragma: export
