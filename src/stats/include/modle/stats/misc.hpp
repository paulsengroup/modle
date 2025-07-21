// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>
#include <vector>

#include "modle/common/common.hpp"

namespace modle::stats {

template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
[[nodiscard]] inline std::vector<FP> compute_gauss_kernel2d(FP sigma, FP truncate = 4.0);

template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
inline void compute_gauss_kernel2d(std::vector<FP>& buff, FP sigma, FP truncate = 4.0);

template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
[[nodiscard]] inline std::vector<FP> compute_gauss_kernel2d(usize size, FP sigma);

template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
inline void compute_gauss_kernel2d(usize radius, std::vector<FP>& buff, FP sigma);

template <class InputIt1, class InputIt2, class N = double,
          class = std::enable_if_t<std::is_arithmetic_v<N>>>
[[nodiscard]] constexpr N cross_correlation(InputIt1 kernel_first, InputIt1 kernel_last,
                                            InputIt2 buff_first);

template <class Rng1, class Rng2, class N = double,
          class = std::enable_if_t<std::is_arithmetic_v<N>>>
[[nodiscard]] constexpr N cross_correlation(const Rng1& kernel, const Rng2& buff);

}  // namespace modle::stats

#include "../../../misc_impl.hpp"
