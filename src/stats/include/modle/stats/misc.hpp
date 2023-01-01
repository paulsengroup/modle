// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>  // for enable_if, is_floating_point
#include <vector>       // for vector

#include "modle/common/common.hpp"  // for usize

namespace modle::stats {

// https://patrickfuller.github.io/gaussian-blur-image-processing-for-scientists-and-engineers-part-4/
template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
[[nodiscard]] inline std::vector<FP> compute_gauss_kernel2d(FP sigma, FP truncate = 4.0);

template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
inline void compute_gauss_kernel2d(std::vector<FP>& buff, FP sigma, FP truncate = 4.0);

template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
[[nodiscard]] inline std::vector<FP> compute_gauss_kernel2d(usize size, FP sigma);

template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
inline void compute_gauss_kernel2d(usize radius, std::vector<FP>& buff, FP sigma);
}  // namespace modle::stats

#include "../../../misc_impl.hpp"
