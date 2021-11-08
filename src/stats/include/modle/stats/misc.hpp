// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>  // for enable_if, is_floating_point
#include <vector>       // for vector

#include "modle/common/common.hpp"  // for usize

namespace modle::stats {

// https://patrickfuller.github.io/gaussian-blur-image-processing-for-scientists-and-engineers-part-4/
template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
[[nodiscard]] inline std::vector<FP> compute_gauss_kernel(FP sigma = 1, FP cutoff = 0.005);

template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
inline void compute_gauss_kernel(std::vector<FP>& buff, FP sigma = 1, FP cutoff = 0.005);

template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
[[nodiscard]] inline std::vector<FP> compute_gauss_kernel(usize size, FP sigma = 1);

template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
inline void compute_gauss_kernel(usize size, std::vector<FP>& buff, FP sigma = 1);
}  // namespace modle::stats

#include "../../../misc_impl.hpp"
