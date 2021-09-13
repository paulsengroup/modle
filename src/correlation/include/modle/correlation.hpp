// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span

#include <cstddef>  // for size_t
#include <cstdint>  // for uint_fast8_t
#include <utility>  // for pair
#include <vector>   // for vector

namespace modle::correlation {

enum Algorithm : uint_fast8_t { pearson, spearman };

[[nodiscard]] inline double compute_pearson_significance(double pcc, size_t n);
[[nodiscard]] inline double compute_spearman_significance(double rho, size_t n);

template <typename N1, typename N2>
[[nodiscard]] inline double compute_pearson(absl::Span<const N1> v1, absl::Span<const N2> v2);
template <typename N1, typename N2>
[[nodiscard]] inline double compute_pearson(const std::vector<N1>& v1, const std::vector<N2>& v2);
template <typename N1, typename N2>
[[nodiscard]] inline double compute_spearman(absl::Span<const N1> v1, absl::Span<const N2> v2);
template <typename N1, typename N2>
[[nodiscard]] inline double compute_spearman(const std::vector<N1>& v1, const std::vector<N2>& v2);
/*
template <typename Rng>
[[nodiscard]] inline double compute_kendall_b(Rng r1, Rng r2);
*/

template <typename N1, typename N2>
[[nodiscard]] inline std::pair<std::vector<double>, std::vector<double>> compute_corr(
    absl::Span<const N1> v1, absl::Span<const N2> v2, Algorithm type, size_t window_span,
    size_t window_overlap);

template <typename N1, typename N2>
[[nodiscard]] inline double compute_sed(absl::Span<const N1> v1, absl::Span<const N2> v2);

namespace utils {
// Return the indices corresponding to the sorted vector
template <typename N>
[[nodiscard]] inline std::vector<size_t> sort_range_by_idx(absl::Span<const N> v);

template <typename N>
[[nodiscard]] inline std::vector<size_t> sort_range_by_idx(const std::vector<N>& v);

// This returns a vector of integers corresponding to the rank of the values of
// vector v
template <typename N>
[[nodiscard]] inline std::vector<double> compute_element_ranks(absl::Span<const N> v);

template <typename N>
[[nodiscard]] inline std::vector<double> compute_element_ranks(const std::vector<N>& v);
}  // namespace utils

}  // namespace modle::correlation

#include "../../correlation_impl.hpp"        // IWYU pragma: export
#include "../../correlation_utils_impl.hpp"  // IWYU pragma: export

// IWYU pragma: private, include "../../correlation_impl.hpp"
// IWYU pragma: private, include "../../correlation_utils_impl.hpp"
