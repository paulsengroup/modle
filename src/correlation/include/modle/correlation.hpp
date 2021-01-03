#pragma once
#include <absl/types/span.h>

#include <cstdint>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace modle::correlation {

[[nodiscard]] inline double compute_spearman_significance(double rho, std::size_t n);
[[nodiscard]] inline double compute_pearson_significance(double pcc, std::size_t n);
//[[nodiscard]] inline double compute_kendall_b_significance(uint32_t nc, uint32_t nd, uint32_t
// size,
//                                                    uint32_t tie1, uint32_t tie2);

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
[[nodiscard]] inline std::pair<std::vector<double>, std::vector<double>> compute_pearson(
    absl::Span<const N1> v1, absl::Span<const N2> v2, std::size_t window_span,
    std::size_t window_overlap = 0);

template <typename N1, typename N2>
[[nodiscard]] inline std::pair<std::vector<double>, std::vector<double>> compute_spearman(
    absl::Span<const N1> v1, absl::Span<const N2> v2, std::size_t window_span,
    std::size_t window_overlap = 0);

/*
template <typename N, typename Rng>
[[nodiscard]] inline std::pair<std::vector<double>, std::vector<double>> compute_kendall_b(
    const std::vector<N>& v1, const std::vector<N>& v2, std::size_t window_span,
    std::size_t window_overlap = 0);
*/
}  // namespace modle::correlation

#include "../../correlation_impl.hpp"
