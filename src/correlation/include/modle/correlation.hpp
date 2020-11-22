#pragma once
#include <cstdint>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

#include "range/v3/range.hpp"
#include "range/v3/view.hpp"

namespace modle::correlation {

[[nodiscard]] double compute_spearman_significance(double rho, std::size_t n);
[[nodiscard]] double compute_pearson_significance(double pcc, std::size_t n);
//[[nodiscard]] double compute_kendall_b_significance(uint32_t nc, uint32_t nd, uint32_t size,
//                                                    uint32_t tie1, uint32_t tie2);

template <typename Rng>
[[nodiscard]] double compute_pearson(Rng r1, Rng r2);
template <typename Rng>
[[nodiscard]] double compute_spearman(Rng r1, Rng r2);
/*
template <typename Rng>
[[nodiscard]] double compute_kendall_b(Rng r1, Rng r2);
*/

template <typename N, typename Rng>
[[nodiscard]] std::pair<std::vector<double>, std::vector<double>> compute_pearson(
    const std::vector<N>& v1, const std::vector<N>& v2, std::size_t window_span,
    std::size_t window_overlap = 0);

template <typename N, typename Rng>
[[nodiscard]] std::pair<std::vector<double>, std::vector<double>> compute_spearman(
    const std::vector<N>& v1, const std::vector<N>& v2, std::size_t window_span,
    std::size_t window_overlap = 0);

/*
template <typename N, typename Rng>
[[nodiscard]] std::pair<std::vector<double>, std::vector<double>> compute_kendall_b(
    const std::vector<N>& v1, const std::vector<N>& v2, std::size_t window_span,
    std::size_t window_overlap = 0);
*/
}  // namespace modle::correlation

#include "../../correlation_impl.hpp"
