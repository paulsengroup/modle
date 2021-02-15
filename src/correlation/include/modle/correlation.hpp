#pragma once
#include <absl/types/span.h>

#include <cstdint>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace modle::correlation {

enum Algorithm : uint8_t { pearson, spearman };

[[nodiscard]] inline double compute_pearson_significance(double pcc, std::size_t n);
[[nodiscard]] inline double compute_spearman_significance(double rho, std::size_t n);

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
    absl::Span<const N1> v1, absl::Span<const N2> v2, Algorithm type, std::size_t window_span,
    std::size_t window_overlap);

template <typename N1, typename N2>
[[nodiscard]] inline double compute_sed(absl::Span<const N1> v1, absl::Span<const N2> v2);
}  // namespace modle::correlation

#include "../../correlation_impl.hpp"
