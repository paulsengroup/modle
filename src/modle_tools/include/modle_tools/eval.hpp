#pragma once

#include <absl/types/span.h>  // for Span

#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for uint64_t, int64_t, uint_fast8_t
#include <string>       // for string
#include <string_view>  // for string_view
#include <utility>      // for pair
#include <vector>       // for vector

namespace modle::tools {

enum Transformation : uint_fast8_t { Linear, Cross };

[[nodiscard]] inline std::vector<std::pair<std::string, int64_t>> select_chromosomes_for_eval(
    std::string_view path_to_cooler1, std::string_view path_to_cooler2, size_t bin_size);

template <typename N1, typename N2>
inline void compute_pearson_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                                       std::vector<double>& pcc_buff,
                                       std::vector<double>& pval_buff, size_t nrows,
                                       Transformation t);
template <typename N1, typename N2>
inline void compute_pearson_over_range(const std::vector<N1>& vin1, const std::vector<N2>& vin2,
                                       std::vector<double>& pcc_buff,
                                       std::vector<double>& pval_buff, size_t nrows,
                                       Transformation t);

template <typename N1, typename N2>
inline void compute_spearman_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                                        std::vector<double>& rho_buff,
                                        std::vector<double>& pval_buff, size_t nrows,
                                        Transformation t);
template <typename N1, typename N2>
inline void compute_spearman_over_range(const std::vector<N1>& vin1, const std::vector<N2>& vin2,
                                        std::vector<double>& rho_buff,
                                        std::vector<double>& pval_buff, size_t nrows,
                                        Transformation t);

template <typename N1, typename N2>
inline void compute_sed_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                                   std::vector<uint64_t>& buff, size_t nrows, size_t ncols,
                                   Transformation t);

template <typename N1, typename N2>
inline void compute_sed_over_range(const std::vector<N1>& vin1, const std::vector<N2>& vin2,
                                   std::vector<uint64_t>& buff, size_t nrows, size_t ncols,
                                   Transformation t);

template <typename N1, typename N2>
void compute_euc_dist_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                                 std::vector<double>& buff, size_t nrows, size_t ncols,
                                 Transformation t);

template <typename N1, typename N2>
void compute_euc_dist_over_range(const std::vector<N1>& vin1, const std::vector<N2>& vin2,
                                 std::vector<double>& buff, size_t nrows, size_t ncols,
                                 Transformation t);
}  // namespace modle::tools

#include "../../eval_impl.hpp"  // IWYU pragma: keep
