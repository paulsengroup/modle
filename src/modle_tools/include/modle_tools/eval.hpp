#pragma once

#include <absl/types/span.h>

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "modle/contacts.hpp"
#include "modle_tools/config.hpp"

namespace modle::tools {

enum Transformation : uint8_t { Linear, Cross };

[[nodiscard]] inline std::vector<std::pair<std::string, int64_t>> select_chromosomes_for_eval(
    std::string_view path_to_cooler1, std::string_view path_to_cooler2, std::size_t bin_size);

template <typename N1, typename N2>
inline void compute_pearson_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                                       std::vector<double>& pcc_buff,
                                       std::vector<double>& pval_buff, std::size_t nrows,
                                       Transformation t);
template <typename N1, typename N2>
inline void compute_pearson_over_range(const std::vector<N1>& vin1, const std::vector<N2>& vin2,
                                       std::vector<double>& pcc_buff,
                                       std::vector<double>& pval_buff, std::size_t nrows,
                                       Transformation t);

template <typename N1, typename N2>
[[nodiscard]] inline std::pair<std::vector<double>, std::vector<double>>
compute_spearman_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                            std::vector<double>& rho_buff, std::vector<double>& pval_buff,
                            std::size_t nrows, Transformation t);
template <typename N1, typename N2>
[[nodiscard]] inline std::pair<std::vector<double>, std::vector<double>>
compute_spearman_over_range(const std::vector<N1>& vin1, const std::vector<N2>& vin2,
                            std::vector<double>& rho_buff, std::vector<double>& pval_buff,
                            std::size_t nrows, Transformation t);

}  // namespace modle::tools

#include "../../eval_impl.hpp"
