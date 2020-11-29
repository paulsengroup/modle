#pragma once

#include <cstdint>
#include <vector>

#include "modle/contacts.hpp"
#include "modle_tools/cli.hpp"

namespace modle::tools {

enum Transformation { Linear, Cross };

template <typename Rng>
[[nodiscard]] std::pair<std::vector<double>, std::vector<double>> compute_pearson_over_range(
    Rng r1, Rng r2, std::size_t nrows, std::size_t ncols, Transformation t);

template <typename Rng>
[[nodiscard]] std::pair<std::vector<double>, std::vector<double>> compute_spearman_over_range(
    Rng r1, Rng r2, std::size_t nrows, std::size_t ncols, Transformation t);
}  // namespace modle::tools

#include "../../eval_impl.hpp"