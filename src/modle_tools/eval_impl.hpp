#pragma once

#include <absl/time/clock.h>
#include <fmt/printf.h>

#include <range/v3/range.hpp>

#include "modle/correlation.hpp"
#include "modle/suppress_compiler_warnings.hpp"

namespace modle::tools {

template <typename Rng, typename N>
void slice_range_w_cross_method(Rng input_rng, std::vector<N> &output_rng, std::size_t nrows,
                                std::size_t offset) {
  static_assert(std::is_arithmetic<N>::value, "N should be a numeric type.");
  static_assert(ranges::input_range<Rng>, "Rng should be an input range.");
  static_assert(std::is_convertible<ranges::range_value_type_t<Rng>, N>::value,
                "The value type of Rng is not convertible to type N.");

  std::size_t idx = offset * nrows;
  for (std::size_t j = 0; j < nrows; ++j) {
    if (idx >= input_rng.size()) {
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      std::fill(output_rng.begin() + j, output_rng.begin() + nrows, 0);
      DISABLE_WARNING_POP
      break;
    }
    output_rng[j] = input_rng[idx];
    idx += nrows + 1;
  }
  idx = offset * nrows - 1;
  for (std::size_t j = nrows; j < 2 * nrows; ++j) {
    if (idx > (offset * nrows) + offset) {
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      std::fill(output_rng.begin() + j, output_rng.end(), 0);
      break;
      DISABLE_WARNING_POP
    }
    output_rng[j] = input_rng[idx++];
  }
}

template <typename Rng, typename N>
void slice_range_w_linear_method(Rng input_rng, std::vector<N> &output_rng, std::size_t nrows,
                                 std::size_t offset) {
  static_assert(std::is_arithmetic<N>::value, "N should be a numeric type.");
  static_assert(ranges::input_range<Rng>, "Rng should be an input range.");
  static_assert(std::is_convertible<ranges::range_value_type_t<Rng>, N>::value,
                "The value type of Rng is not convertible to type N.");

  std::size_t idx = offset * nrows;
  for (std::size_t j = 0; j < 2 * nrows; ++j) {
    if (idx >= input_rng.size() || idx < offset * nrows) {
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      std::fill(output_rng.begin() + j, output_rng.end(), 0);
      break;
      DISABLE_WARNING_POP
    }
    output_rng[j] = input_rng[idx];
    idx += (nrows * (j % 2 == 0)) + 1;
  }
}

template <typename Rng, typename N>
void slice_range(Rng input_rng, std::vector<N> &output_rng, std::size_t nrows, Transformation t,
                 std::size_t offset) {
  output_rng.resize(2 * nrows - 1);
  switch (t) {
    case Transformation::Cross:
      slice_range_w_cross_method(input_rng, output_rng, nrows, offset);
      break;
    case Transformation::Linear:
      slice_range_w_linear_method(input_rng, output_rng, nrows, offset);
      break;
    default:
      assert(false);  // This code should be unreachable
  }
}

template <typename Rng>
std::pair<std::vector<double>, std::vector<double>> compute_pearson_over_range(Rng r1, Rng r2,
                                                                               std::size_t nrows,
                                                                               std::size_t ncols,
                                                                               Transformation t) {
  assert(r1.size() == r2.size());
  std::vector<double> pcc_vals(ncols);
  std::vector<double> pvals(ncols);
  std::vector<ranges::range_value_type_t<Rng>> v1(2 * nrows - 1);
  std::vector<ranges::range_value_type_t<Rng>> v2(2 * nrows - 1);
  for (std::size_t i = 0; i < ncols; ++i) {
    slice_range(r1, v1, nrows, t, i);
    slice_range(r2, v2, nrows, t, i);
    pcc_vals[i] = correlation::compute_pearson(v1, v2);
    pvals[i] = correlation::compute_pearson_significance(pcc_vals[i], v1.size());
  }

  return std::make_pair(pcc_vals, pvals);
}

template <typename Rng>
std::pair<std::vector<double>, std::vector<double>> compute_spearman_over_range(Rng r1, Rng r2,
                                                                                std::size_t nrows,
                                                                                std::size_t ncols,
                                                                                Transformation t) {
  assert(r1.size() == r2.size());
  std::vector<double> rho_vals(ncols);
  std::vector<double> pvals(ncols);
  std::vector<ranges::range_value_type_t<Rng>> v1(2 * nrows - 1);
  std::vector<ranges::range_value_type_t<Rng>> v2(2 * nrows - 1);
  for (std::size_t i = 0; i < ncols; ++i) {
    slice_range(r1, v1, nrows, t, i);
    slice_range(r2, v2, nrows, t, i);
    rho_vals[i] = correlation::compute_spearman(v1, v2);
    pvals[i] = correlation::compute_spearman_significance(rho_vals[i], v1.size());
  }

  return std::make_pair(rho_vals, pvals);
}
}  // namespace modle::tools