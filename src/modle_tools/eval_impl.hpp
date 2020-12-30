#pragma once

#include <absl/time/clock.h>
#include <absl/types/span.h>
#include <fmt/printf.h>

#include "modle/correlation.hpp"
#include "modle/suppress_compiler_warnings.hpp"

namespace modle::tools {

template <typename N>
void slice_range_w_cross_method(absl::Span<const N> vin, std::vector<N> &vout, std::size_t nrows,
                                std::size_t offset) {
  static_assert(std::is_arithmetic<N>::value, "N should be a numeric type.");

  std::size_t idx = offset * nrows;
  for (std::size_t j = 0; j < nrows; ++j) {
    if (idx >= vin.size()) {
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      std::fill(vout.begin() + j, vout.begin() + nrows, 0);
      DISABLE_WARNING_POP
      break;
    }
    vout[j] = vin[idx];
    idx += nrows + 1;
  }
  idx = offset * nrows - 1;
  for (std::size_t j = nrows; j < 2 * nrows; ++j) {
    if (idx > (offset * nrows) + offset) {
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      std::fill(vout.begin() + j, vout.end(), 0);
      break;
      DISABLE_WARNING_POP
    }
    vout[j] = vin[idx++];
  }
}

template <typename N>
void slice_range_w_linear_method(absl::Span<const N> vin, std::vector<N> &vout, std::size_t nrows,
                                 std::size_t offset) {
  static_assert(std::is_arithmetic<N>::value, "N should be a numeric type.");

  std::size_t idx = offset * nrows;
  for (std::size_t j = 0; j < vout.size(); ++j) {
    if (idx >= vin.size() || idx < offset * nrows) {
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      std::fill(vout.begin() + j, vout.end(), 0);
      break;
      DISABLE_WARNING_POP
    }
    vout[j] = vin[idx];
    idx += (nrows * (j % 2 == 0)) + 1;
  }
}

template <typename N>
void slice_range(absl::Span<const N> vin, std::vector<N> &vout, std::size_t nrows, Transformation t,
                 std::size_t offset) {
  vout.resize(2 * nrows - 1);
  switch (t) {
    case Transformation::Cross:
      slice_range_w_cross_method(vin, vout, nrows, offset);
      break;
    case Transformation::Linear:
      slice_range_w_linear_method(vin, vout, nrows, offset);
      break;
    default:
      assert(false);  // This code should be unreachable
  }
}

template <typename N1, typename N2>
std::pair<std::vector<double>, std::vector<double>> compute_pearson_over_range(
    absl::Span<const N1> vin1, absl::Span<const N2> vin2, std::size_t nrows, std::size_t ncols,
    Transformation t) {
  assert(vin1.size() == vin2.size());
  std::vector<double> pcc_vals(ncols);
  std::vector<double> pvals(ncols);
  std::vector<N1> sub_vin1(2 * nrows - 1);
  std::vector<N2> sub_vin2(2 * nrows - 1);
  const auto step = ncols / 25;
  auto t0 = absl::Now();

  for (std::size_t i = 0; i < ncols; ++i) {
    slice_range(vin1, sub_vin1, nrows, t, i);
    slice_range(vin2, sub_vin2, nrows, t, i);
    pcc_vals[i] = correlation::compute_pearson(sub_vin1, sub_vin2);
    pvals[i] = correlation::compute_pearson_significance(pcc_vals[i], sub_vin1.size());

    if ((i + 1) % step == 0) {
      fmt::print(stderr, FMT_STRING("Processed {}/{} bins ({:.4f}%) ({:.2f} corr/s)\n"), i, ncols,
                 100.0 * static_cast<double>(i) / static_cast<double>(ncols),
                 static_cast<double>(step) / absl::ToDoubleSeconds(absl::Now() - t0));
      t0 = absl::Now();
    }
  }

  return {std::move(pcc_vals), std::move(pvals)};
}

template <typename N1, typename N2>
std::pair<std::vector<double>, std::vector<double>> compute_pearson_over_range(
    const std::vector<N1> &vin1, const std::vector<N2> &vin2, std::size_t nrows, std::size_t ncols,
    Transformation t) {
  return compute_pearson_over_range(absl::MakeConstSpan(vin1), absl::MakeConstSpan(vin2), nrows,
                                    ncols, t);
}

template <typename N1, typename N2>
std::pair<std::vector<double>, std::vector<double>> compute_spearman_over_range(
    absl::Span<const N1> vin1, absl::Span<const N2> vin2, std::size_t nrows, std::size_t ncols,
    Transformation t) {
  assert(vin1.size() == vin2.size());
  std::vector<double> rho_vals(ncols);
  std::vector<double> pvals(ncols);
  std::vector<N1> sub_vin1(2 * nrows - 1);
  std::vector<N2> sub_vin2(2 * nrows - 1);
  const auto step = ncols / 25;
  auto t0 = absl::Now();

  for (std::size_t i = 0; i < ncols; ++i) {
    slice_range(vin1, sub_vin1, nrows, t, i);
    slice_range(vin2, sub_vin2, nrows, t, i);
    rho_vals[i] = correlation::compute_spearman(sub_vin1, sub_vin2);
    pvals[i] = correlation::compute_spearman_significance(rho_vals[i], sub_vin1.size());

    if ((i + 1) % step == 0) {
      fmt::print(stderr, FMT_STRING("Processed {}/{} bins ({:.4f}%) ({:.2f} corr/s)\n"), i, ncols,
                 100.0 * static_cast<double>(i) / static_cast<double>(ncols),
                 static_cast<double>(step) / absl::ToDoubleSeconds(absl::Now() - t0));
      t0 = absl::Now();
    }
  }

  return std::make_pair(rho_vals, pvals);
}

template <typename N1, typename N2>
std::pair<std::vector<double>, std::vector<double>> compute_spearman_over_range(
    const std::vector<N1> &vin1, const std::vector<N2> &vin2, std::size_t nrows, std::size_t ncols,
    Transformation t) {
  return compute_spearman_over_range(absl::MakeConstSpan(vin1), absl::MakeConstSpan(vin2), nrows,
                                     ncols, t);
}
}  // namespace modle::tools