#pragma once

#include <absl/algorithm/container.h>
#include <absl/container/btree_set.h>
#include <absl/strings/match.h>
#include <absl/time/clock.h>
#include <absl/types/span.h>
#include <fmt/format.h>

#include "modle/cooler.hpp"
#include "modle/correlation.hpp"
#include "modle/suppress_compiler_warnings.hpp"
#include "modle/utils.hpp"

namespace modle::tools {
std::vector<std::pair<std::string, int64_t>> select_chromosomes_for_eval(
    std::string_view path_to_cooler1, std::string_view path_to_cooler2, std::size_t bin_size) {
  std::vector<std::string> str_buff;
  std::vector<int64_t> int_buff;

  auto build_chr_set = [&](std::string_view path_to_cooler) {
    try {
      cooler::Cooler c(path_to_cooler, cooler::Cooler::READ_ONLY, bin_size);
      c.get_chr_names(str_buff);
      c.get_chr_sizes(int_buff);

      assert(str_buff.size() == c.get_nchroms());
      assert(int_buff.size() == str_buff.size());

      absl::btree_set<std::pair<std::string, int64_t>> chr_set;

      for (auto i = 0UL; i < str_buff.size(); ++i) {
        chr_set.emplace(str_buff[i], int_buff[i]);
      }
      return chr_set;

    } catch (const std::runtime_error &e) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("An error occurred while reading file '{}': {}"), path_to_cooler, e.what()));
    }
  };

  std::vector<std::pair<std::string, int64_t>> chr_intersection;
  const auto chr_set1 = build_chr_set(path_to_cooler1);
  const auto chr_set2 = build_chr_set(path_to_cooler2);

  absl::c_set_intersection(
      chr_set1, chr_set2, std::back_inserter(chr_intersection),
      [](const auto &c1, const auto &c2) { return modle::utils::chr_less_than_operator(c1, c2); });

  return chr_intersection;
}

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
void compute_pearson_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                                std::vector<double> &pcc_buff, std::vector<double> &pval_buff,
                                std::size_t nrows, std::size_t ncols, Transformation t) {
  assert(vin1.size() == vin2.size());
  pcc_buff.resize(ncols);
  pval_buff.resize(ncols);
  std::vector<N1> sub_vin1(2 * nrows - 1);
  std::vector<N2> sub_vin2(2 * nrows - 1);
  // const auto step = ncols / 25;
  // auto t0 = absl::Now();

  for (std::size_t i = 0; i < ncols; ++i) {
    slice_range(vin1, sub_vin1, nrows, t, i);
    slice_range(vin2, sub_vin2, nrows, t, i);
    pcc_buff[i] = correlation::compute_pearson(sub_vin1, sub_vin2);
    pval_buff[i] = correlation::compute_pearson_significance(pcc_buff[i], sub_vin1.size());

    /*
    if ((i + 1) % step == 0) {
      fmt::print(stderr, FMT_STRING("Processed {}/{} bins ({:.4f}%) ({:.2f} corr/s)\n"), i, ncols,
                 100.0 * static_cast<double>(i) / static_cast<double>(ncols),
                 static_cast<double>(step) / absl::ToDoubleSeconds(absl::Now() - t0));
      t0 = absl::Now();
    }
     */
  }
}

template <typename N1, typename N2>
void compute_pearson_over_range(const std::vector<N1> &vin1, const std::vector<N2> &vin2,
                                std::vector<double> &pcc_buff, std::vector<double> &pval_buff,
                                std::size_t nrows, std::size_t ncols, Transformation t) {
  return compute_pearson_over_range(absl::MakeConstSpan(vin1), absl::MakeConstSpan(vin2), pcc_buff,
                                    pval_buff, nrows, ncols, t);
}

template <typename N1, typename N2>
void compute_spearman_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                                 std::vector<double> &rho_buff, std::vector<double> &pval_buff,
                                 std::size_t nrows, std::size_t ncols, Transformation t) {
  assert(vin1.size() == vin2.size());
  rho_buff.resize(ncols);
  pval_buff.resize(ncols);
  std::vector<N1> sub_vin1(2 * nrows - 1);
  std::vector<N2> sub_vin2(2 * nrows - 1);
  // const auto step = ncols / 25;
  // auto t0 = absl::Now();

  for (std::size_t i = 0; i < ncols; ++i) {
    slice_range(vin1, sub_vin1, nrows, t, i);
    slice_range(vin2, sub_vin2, nrows, t, i);
    rho_buff[i] = correlation::compute_spearman(sub_vin1, sub_vin2);
    pval_buff[i] = correlation::compute_spearman_significance(rho_buff[i], sub_vin1.size());
    /*
        if ((i + 1) % step == 0) {
          fmt::print(stderr, FMT_STRING("Processed {}/{} bins ({:.4f}%) ({:.2f} corr/s)\n"), i,
       ncols, 100.0 * static_cast<double>(i) / static_cast<double>(ncols), static_cast<double>(step)
       / absl::ToDoubleSeconds(absl::Now() - t0)); t0 = absl::Now();
        }
        */
  }
}

template <typename N1, typename N2>
void compute_spearman_over_range(const std::vector<N1> &vin1, const std::vector<N2> &vin2,
                                 std::vector<double> &rho_buff, std::vector<double> &pval_buff,
                                 std::size_t nrows, std::size_t ncols, Transformation t) {
  return compute_spearman_over_range(absl::MakeConstSpan(vin1), absl::MakeConstSpan(vin2), rho_buff,
                                     pval_buff, nrows, ncols, t);
}

}  // namespace modle::tools
