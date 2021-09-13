// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/algorithm/container.h>      // for c_set_intersection
#include <absl/container/btree_set.h>      // for btree_set
#include <absl/container/flat_hash_map.h>  // for flat_hash_map, BitMask, operator!=
#include <absl/meta/type_traits.h>         // for remove_reference_t
#include <absl/strings/match.h>            // for StartsWith
#include <absl/strings/str_cat.h>          // for StrCat
#include <absl/strings/string_view.h>      // for string_view
#include <absl/strings/strip.h>            // for StripPrefix
#include <absl/time/clock.h>               // for Now
#include <absl/time/time.h>                // for FormatDuration, operator-, Time
#include <absl/types/span.h>               // for Span, MakeConstSpan
#include <fmt/format.h>                    // for make_format_args, vformat_to, FMT...
#include <fmt/ostream.h>                   // for formatbuf<>::int_type
#include <spdlog/spdlog.h>                 // for info, warn

#include <algorithm>                        // for fill, max, transform
#include <array>                            // for array, array<>::value_type
#include <boost/filesystem/operations.hpp>  // for create_directories
#include <boost/filesystem/path.hpp>        // for path, operator<<
#include <cassert>                          // for assert
#include <cmath>                            // for sqrt
#include <cstdint>                          // for int64_t, uint32_t, uint_fast8_t
#include <cstdio>                           // for size_t
#include <exception>                        // for exception
#include <functional>                       // for ref
#include <initializer_list>                 // for initializer_list
#include <iosfwd>                           // for streamsize
#include <iterator>                         // for back_insert_iterator, insert_iter...
#include <memory>                           // for unique_ptr
#include <stdexcept>                        // for runtime_error, logic_error
#include <string>                           // for string, basic_string
#include <thread>                           // for thread
#include <type_traits>                      // for is_arithmetic
#include <utility>                          // for pair, move, make_pair
#include <vector>                           // for vector

#include "modle/bed.hpp"                                // for Parser
#include "modle/bigwig.hpp"                             // for write_range, init_bigwig_file
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...
#include "modle/common/utils.hpp"                       // for chrom_less_than_operator
#include "modle/contacts.hpp"                           // for ContactMatrix
#include "modle/cooler.hpp"                             // for Cooler, Cooler::READ_ONLY
#include "modle/correlation.hpp"                        // for compute_pearson_significance, com...
#include "modle/hdf5.hpp"                               // for read_attribute_int
#include "modle_tools/config.hpp"                       // for eval_config
#include "modle_tools/tools.hpp"                        // for eval_subcmd

namespace modle::tools {

enum Transformation : uint_fast8_t { Linear, Cross };

std::vector<std::pair<std::string, int64_t>> select_chromosomes_for_eval(
    std::string_view path_to_cooler1, std::string_view path_to_cooler2, size_t bin_size) {
  std::vector<std::string> str_buff;
  std::vector<int64_t> int_buff;

  auto build_chrom_set = [&](std::string_view path_to_cooler) {
    try {
      cooler::Cooler c(std::string{path_to_cooler}, cooler::Cooler::READ_ONLY, bin_size);
      c.get_chrom_names(str_buff);
      c.get_chrom_sizes(int_buff);

      assert(str_buff.size() == c.get_nchroms());
      assert(int_buff.size() == str_buff.size());

      absl::btree_set<std::pair<std::string, int64_t>> chrom_set;

      for (auto i = 0UL; i < str_buff.size(); ++i) {
        chrom_set.emplace(str_buff[i], int_buff[i]);
      }
      return chrom_set;

    } catch (const std::runtime_error &e) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("An error occurred while reading file '{}': {}"), path_to_cooler, e.what()));
    }
  };

  std::vector<std::pair<std::string, int64_t>> chrom_intersection;
  const auto chrom_set1 = build_chrom_set(path_to_cooler1);
  const auto chrom_set2 = build_chrom_set(path_to_cooler2);

  absl::c_set_intersection(chrom_set1, chrom_set2, std::back_inserter(chrom_intersection),
                           [](const auto &c1, const auto &c2) {
                             return modle::utils::chrom_less_than_operator(c1, c2);
                           });

  return chrom_intersection;
}

template <typename N>
void slice_range_w_cross_method(absl::Span<const N> vin, std::vector<N> &vout, size_t nrows,
                                size_t offset) {
  static_assert(std::is_arithmetic<N>::value, "N should be a numeric type.");

  size_t idx = offset * nrows;
  for (size_t j = 0; j < nrows; ++j) {
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
  for (size_t j = nrows; j < vout.size(); ++j) {
    if (idx > (offset * nrows) + offset) {
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      std::fill(vout.begin() + j, vout.end(), 0);
      DISABLE_WARNING_POP
      break;
    }
    vout[j] = vin[idx++];
  }
}

template <typename N>
void slice_range_w_linear_method(absl::Span<const N> vin, std::vector<N> &vout, size_t nrows,
                                 size_t offset) {
  static_assert(std::is_arithmetic<N>::value, "N should be a numeric type.");

  size_t idx = offset * nrows;
  for (size_t j = 0; j < vout.size(); ++j) {
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
void slice_range(absl::Span<const N> vin, std::vector<N> &vout, size_t nrows, Transformation t,
                 size_t offset) {
  vout.resize(2 * nrows - 1);
  switch (t) {
    case Transformation::Cross:
      slice_range_w_cross_method(vin, vout, nrows, offset);
      break;
    case Transformation::Linear:
      slice_range_w_linear_method(vin, vout, nrows, offset);
      break;
    default:
      throw std::logic_error("Unreachable code");
  }
}

template <typename N1, typename N2>
void compute_pearson_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                                std::vector<double> &pcc_buff, std::vector<double> &pval_buff,
                                size_t nrows, size_t ncols, Transformation t) {
  assert(vin1.size() == vin2.size());
  pcc_buff.resize(ncols);
  pval_buff.resize(ncols);
  std::vector<N1> sub_vin1(2 * nrows - 1);
  std::vector<N2> sub_vin2(2 * nrows - 1);
  // const auto step = ncols / 25;
  // auto t0 = absl::Now();

  for (size_t i = 0; i < ncols; ++i) {
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
                                size_t nrows, size_t ncols, Transformation t) {
  compute_pearson_over_range(absl::MakeConstSpan(vin1), absl::MakeConstSpan(vin2), pcc_buff,
                             pval_buff, nrows, ncols, t);
}

template <typename N1, typename N2>
void compute_spearman_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                                 std::vector<double> &rho_buff, std::vector<double> &pval_buff,
                                 size_t nrows, size_t ncols, Transformation t) {
  assert(vin1.size() == vin2.size());
  rho_buff.resize(ncols);
  pval_buff.resize(ncols);
  std::vector<N1> sub_vin1(2 * nrows - 1);
  std::vector<N2> sub_vin2(2 * nrows - 1);

  for (size_t i = 0; i < ncols; ++i) {
    slice_range(vin1, sub_vin1, nrows, t, i);
    slice_range(vin2, sub_vin2, nrows, t, i);
    rho_buff[i] = correlation::compute_spearman(sub_vin1, sub_vin2);
    pval_buff[i] = correlation::compute_spearman_significance(rho_buff[i], sub_vin1.size());
  }
}

template <typename N1, typename N2>
void compute_spearman_over_range(const std::vector<N1> &vin1, const std::vector<N2> &vin2,
                                 std::vector<double> &rho_buff, std::vector<double> &pval_buff,
                                 size_t nrows, size_t ncols, Transformation t) {
  compute_spearman_over_range(absl::MakeConstSpan(vin1), absl::MakeConstSpan(vin2), rho_buff,
                              pval_buff, nrows, ncols, t);
}

template <typename N1, typename N2>
void compute_sed_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                            std::vector<double> &buff, size_t nrows, size_t ncols,
                            Transformation t) {
  assert(vin1.size() == vin2.size());
  buff.resize(ncols);
  std::vector<N1> sub_vin1(2 * nrows - 1);
  std::vector<N2> sub_vin2(2 * nrows - 1);

  for (size_t i = 0; i < ncols; ++i) {
    slice_range(vin1, sub_vin1, nrows, t, i);
    slice_range(vin2, sub_vin2, nrows, t, i);
    buff[i] = static_cast<double>(correlation::compute_sed(sub_vin1, sub_vin2));
  }
}

template <typename N1, typename N2>
void compute_sed_over_range(const std::vector<N1> &vin1, const std::vector<N2> &vin2,
                            std::vector<double> &buff, size_t nrows, size_t ncols,
                            Transformation t) {
  compute_sed_over_range(absl::MakeConstSpan(vin1), absl::MakeConstSpan(vin2), buff, nrows, ncols,
                         t);
}

template <typename N1, typename N2>
void compute_euc_dist_over_range(absl::Span<const N1> vin1, absl::Span<const N2> vin2,
                                 std::vector<double> &buff, size_t nrows, size_t ncols,
                                 Transformation t) {
  compute_sed_over_range(vin1, vin2, buff, nrows, ncols, t);
  std::transform(buff.begin(), buff.end(), buff.begin(), [](const auto n) { return std::sqrt(n); });
}

template <typename N1, typename N2>
void compute_euc_dist_over_range(const std::vector<N1> &vin1, const std::vector<N2> &vin2,
                                 std::vector<double> &buff, size_t nrows, size_t ncols,
                                 Transformation t) {
  compute_euc_dis_over_range(absl::MakeConstSpan(vin1), absl::MakeConstSpan(vin2), buff, nrows,
                             ncols, t);
}

void eval_subcmd(const modle::tools::eval_config &c) {
  assert(c.compute_spearman || c.compute_pearson);  // NOLINT
  const auto bin_size =
      static_cast<size_t>(hdf5::read_attribute_int(c.path_to_input_matrix.string(), "bin-size"));

  auto chrom_list =  // This cannot be made const
      select_chromosomes_for_eval(c.path_to_input_matrix.string(),
                                  c.path_to_reference_matrix.string(), bin_size);
  if (chrom_list.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Files {} and {} have 0 chromosomes in common. Make sure you are not trying "
                   "to compare different genome assemblies (chromosomes needs to have the same "
                   "name and size in order to qualify for comparison)"),
        c.path_to_input_matrix, c.path_to_reference_matrix));
  }

  absl::flat_hash_map<std::string, std::pair<size_t, size_t>> chrom_subranges;
  if (!c.path_to_chrom_subranges.empty()) {
    const auto records = bed::Parser(c.path_to_chrom_subranges).parse_all();
    chrom_subranges.reserve(records.size());
    std::transform(records.begin(), records.end(),
                   std::inserter(chrom_subranges, chrom_subranges.end()),
                   [](const auto &r) -> std::pair<std::string, std::pair<size_t, size_t>> {
                     return {r.chrom, {r.chrom_start, r.chrom_end}};
                   });
  }

  {
    std::vector<std::string_view> chrom_names(chrom_list.size());
    std::transform(chrom_list.begin(), chrom_list.end(), chrom_names.begin(), [](const auto &p) {
      return std::string_view{p.first.data(), p.first.size()};
    });
    if (chrom_list.size() == 1) {
      spdlog::info(FMT_STRING("Computing correlation for chromosome: '{}'"), chrom_names.front());
    } else {
      spdlog::info(FMT_STRING("Computing correlation for the following {} chromosomes: '{}'"),
                   chrom_list.size(), fmt::join(chrom_names, "', '"));
    }
  }

  const auto bn = c.output_base_name.string();

  auto create_bwig_file = [&](std::string_view fname_suffix, bool skip) -> bigwig::file {
    if (skip) {
      return {nullptr, bigwig::close_bigwig_file};
    }
    return bigwig::init_bigwig_file(absl::StrCat(bn, "_", absl::StripPrefix(fname_suffix, "_")),
                                    chrom_list);
  };

  boost::filesystem::create_directories(c.output_base_name.parent_path());
  // Init files and write bw header
  auto bw_corr_linear_pearson = create_bwig_file("pearson_r_linear.bw", !c.compute_pearson);
  auto bw_pv_linear_pearson = create_bwig_file("pearson_pv_linear.bw", !c.compute_pearson);
  auto bw_corr_cross_pearson = create_bwig_file("pearson_r_cross.bw", !c.compute_pearson);
  auto bw_pv_cross_pearson = create_bwig_file("pearson_pv_cross.bw", !c.compute_pearson);
  auto bw_corr_linear_spearman = create_bwig_file("spearman_r_linear.bw", !c.compute_spearman);
  auto bw_pv_linear_spearman = create_bwig_file("spearman_pv_linear.bw", !c.compute_spearman);
  auto bw_corr_cross_spearman = create_bwig_file("spearman_r_cross.bw", !c.compute_spearman);
  auto bw_pv_cross_spearman = create_bwig_file("spearman_pv_cross.bw", !c.compute_spearman);
  auto bw_linear_sed = create_bwig_file("eucl_dist_linear.bw", !c.compute_edist);
  auto bw_cross_sed = create_bwig_file("eucl_dist_cross.bw", !c.compute_edist);

  auto ref_cooler = cooler::Cooler(c.path_to_reference_matrix, cooler::Cooler::READ_ONLY, bin_size);
  auto input_cooler = cooler::Cooler(c.path_to_input_matrix, cooler::Cooler::READ_ONLY, bin_size);
  const auto nrows = (c.diagonal_width / bin_size) + (c.diagonal_width % bin_size != 0);  // NOLINT
  assert(nrows != 0);                                                                     // NOLINT

  std::vector<double> pc_linear_corr_buff;
  std::vector<double> pc_linear_pval_buff;

  std::vector<double> pc_cross_corr_buff;
  std::vector<double> pc_cross_pval_buff;

  std::vector<double> sc_linear_corr_buff;
  std::vector<double> sc_linear_pval_buff;

  std::vector<double> ed_linear_buff;
  std::vector<double> ed_cross_buff;

  std::vector<double> sc_cross_corr_buff;
  std::vector<double> sc_cross_pval_buff;

  auto pcc = [&](std::string_view chrom_name, absl::Span<const uint32_t> v1,
                 absl::Span<const uint32_t> v2, size_t ncols, size_t offset,
                 std::vector<double> &corr_buff, std::vector<double> &pval_buff, Transformation t) {
    if (!c.compute_pearson) {
      return;
    }
    const auto t0 = absl::Now();
    compute_pearson_over_range(v1, v2, corr_buff, pval_buff, nrows, ncols, t);
    switch (t) {
      case Transformation::Linear:
        bigwig::write_range(std::string{chrom_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_linear_pearson);
        bigwig::write_range(std::string{chrom_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_linear_pearson);
        spdlog::info(FMT_STRING("Pearson \"linear\" correlation calculation completed in {}."),
                     absl::FormatDuration(absl::Now() - t0));
        break;
      case Transformation::Cross:
        bigwig::write_range(std::string{chrom_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_cross_pearson);
        bigwig::write_range(std::string{chrom_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_cross_pearson);
        spdlog::info(FMT_STRING("Pearson \"cross\" correlation calculation completed in {}."),
                     absl::FormatDuration(absl::Now() - t0));
        break;
    }
  };

  auto src = [&](std::string_view chrom_name, absl::Span<const uint32_t> v1,
                 absl::Span<const uint32_t> v2, size_t ncols, size_t offset,
                 std::vector<double> &corr_buff, std::vector<double> &pval_buff, Transformation t) {
    if (!c.compute_spearman) {
      return;
    }
    const auto t0 = absl::Now();
    compute_spearman_over_range(v1, v2, corr_buff, pval_buff, nrows, ncols, t);
    switch (t) {
      case Transformation::Linear:
        bigwig::write_range(std::string{chrom_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_linear_spearman);
        bigwig::write_range(std::string{chrom_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_linear_spearman);
        spdlog::info(FMT_STRING("Spearman \"linear\" correlation calculation completed in {}."),
                     absl::FormatDuration(absl::Now() - t0));
        break;
      case Transformation::Cross:
        bigwig::write_range(std::string{chrom_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_cross_spearman);
        bigwig::write_range(std::string{chrom_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_cross_spearman);
        spdlog::info(FMT_STRING("Spearman \"cross\" correlation calculation completed in {}."),
                     absl::FormatDuration(absl::Now() - t0));
        break;
    }
  };

  auto edist = [&](std::string_view chrom_name, absl::Span<const uint32_t> v1,
                   absl::Span<const uint32_t> v2, size_t ncols, size_t offset,
                   std::vector<double> &sed_buff, Transformation t) {
    if (!c.compute_edist) {
      return;
    }
    const auto t0 = absl::Now();
    compute_euc_dist_over_range(v1, v2, sed_buff, nrows, ncols, t);
    switch (t) {
      case Transformation::Linear:
        bigwig::write_range(std::string{chrom_name}, sed_buff, offset, bin_size, bin_size,
                            bw_linear_sed);
        spdlog::info(FMT_STRING("Euclidean dist. \"linear\" calculation completed in {}."),
                     absl::FormatDuration(absl::Now() - t0));
        break;
      case Transformation::Cross:
        bigwig::write_range(std::string{chrom_name}, sed_buff, offset, bin_size, bin_size,
                            bw_cross_sed);
        spdlog::info(FMT_STRING("Euclidean dist. \"cross\" calculation completed in {}."),
                     absl::FormatDuration(absl::Now() - t0));
        break;
    }
  };

  absl::flat_hash_map<std::string, size_t> ref_chrom_idxes;
  absl::flat_hash_map<std::string, size_t> inp_chrom_idxes;

  size_t i = 0;
  for (auto &name : ref_cooler.get_chrom_names()) {
    ref_chrom_idxes.emplace(std::move(name), i++);
  }
  i = 0;
  for (auto &name : input_cooler.get_chrom_names()) {
    inp_chrom_idxes.emplace(std::move(name), i++);
  }

  for (const auto &chrom : chrom_list) {
    std::string q;
    for (const auto &prefix : {"chrom", "chrom", "chrom"}) {
      if (absl::StartsWith(chrom.first, prefix)) {
        q = absl::StripPrefix(chrom.first, prefix);
      } else {
        q = absl::StrCat(prefix, chrom.first);
      }
      if (auto it = ref_chrom_idxes.find(q); it != ref_chrom_idxes.end()) {
        const auto idx = it->second;
        ref_chrom_idxes.emplace(chrom.first, idx);
        ref_chrom_idxes.erase(q);
      }
      if (auto it = inp_chrom_idxes.find(q); it != inp_chrom_idxes.end()) {
        const auto idx = it->second;
        inp_chrom_idxes.emplace(chrom.first, idx);
        inp_chrom_idxes.erase(q);
      }
    }
  }

  std::array<std::thread, 6> threads;  // NOLINT
  for (const auto &chrom : chrom_list) {
    const auto &chrom_name = chrom.first;
    auto chrom_subrange = std::make_pair(0UL, static_cast<size_t>(chrom.second));
    if (!chrom_subranges.empty()) {
      auto it = chrom_subranges.find(chrom_name);

      if (it != chrom_subranges.end()) {
        chrom_subrange = it->second;
      } else {  // Try common prefixes
        for (const auto &prefix : {"chrom", "chrom", "chrom"}) {
          if (absl::StartsWith(chrom_name, prefix)) {
            it = chrom_subranges.find(absl::StripPrefix(chrom_name, prefix));
          } else {
            it = chrom_subranges.find(absl::StrCat(prefix, chrom_name));
          }
          if (it != chrom_subranges.end()) {
            chrom_subrange = it->second;
            break;
          }
        }
      }
    }

    auto t0 = absl::Now();
    spdlog::info(FMT_STRING("Reading contacts for '{}' into memory..."), chrom_name);
    if (!ref_cooler.has_contacts_for_chrom(ref_chrom_idxes.at(chrom_name))) {
      spdlog::warn(FMT_STRING("WARNING: reference contact matrix doesn't have any contacts for "
                              "'{}'. SKIPPING!"),
                   chrom_name);
      continue;
    }

    if (!input_cooler.has_contacts_for_chrom(inp_chrom_idxes.at(chrom_name))) {
      spdlog::warn(FMT_STRING("WARNING: contact matrix doesn't have any contacts "
                              "for '{}'. SKIPPING!"),
                   chrom_name);
      continue;
    }

    auto cmatrix1 = ref_cooler.cooler_to_cmatrix(chrom_name, nrows, chrom_subrange);
    if (c.deplete_contacts_from_reference) {
      cmatrix1.deplete_contacts(c.depletion_multiplier);
    }
    spdlog::info(FMT_STRING("Read {:.2f}M contacts for a {}x{} reference matrix in {} using "
                            "{:.2f} MB of RAM."),
                 static_cast<double>(cmatrix1.get_tot_contacts()) / 1.0e6, cmatrix1.nrows(),
                 cmatrix1.ncols(), absl::FormatDuration(absl::Now() - t0),
                 cmatrix1.get_matrix_size_in_mb());
    t0 = absl::Now();

    const auto cmatrix2 = input_cooler.cooler_to_cmatrix(chrom_name, nrows, chrom_subrange);
    spdlog::info(FMT_STRING("Read {:.2f}M contacts for a {}x{} input matrix in {} using "
                            "{:.2f} MB of RAM."),
                 static_cast<double>(cmatrix2.get_tot_contacts()) / 1.0e6, cmatrix2.nrows(),
                 cmatrix2.ncols(), absl::FormatDuration(absl::Now() - t0),
                 cmatrix2.get_matrix_size_in_mb());

    if (cmatrix1.ncols() != cmatrix2.ncols() || cmatrix1.nrows() != cmatrix2.nrows()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("An error occurred while computing the correlation for "
                                 "'{}' between files "
                                 "{} and {}: Contact matrices should have the same shape "
                                 "m1=[{}][{}], m2=[{}][{}]"),
                      chrom_name, c.path_to_reference_matrix, c.path_to_input_matrix,
                      cmatrix1.nrows(), cmatrix1.ncols(), cmatrix2.nrows(), cmatrix2.ncols()));
    }

    const auto ncols = cmatrix1.ncols();

    const auto &v1 = cmatrix1.get_raw_count_vector();
    const auto &v2 = cmatrix2.get_raw_count_vector();
    assert(v1.size() == v2.size());  // NOLINT

    spdlog::info(FMT_STRING("Computing correlation(s) for '{}'..."), chrom_name);

    threads[0] = std::thread(pcc, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(pc_linear_corr_buff), std::ref(pc_linear_pval_buff),
                             Transformation::Linear);
    threads[1] = std::thread(pcc, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(pc_cross_corr_buff), std::ref(pc_cross_pval_buff),
                             Transformation::Cross);
    threads[2] = std::thread(src, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(sc_linear_corr_buff), std::ref(sc_linear_pval_buff),
                             Transformation::Linear);
    threads[3] = std::thread(src, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(sc_cross_corr_buff), std::ref(sc_cross_pval_buff),
                             Transformation::Cross);

    threads[4] = std::thread(edist, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(ed_linear_buff), Transformation::Linear);

    threads[5] = std::thread(edist, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(ed_cross_buff), Transformation::Cross);
    for (auto &t : threads) {
      t.join();
    }
  }
}

}  // namespace modle::tools
