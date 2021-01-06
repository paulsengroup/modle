#pragma once

#include <absl/strings/match.h>  // for EndsWith
#include <absl/time/clock.h>     // for Now
#include <absl/time/time.h>      // for FormatDuration, operator-, Time
#include <fmt/format.h>          // for format, print

#include <array>
#include <boost/asio/thread_pool.hpp>
#include <cassert>
#include <cstdint>     // for uint32_t
#include <cstdio>      // for stderr
#include <filesystem>  // for create_directories
#include <mutex>

#include "modle/bigwig.hpp"    // for write_range
#include "modle/contacts.hpp"  // for ContactMatrix<>::Header, ContactMatrix
#include "modle/cooler.hpp"
#include "modle/hdf5.hpp"
#include "modle_tools/config.hpp"   // for config
#include "modle_tools/convert.hpp"  // for convert_to_hic, convert_to_tsv
#include "modle_tools/eval.hpp"     // for Transformation, Cross, Linear

namespace modle::tools {

void convert_subcmd(const modle::tools::config& c) {
  auto argv = modle::utils::init_juicer_tools_argv(c.path_to_juicer_tools, c.juicer_tools_mem);
  std::filesystem::create_directories(c.output_base_name);
  switch (c.output_format) {
    case config::output_format::hic:
      modle::tools::convert_to_hic(c, argv);
      break;
    case config::output_format::cooler:
      // TODO: Handle this case upstream of this branch.
      break;
    case config::output_format::tsv:
      boost::asio::thread_pool tpool(std::min(c.path_to_input_matrices.size(), c.nthreads));
      modle::tools::convert_to_tsv(tpool, c);
      tpool.join();
      break;
  }
}

void eval_subcmd(const modle::tools::config& c) {
  assert(c.compute_spearman || c.compute_pearson);  // NOLINT
  assert(c.path_to_input_matrices.size() == 1);     // NOLINT
  const auto& path_to_input_cmatrix = c.path_to_input_matrices.front();

  const auto bin_size =
      static_cast<std::size_t>(hdf5::read_attribute_int(path_to_input_cmatrix, "bin-size"));

  auto chr_list =  // This cannot be made const
      select_chromosomes_for_eval(path_to_input_cmatrix, c.path_to_reference_matrix, bin_size);
  if (chr_list.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Files '{}' and '{}' have 0 chromosomes in common. Make sure you are not trying "
                   "to compare different genome assemblies (chromosomes needs to have the same "
                   "name and size in order to qualify for comparison)"),
        path_to_input_cmatrix, c.path_to_reference_matrix));
  }
  {
    std::vector<std::string_view> chr_names(chr_list.size());
    std::transform(chr_list.begin(), chr_list.end(), chr_names.begin(), [](const auto& p) {
      return std::string_view{p.first.data(), p.first.size()};
    });
    if (chr_list.size() == 1) {
      fmt::print(stderr, FMT_STRING("Computing correlation for chromosome: '{}'\n"),
                 chr_names.front());
    } else {
      fmt::print(stderr,
                 FMT_STRING("Computing correlation for the following {} chromosomes: '{}'\n"),
                 chr_list.size(), absl::StrJoin(chr_names, "', '"));
    }
  }

  c.compute_pearson ? fmt::format(FMT_STRING("{}_spearman_linear_r.bw"), c.output_base_name) : "";
  const auto out_path_pv_linear_spearman =
      c.compute_pearson ? fmt::format(FMT_STRING("{}_spearman_linear_pv.bw"), c.output_base_name)
                        : "";
  const auto out_path_corr_cross_spearman =
      c.compute_pearson ? fmt::format(FMT_STRING("{}_spearman_cross_r.bw"), c.output_base_name)
                        : "";
  const auto out_path_pv_cross_spearman =
      c.compute_pearson ? fmt::format(FMT_STRING("{}_spearman_cross_pv.bw"), c.output_base_name)
                        : "";
  const auto& bn = c.output_base_name;
  // Init files and write bw header
  auto* bw_corr_linear_pearson =
      c.compute_pearson
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_pearson_linear_r.bw"), chr_list)
          : nullptr;
  auto* bw_pv_linear_pearson =
      c.compute_pearson
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_pearson_linear_pv.bw"), chr_list)
          : nullptr;
  auto* bw_corr_cross_pearson =
      c.compute_pearson
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_pearson_cross_r.bw"), chr_list)
          : nullptr;
  auto* bw_pv_cross_pearson =
      c.compute_pearson
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_pearson_cross_pv.bw"), chr_list)
          : nullptr;
  auto* bw_corr_linear_spearman =
      c.compute_spearman
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_spearman_linear_r.bw"), chr_list)
          : nullptr;
  auto* bw_pv_linear_spearman =
      c.compute_spearman
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_spearman_linear_pv.bw"), chr_list)
          : nullptr;
  auto* bw_corr_cross_spearman =
      c.compute_spearman
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_spearman_cross_r.bw"), chr_list)
          : nullptr;
  auto* bw_pv_cross_spearman =
      c.compute_spearman
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_spearman_cross_pv.bw"), chr_list)
          : nullptr;

  auto ref_cooler = cooler::Cooler(c.path_to_reference_matrix, cooler::Cooler::READ_ONLY, bin_size);
  auto input_cooler = cooler::Cooler(path_to_input_cmatrix, cooler::Cooler::READ_ONLY, bin_size);
  const auto nrows = (c.diagonal_width / bin_size) + (c.diagonal_width % bin_size != 0);  // NOLINT
  assert(nrows != 0);                                                                     // NOLINT

  std::vector<double> pc_linear_corr_buff;
  std::vector<double> pc_linear_pval_buff;

  std::vector<double> pc_cross_corr_buff;
  std::vector<double> pc_cross_pval_buff;

  std::vector<double> sc_linear_corr_buff;
  std::vector<double> sc_linear_pval_buff;

  std::vector<double> sc_cross_corr_buff;
  std::vector<double> sc_cross_pval_buff;

  auto pcc = [&](std::string_view chr_name, absl::Span<const uint32_t> v1,
                 absl::Span<const uint32_t> v2, std::size_t ncols, std::vector<double>& corr_buff,
                 std::vector<double>& pval_buff, Transformation t) {
    if (!c.compute_pearson) {
      return;
    }
    const auto t0 = absl::Now();
    compute_pearson_over_range(v1, v2, corr_buff, pval_buff, nrows, ncols, t);
    switch (t) {
      case Transformation::Linear:
        bigwig::write_range(std::string{chr_name}, corr_buff, 0, bin_size, bin_size,
                            bw_corr_linear_pearson);
        bigwig::write_range(std::string{chr_name}, pval_buff, 0, bin_size, bin_size,
                            bw_pv_linear_pearson);
        fmt::print(stderr,
                   FMT_STRING("Pearson \"linear\" correlation calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
      case Transformation::Cross:
        bigwig::write_range(std::string{chr_name}, corr_buff, 0, bin_size, bin_size,
                            bw_corr_cross_pearson);
        bigwig::write_range(std::string{chr_name}, pval_buff, 0, bin_size, bin_size,
                            bw_pv_cross_pearson);
        fmt::print(stderr,
                   FMT_STRING("Pearson \"cross\" correlation calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
    }
  };

  auto src = [&](std::string_view chr_name, absl::Span<const uint32_t> v1,
                 absl::Span<const uint32_t> v2, std::size_t ncols, std::vector<double>& corr_buff,
                 std::vector<double>& pval_buff, Transformation t) {
    if (!c.compute_spearman) {
      return;
    }
    const auto t0 = absl::Now();
    compute_spearman_over_range(v1, v2, corr_buff, pval_buff, nrows, ncols, t);
    switch (t) {
      case Transformation::Linear:
        bigwig::write_range(std::string{chr_name}, corr_buff, 0, bin_size, bin_size,
                            bw_corr_linear_spearman);
        bigwig::write_range(std::string{chr_name}, pval_buff, 0, bin_size, bin_size,
                            bw_pv_linear_spearman);
        fmt::print(stderr,
                   FMT_STRING("Spearman \"linear\" correlation calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
      case Transformation::Cross:
        bigwig::write_range(std::string{chr_name}, corr_buff, 0, bin_size, bin_size,
                            bw_corr_cross_spearman);
        bigwig::write_range(std::string{chr_name}, pval_buff, 0, bin_size, bin_size,
                            bw_pv_cross_spearman);
        fmt::print(stderr,
                   FMT_STRING("Spearman \"cross\" correlation calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
    }
  };

  std::array<std::thread, 4> threads;

  for (const auto& chr : chr_list) {
    const auto& chr_name = chr.first;
    auto t0 = absl::Now();
    fmt::print(stderr, FMT_STRING("Reading contacts for '{}' into memory...\n"), chr_name);
    const auto cmatrix1 = ref_cooler.cooler_to_cmatrix(chr_name, nrows);
    fmt::print(stderr, FMT_STRING("Read {}x{} reference matrix in {} using {:.2f} MB of RAM.\n"),
               cmatrix1.n_rows(), cmatrix1.n_cols(), absl::FormatDuration(absl::Now() - t0),
               cmatrix1.get_matrix_size_in_mb());
    t0 = absl::Now();
    const auto cmatrix2 = input_cooler.cooler_to_cmatrix(chr_name, nrows);
    fmt::print(stderr, FMT_STRING("Read {}x{} input matrix in {} using {:.2f} MB of RAM.\n"),
               cmatrix2.n_rows(), cmatrix2.n_cols(), absl::FormatDuration(absl::Now() - t0),
               cmatrix2.get_matrix_size_in_mb());

    if (cmatrix1.n_cols() != cmatrix2.n_cols() || cmatrix1.n_rows() != cmatrix2.n_rows()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("An error occurred while computing the correlation for '{}' between files "
                     "'{}' and '{}': Contact matrices should have the same shape "
                     "m1=[{}][{}], m2=[{}][{}]"),
          chr_name, c.path_to_reference_matrix, path_to_input_cmatrix, cmatrix1.n_rows(),
          cmatrix1.n_cols(), cmatrix2.n_rows(), cmatrix2.n_cols()));
    }

    const auto ncols = cmatrix1.n_cols();

    const auto& v1 = cmatrix1.get_raw_count_vector();
    const auto& v2 = cmatrix2.get_raw_count_vector();
    assert(v1.size() == v2.size());  // NOLINT

    fmt::print(stderr, FMT_STRING("Computing correlation(s) for '{}'...\n"), chr_name);

    threads[0] = std::thread(pcc, chr_name, v1, v2, ncols, std::ref(pc_linear_corr_buff),
                             std::ref(pc_linear_pval_buff), Transformation::Linear);
    threads[1] = std::thread(pcc, chr_name, v1, v2, ncols, std::ref(pc_cross_corr_buff),
                             std::ref(pc_cross_pval_buff), Transformation::Cross);
    threads[2] = std::thread(src, chr_name, v1, v2, ncols, std::ref(sc_linear_corr_buff),
                             std::ref(sc_linear_pval_buff), Transformation::Linear);
    threads[3] = std::thread(src, chr_name, v1, v2, ncols, std::ref(sc_cross_corr_buff),
                             std::ref(sc_cross_pval_buff), Transformation::Cross);
    for (auto& t : threads) {
      t.join();
    }
  }

  bigwig::close_bigwig_file(bw_corr_linear_pearson);
  bigwig::close_bigwig_file(bw_pv_linear_pearson);
  bigwig::close_bigwig_file(bw_corr_cross_pearson);
  bigwig::close_bigwig_file(bw_pv_cross_pearson);
  bigwig::close_bigwig_file(bw_corr_linear_spearman);
  bigwig::close_bigwig_file(bw_pv_linear_spearman);
  bigwig::close_bigwig_file(bw_corr_cross_spearman);
  bigwig::close_bigwig_file(bw_pv_cross_spearman);
}

}  // namespace modle::tools
