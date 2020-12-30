#include "modle_tools/tools.hpp"

#include <absl/strings/match.h>  // for EndsWith
#include <absl/time/clock.h>     // for Now
#include <absl/time/time.h>      // for FormatDuration, operator-, Time
#include <fmt/format.h>          // for format, print

#include <boost/asio/thread_pool.hpp>
#include <cassert>
#include <cstdint>     // for uint32_t
#include <cstdio>      // for stderr
#include <filesystem>  // for create_directories
#include <mutex>

#include "modle/bigwig.hpp"    // for write_range
#include "modle/contacts.hpp"  // for ContactMatrix<>::Header, ContactMatrix
#include "modle/cooler.hpp"
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

  auto chr_list =  // This cannot be made const
      select_chromosomes_for_eval(path_to_input_cmatrix, c.path_to_reference_matrix);
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
  std::filesystem::create_directories(c.output_base_name);

  absl::flat_hash_map<std::pair<std::string, int64_t>, std::vector<double>> corr_linear;
  corr_linear.reserve(chr_list.size());
  std::transform(chr_list.begin(), chr_list.end(), std::inserter(corr_linear, corr_linear.end()),
                 [](const auto& p) { return std::make_pair(p, std::vector<double>{}); });
  auto pv_linear = corr_linear;
  auto corr_cross = corr_linear;
  auto pv_cross = corr_linear;

  if (c.compute_pearson) {  // TODO: Deduplicate this code
    boost::asio::thread_pool tpool(std::min(chr_list.size(), c.nthreads));
    // Generate file names
    const auto out_path_corr_linear =
        fmt::format(FMT_STRING("{}_pearson_linear_r.bw"), c.output_base_name);
    const auto out_path_pv_linear =
        fmt::format(FMT_STRING("{}_pearson_linear_pv.bw"), c.output_base_name);
    const auto out_path_corr_cross =
        fmt::format(FMT_STRING("{}_pearson_cross_r.bw"), c.output_base_name);
    const auto out_path_pv_cross =
        fmt::format(FMT_STRING("{}_pearson_cross_pv.bw"), c.output_base_name);

    // Init files and write bw header
    auto* bw_corr_linear = bigwig::init_bigwig_file(out_path_corr_linear, chr_list);
    auto* bw_pv_linear = bigwig::init_bigwig_file(out_path_pv_linear, chr_list);
    auto* bw_corr_cross = bigwig::init_bigwig_file(out_path_corr_cross, chr_list);
    auto* bw_pv_cross = bigwig::init_bigwig_file(out_path_pv_cross, chr_list);

    const auto bin_size =
        static_cast<std::size_t>(cooler::read_attribute_int(path_to_input_cmatrix, "bin-size"));

    for (const auto& chr : chr_list) {
      boost::asio::post(tpool, [f1 = c.path_to_reference_matrix, f2 = path_to_input_cmatrix,
                                chr_name = chr.first, &corr_linear = corr_linear.at(chr),
                                &pv_linear = pv_linear.at(chr), &corr_cross = corr_cross.at(chr),
                                &pv_cross = pv_cross.at(chr), &c, &bin_size]() {
        auto t0 = absl::Now();
        fmt::print(
            stderr,
            FMT_STRING("Reading contacts for '{}' from the reference matrix '{}' into memory..."),
            chr_name, f1);

        auto cmatrix1 = cooler::cooler_to_cmatrix(f1, chr_name, c.diagonal_width, bin_size);
        fmt::print(stderr,
                   FMT_STRING(" DONE! Read a {}x{} matrix in {} using {:.2f} MB of RAM.\nReading "
                              "contacts for '{}' from file '{}' into memory..."),
                   cmatrix1.n_rows(), cmatrix1.n_cols(), absl::FormatDuration(absl::Now() - t0),
                   cmatrix1.get_matrix_size_in_mb(), chr_name, f2);

        t0 = absl::Now();
        auto cmatrix2 = cooler::cooler_to_cmatrix(f2, chr_name, c.diagonal_width, bin_size);
        fmt::print(stderr, FMT_STRING(" DONE! Read a {}x{} matrix in {} using {:.2f} MB of RAM.\n"),
                   cmatrix1.n_rows(), cmatrix1.n_cols(), absl::FormatDuration(absl::Now() - t0),
                   cmatrix1.get_matrix_size_in_mb());

        if (cmatrix1.n_cols() != cmatrix2.n_cols() || cmatrix1.n_rows() != cmatrix2.n_rows()) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("An error occurred while computing the correlation for '{}' between files "
                         "'{}' and '{}': Contact matrices should have the same shape "
                         "m1=[{}][{}], m2=[{}][{}]"),
              chr_name, f1, f2, cmatrix1.n_rows(), cmatrix1.n_cols(), cmatrix2.n_rows(),
              cmatrix2.n_cols()));
        }
        const auto nrows = cmatrix1.n_rows();
        const auto ncols = cmatrix1.n_cols();

        const auto& v1 = cmatrix1.get_raw_count_vector();
        const auto& v2 = cmatrix2.get_raw_count_vector();
        assert(v1.size() == v2.size());  // NOLINT
        t0 = absl::Now();
        {
          auto pearson = compute_pearson_over_range(v1, v2, nrows, ncols, Transformation::Linear);
          corr_linear = std::move(pearson.first);
          pv_linear = std::move(pearson.second);
        }
        fmt::print(stderr, FMT_STRING("Pearson \"Linear\" for '{}' computed in {}!\n"), chr_name,
                   absl::FormatDuration(absl::Now() - t0));

        t0 = absl::Now();
        {
          auto pearson = compute_pearson_over_range(v1, v2, nrows, ncols, Transformation::Cross);
          corr_cross = std::move(pearson.first);
          pv_cross = std::move(pearson.second);
        }
        fmt::print(stderr, FMT_STRING("Pearson \"Cross\" for '{}' computed in {}!\n"), chr_name,
                   absl::FormatDuration(absl::Now() - t0));
      });
    }
    tpool.join();

    const auto t0 = absl::Now();
    fmt::print(stderr, "Writing correlation values and significance to file...");
    bigwig::write_range(corr_linear, 0, bin_size, bin_size, bw_corr_linear);
    bigwig::write_range(pv_linear, 0, bin_size, bin_size, bw_pv_linear);
    bigwig::write_range(corr_cross, 0, bin_size, bin_size, bw_corr_cross);
    bigwig::write_range(pv_cross, 0, bin_size, bin_size, bw_pv_cross);
    fmt::print(stderr, FMT_STRING(" DONE! Writing took {}!\n"),
               absl::FormatDuration(absl::Now() - t0));
  }

  if (c.compute_spearman) {
    boost::asio::thread_pool tpool(std::min(chr_list.size(), c.nthreads));
    // Generate file names
    const auto out_path_corr_linear =
        fmt::format(FMT_STRING("{}_spearman_linear_rho.bw"), c.output_base_name);
    const auto out_path_pv_linear =
        fmt::format(FMT_STRING("{}_spearman_linear_pv.bw"), c.output_base_name);
    const auto out_path_corr_cross =
        fmt::format(FMT_STRING("{}_spearman_cross_rho.bw"), c.output_base_name);
    const auto out_path_pv_cross =
        fmt::format(FMT_STRING("{}_spearman_cross_pv.bw"), c.output_base_name);

    // Init files and write bw header
    auto* bw_corr_linear = bigwig::init_bigwig_file(out_path_corr_linear, chr_list);
    auto* bw_pv_linear = bigwig::init_bigwig_file(out_path_pv_linear, chr_list);
    auto* bw_corr_cross = bigwig::init_bigwig_file(out_path_corr_cross, chr_list);
    auto* bw_pv_cross = bigwig::init_bigwig_file(out_path_pv_cross, chr_list);

    const auto bin_size =
        static_cast<std::size_t>(cooler::read_attribute_int(path_to_input_cmatrix, "bin-size"));

    for (const auto& chr : chr_list) {
      boost::asio::post(tpool, [f1 = c.path_to_reference_matrix, f2 = path_to_input_cmatrix,
                                chr_name = chr.first, &corr_linear = corr_linear.at(chr),
                                &pv_linear = pv_linear.at(chr), &corr_cross = corr_cross.at(chr),
                                &pv_cross = pv_cross.at(chr), &c, &bin_size]() {
        auto t0 = absl::Now();
        fmt::print(
            stderr,
            FMT_STRING("Reading contacts for '{}' from the reference matrix '{}' into memory..."),
            chr_name, f1);

        auto cmatrix1 = cooler::cooler_to_cmatrix(f1, chr_name, c.diagonal_width, bin_size);
        fmt::print(stderr,
                   FMT_STRING(" DONE! Read a {}x{} matrix in {} using {:.2f} MB of RAM.\nReading "
                              "contacts for '{}' from file '{}' into memory..."),
                   cmatrix1.n_rows(), cmatrix1.n_rows(), absl::FormatDuration(absl::Now() - t0),
                   cmatrix1.get_matrix_size_in_mb(), chr_name, f2);

        t0 = absl::Now();
        auto cmatrix2 = cooler::cooler_to_cmatrix(f2, chr_name, c.diagonal_width, bin_size);
        fmt::print(stderr, FMT_STRING(" DONE! Read a {}x{} matrix in {} using {:.2f} MB of RAM.\n"),
                   cmatrix1.n_rows(), cmatrix1.n_rows(), absl::FormatDuration(absl::Now() - t0),
                   cmatrix1.get_matrix_size_in_mb());

        if (cmatrix1.n_cols() != cmatrix2.n_cols() || cmatrix1.n_rows() != cmatrix2.n_rows()) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("An error occurred while computing the correlation for '{}' between files "
                         "'{}' and '{}': Contact matrices should have the same shape "
                         "m1=[{}][{}], m2=[{}][{}]"),
              chr_name, f1, f2, cmatrix1.n_rows(), cmatrix1.n_cols(), cmatrix2.n_rows(),
              cmatrix2.n_cols()));
        }
        const auto nrows = cmatrix1.n_rows();
        const auto ncols = cmatrix1.n_cols();

        const auto& v1 = cmatrix1.get_raw_count_vector();
        const auto& v2 = cmatrix2.get_raw_count_vector();
        assert(v1.size() == v2.size());  // NOLINT
        t0 = absl::Now();
        {
          auto pearson = compute_spearman_over_range(v1, v2, nrows, ncols, Transformation::Linear);
          corr_linear = std::move(pearson.first);
          pv_linear = std::move(pearson.second);
        }
        fmt::print(stderr, FMT_STRING("Spearman \"Linear\" for '{}' computed in {}!\n"), chr_name,
                   absl::FormatDuration(absl::Now() - t0));

        t0 = absl::Now();
        {
          auto pearson = compute_spearman_over_range(v1, v2, nrows, ncols, Transformation::Cross);
          corr_cross = std::move(pearson.first);
          pv_cross = std::move(pearson.second);
        }
        fmt::print(stderr, FMT_STRING("Spearman \"Cross\" for '{}' computed in {}!\n"), chr_name,
                   absl::FormatDuration(absl::Now() - t0));
      });
    }
    tpool.join();

    const auto t0 = absl::Now();
    fmt::print(stderr, "Writing correlation values and significance to file...");
    bigwig::write_range(corr_linear, 0, bin_size, bin_size, bw_corr_linear);
    bigwig::write_range(pv_linear, 0, bin_size, bin_size, bw_pv_linear);
    bigwig::write_range(corr_cross, 0, bin_size, bin_size, bw_corr_cross);
    bigwig::write_range(pv_cross, 0, bin_size, bin_size, bw_pv_cross);
    fmt::print(stderr, FMT_STRING(" DONE! Writing took {}!\n"),
               absl::FormatDuration(absl::Now() - t0));
  }
}

}  // namespace modle::tools