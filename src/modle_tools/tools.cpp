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

#include "modle/bigwig.hpp"           // for write_range
#include "modle/contacts.hpp"         // for ContactMatrix<>::Header, ContactMatrix
#include "modle/juicer_contacts.hpp"  // for run_juicer_dump_and_parse_contacts
#include "modle/suppress_compiler_warnings.hpp"
#include "modle_tools/config.hpp"   // for config
#include "modle_tools/convert.hpp"  // for convert_to_hic, convert_to_tsv
#include "modle_tools/eval.hpp"     // for Transformation, Cross, Linear

namespace modle::tools {

void convert_subcmd(const modle::tools::config& c) {
  assert(c.convert_to_hic || c.convert_to_tsv);
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_CONVERSION
  boost::asio::thread_pool tpool(  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
      std::min((c.convert_to_hic + c.convert_to_tsv) * c.path_to_input_matrices.size(),
               c.nthreads));
  DISABLE_WARNING_POP
  auto argv = modle::utils::init_juicer_tools_argv(c.path_to_juicer_tools, c.juicer_tools_mem);
  std::filesystem::create_directories(c.out_dir);
  if (c.convert_to_hic) {
    modle::tools::convert_to_hic(tpool, c, argv);
  }
  if (c.convert_to_tsv) {
    modle::tools::convert_to_tsv(tpool, c);
  }
  tpool.join();
}

void eval_subcmd(const modle::tools::config& c) {
  assert(c.compute_spearman || c.compute_pearson);
  boost::asio::thread_pool tpool(
      std::min(c.path_to_input_matrices.size(), c.nthreads));

  auto argv = modle::utils::init_juicer_tools_argv(c.path_to_juicer_tools, c.juicer_tools_mem);
  std::filesystem::create_directories(c.out_dir);

  for (auto i = 0U; i < c.path_to_input_matrices.size(); ++i) {
    boost::asio::post(tpool, [&, i]() {
      modle::ContactMatrix<uint32_t> cmatrix(c.path_to_input_matrices[i]);
      auto reference_cmatrix =
          absl::EndsWith(c.path_to_reference_matricx, ".hic")
              ? modle::juicer_contacts::run_juicer_dump_and_parse_contacts(
                    c.path_to_input_matrices[i], c.path_to_reference_matricx, c.chr_name_hic,
                    c.chr_offset_hic, c.path_to_juicer_tools, c.juicer_tools_mem)
              : modle::ContactMatrix<uint32_t>(c.path_to_reference_matricx);
      assert(cmatrix.n_cols() == reference_cmatrix.n_cols() &&
             cmatrix.n_rows() == reference_cmatrix.n_rows());
      const auto& v1 = cmatrix.get_raw_count_vector();
      const auto& v2 = reference_cmatrix.get_raw_count_vector();
      const auto header = ContactMatrix<uint32_t>::parse_header(c.path_to_input_matrices[i]);
      assert(v1.size() == v2.size());
      if (c.compute_pearson) {
        auto t0 = absl::Now();
        auto pearson = compute_pearson_over_range(v1, v2, cmatrix.n_rows(), cmatrix.n_cols(),
                                                  Transformation::Linear);
        fmt::print(stderr, "Pearson \"Linear\" computed in {}!\n",
                   absl::FormatDuration(absl::Now() - t0));
        auto out_path_pearson =
            fmt::format("{}/{}_pearson_linear_r.bw", c.out_dir, header.chr_name);
        bigwig::write_range(header.chr_name, header.end, pearson.first, header.start,
                            header.bin_size, header.bin_size, out_path_pearson);
        out_path_pearson = fmt::format("{}/{}_pearson_linear_pv.bw", c.out_dir, header.chr_name);
        bigwig::write_range(header.chr_name, header.end, pearson.second, header.start,
                            header.bin_size, header.bin_size, out_path_pearson);

        t0 = absl::Now();
        pearson = compute_pearson_over_range(v1, v2, cmatrix.n_rows(), cmatrix.n_cols(),
                                             Transformation::Cross);
        fmt::print(stderr, "Pearson \"Cross\" computed in {}!\n",
                   absl::FormatDuration(absl::Now() - t0));
        out_path_pearson = fmt::format("{}/{}_pearson_cross_r.bw", c.out_dir, header.chr_name);
        bigwig::write_range(header.chr_name, header.end, pearson.first, header.start,
                            header.bin_size, header.bin_size, out_path_pearson);
        out_path_pearson = fmt::format("{}/{}_pearson_cross_pv.bw", c.out_dir, header.chr_name);
        bigwig::write_range(header.chr_name, header.end, pearson.second, header.start,
                            header.bin_size, header.bin_size, out_path_pearson);
      }
      if (c.compute_spearman) {
        auto t0 = absl::Now();
        auto pearson = compute_spearman_over_range(v1, v2, cmatrix.n_rows(), cmatrix.n_cols(),
                                                   Transformation::Linear);
        fmt::print(stderr, "Spearman \"Linear\" computed in {}!\n",
                   absl::FormatDuration(absl::Now() - t0));
        auto out_path_pearson =
            fmt::format("{}/{}_spearman_linear_r.bw", c.out_dir, header.chr_name);
        bigwig::write_range(header.chr_name, header.end, pearson.first, header.start,
                            header.bin_size, header.bin_size, out_path_pearson);
        out_path_pearson = fmt::format("{}/{}_spearman_linear_pv.bw", c.out_dir, header.chr_name);
        bigwig::write_range(header.chr_name, header.end, pearson.second, header.start,
                            header.bin_size, header.bin_size, out_path_pearson);

        t0 = absl::Now();
        pearson = compute_spearman_over_range(v1, v2, cmatrix.n_rows(), cmatrix.n_cols(),
                                              Transformation::Cross);
        fmt::print(stderr, "Spearman \"Cross\" computed in {}!\n",
                   absl::FormatDuration(absl::Now() - t0));
        out_path_pearson = fmt::format("{}/{}_spearman_cross_r.bw", c.out_dir, header.chr_name);
        bigwig::write_range(header.chr_name, header.end, pearson.first, header.start,
                            header.bin_size, header.bin_size, out_path_pearson);
        out_path_pearson = fmt::format("{}/{}_spearman_cross_pv.bw", c.out_dir, header.chr_name);
        bigwig::write_range(header.chr_name, header.end, pearson.second, header.start,
                            header.bin_size, header.bin_size, out_path_pearson);
      }
    });
  }
  tpool.join();
}
}  // namespace modle::tools