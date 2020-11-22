#include "modle/bigwig.hpp"
#include "modle/contacts.hpp"
#include "modle/correlation.hpp"
#include "modle_tools/convert.hpp"
#include "modle_tools/eval.hpp"
#include "modle_tools/utils.hpp"
#include "range/v3/view.hpp"

namespace modle::tools {

void convert(const modle::tools::config& c) {
  auto argv = modle::tools::utils::init_juicer_tools_argv(c);
  std::filesystem::create_directories(c.out_dir);
  if (c.convert_to_hic) {
    modle::tools::convert_to_hic(c, argv);
  }
  if (c.convert_to_tsv) {
    modle::tools::convert_to_tsv(c);
  }
}

void eval(const modle::tools::config& c) {
  modle::ContactMatrix<uint32_t> cmatrix(c.path_to_input_matrix);
  auto reference_cmatrix = absl::EndsWith(c.path_to_reference_cmatrix, ".hic")
                               ? modle::tools::parse_hic_matrix(c)
                               : modle::ContactMatrix<uint32_t>(c.path_to_reference_cmatrix);
  assert(cmatrix.n_cols() == reference_cmatrix.n_cols() &&
         cmatrix.n_rows() == reference_cmatrix.n_rows());
  const auto& v1 = cmatrix.get_raw_count_vector();
  const auto& v2 = reference_cmatrix.get_raw_count_vector();
  assert(c.compute_spearman || c.compute_pearson);
  const auto header = ContactMatrix<uint32_t>::parse_header(c.path_to_input_matrix);
  const auto window_size = c.sliding_window_size == 0 ? cmatrix.n_rows() : c.sliding_window_size;
  assert(v1.size() == v2.size());
  if (c.compute_pearson) {
    const auto spearman =
        correlation::compute_pearson(v1, v2, window_size, c.sliding_window_overlap);
    auto out_path_spearman = fmt::format("{}/{}_pearson_r.bw", c.out_dir, header.chr_name);
    bigwig::write_range(header.chr_name, header.end, spearman.first, header.start, header.bin_size,
                        header.bin_size, out_path_spearman);
    out_path_spearman = fmt::format("{}/{}_pearson_pv.bw", c.out_dir, header.chr_name);
    bigwig::write_range(header.chr_name, header.end, spearman.second, header.start, header.bin_size,
                        header.bin_size, out_path_spearman);
  }
  if (c.compute_spearman) {
    const auto spearman =
        correlation::compute_spearman(v1, v2, window_size, c.sliding_window_overlap);
    auto out_path_spearman = fmt::format("{}/{}_spearman_rho.bw", c.out_dir, header.chr_name);
    bigwig::write_range(header.chr_name, header.end, spearman.first, header.start, header.bin_size,
                        header.bin_size, out_path_spearman);
    out_path_spearman = fmt::format("{}/{}_spearman_pv.bw", c.out_dir, header.chr_name);
    bigwig::write_range(header.chr_name, header.end, spearman.second, header.start, header.bin_size,
                        header.bin_size, out_path_spearman);
  }
}
}  // namespace modle::tools