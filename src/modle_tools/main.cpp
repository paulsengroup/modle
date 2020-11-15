#include "modle/contacts.hpp"
#include "modle/correlation.hpp"
#include "modle_tools/convert.hpp"
#include "modle_tools/eval.hpp"
#include "modle_tools/utils.hpp"

void convert(const modle::tools::config& c) {
  auto argv = modle::tools::utils::init_juicer_tools_argv(c);
  std::filesystem::create_directories(c.out_dir);
  if (c.convert_to_hic) modle::tools::convert::convert_to_hic(c, argv);
  if (c.convert_to_tsv) modle::tools::convert::convert_to_tsv(c);
}

void eval(const modle::tools::config& c) {
  modle::ContactMatrix<uint32_t> cmatrix(c.path_to_input_matrix);
  auto reference_cmatrix = absl::EndsWith(c.path_to_reference_cmatrix, ".hic")
                               ? modle::tools::eval::parse_hic_matrix(c)
                               : modle::ContactMatrix<uint32_t>(c.path_to_reference_cmatrix);
  assert(cmatrix.n_cols() == reference_cmatrix.n_cols() &&
         cmatrix.n_rows() == reference_cmatrix.n_rows());
  modle::CorrelationTest<uint32_t> corr(cmatrix.get_raw_count_vector(),
                                        reference_cmatrix.get_raw_count_vector());
  if (c.compute_spearman) {
    const auto max_size = cmatrix.n_cols() * cmatrix.n_rows();
    const auto spearman = corr.compute_spearman(std::min(max_size, c.sliding_window_size),
                                                std::min(max_size - 1, c.sliding_window_overlap));
    for (const auto& [rho, pv] : spearman) {
      fmt::fprintf(stderr, "rho=%.4G;\tpv=%.4G\n", rho, pv);
    }
  }
}

int main(int argc, char** argv) {
  modle::tools::Cli cli(argc, argv);
  auto c = cli.parse_arguments();
  if (!cli.is_ok()) return c.exit_code;

  std::filesystem::create_directories(c.tmp_dir);
  try {
    switch (cli.get_subcommand()) {
      case modle::tools::Cli::subcommand::convert:
        convert(c);
        break;
      case modle::tools::Cli::subcommand::eval:
        eval(c);
        break;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in main.cpp should be unreachable! If you see "
            "this message, please open an issue on GitHub");
    }
  } catch (const std::runtime_error& err) {
    fmt::fprintf(stderr, "FAILURE: %s.\n", err.what());
    if (!c.keep_tmp_files && std::filesystem::is_empty(c.tmp_dir))
      std::filesystem::remove_all(c.tmp_dir);
    return 1;
  }
  if (!c.keep_tmp_files && std::filesystem::is_empty(c.tmp_dir))
    std::filesystem::remove_all(c.tmp_dir);
  return 0;
}
