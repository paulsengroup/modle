#include <algorithm>
#include <filesystem>
#include <string>

#include "CLI/CLI.hpp"
#include "absl/strings/str_format.h"
#include "modle/contacts.hpp"
#include "modle/correlation_test.hpp"

namespace modle::helpers::cmatrix_comp {

struct config {
  std::string m1;
  std::string m2;
  std::string out_dir;
  bool force{false};
  uint64_t sliding_window_size{3000};
  uint64_t sliding_window_overlap{2400};
  std::string base_name;
  int exit_code{0};
};

[[nodiscard]] config parse_cli(int argc, char** argv) {
  CLI::App cli;
  cli.name(argv[0]);
  config c;
  // clang-format off
  cli.description("Modle's helper tool. This tool aids the comparison between modle contact matrices and the contact file produced by juicer dump.");
  cli.add_option("-1, --m1", c.m1, "Path to modle's contact matrix")->check(CLI::ExistingFile)->required();
  cli.add_option("-2, --m2", c.m2, "Path to the contact file produced by juicer dump")->check(CLI::ExistingFile)->required();
  cli.add_option("-o, --output-dir", c.out_dir, "Path where to output the wig and CSV files.")->required();
  cli.add_option("-w, --sliding-window-size", c.sliding_window_size, "Sliding window size.")->check(CLI::PositiveNumber)->capture_default_str();
  cli.add_option("--sliding-window-overlap", c.sliding_window_overlap, "Overlap between consecutive sliding-windows.")->check(CLI::PositiveNumber)->capture_default_str();
  cli.add_flag("-f, --force", c.force, "Overwrite existing file(s).")->capture_default_str();
  // clang-format on
  try {
    cli.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    c.exit_code =
        cli.exit(e);  //  This takes care of formatting and printing error messages (if any)
    return c;
  }
  if (!std::filesystem::is_directory(c.out_dir) && std::filesystem::exists(c.out_dir)) {
    throw std::runtime_error(
        absl::StrFormat("'%s' should point to a directory or a non-existing path.", c.out_dir));
  }

  if (c.m1.rfind(".tsv.bz2") == std::string::npos) {
    throw std::runtime_error("The matrix passed to --m1 should be in .tsv.bz2 format.");
  }

  c.base_name =
      absl::StrFormat("%s/%s", c.out_dir, std::string_view(c.m1.data(), c.m1.rfind(".tsv.bz2")));
  if (!c.force) {
    std::vector<std::string> collisions;
    for (std::string_view suffix : {"rho", "tau", "rho_pv", "tau_pv"}) {
      for (std::string_view ext : {"tsv.bz2", "wig"}) {
        if (auto file = absl::StrFormat("%s_%s.%s", c.base_name, suffix, ext);
            std::filesystem::exists(file)) {
          collisions.emplace_back(std::move(file));
        }
      }
    }
    if (!collisions.empty()) {
      throw std::runtime_error(
          absl::StrFormat("Detected %lu file name collisions: refusing to proceed. Pass --force to "
                          "overwrite any existing file.\nColliding file(s):\n - %s",
                          collisions.size(), absl::StrJoin(collisions, "\n - ")));
    }
  }

  if (c.sliding_window_size <= c.sliding_window_overlap) {
    throw std::runtime_error(
        absl::StrFormat("--sliding-window-size should be > --sliding-window-overlap: "
                        "--sliding-window-size=%lu --sliding-window-overlap=%lu.",
                        c.sliding_window_size, c.sliding_window_overlap));
  }

  return c;
}

}  // namespace modle::helpers::cmatrix_comp

int main(int argc, char** argv) {
  auto c = modle::helpers::cmatrix_comp::parse_cli(argc, argv);
  if (c.exit_code != 0) return c.exit_code;
  auto m1 = modle::ContactMatrix<uint32_t>(c.m1);
  auto m2 = modle::ContactMatrix<uint32_t>(c.m2, m1.n_rows(), m1.n_cols());

  modle::CorrelationTest<uint32_t> corr(m1.get_raw_count_vector(), m2.get_raw_count_vector());
  const auto max_size = m1.n_cols() * m1.n_rows();
  const auto spearman = corr.compute_spearman(std::min(max_size, c.sliding_window_size),
                                              std::min(max_size - 1, c.sliding_window_overlap));
  const auto kendall = corr.compute_kendall(std::min(max_size, c.sliding_window_size),
                                            std::min(max_size - 1, c.sliding_window_overlap));

  assert(spearman.size() == kendall.size());
  absl::PrintF("idx\trho\trho_pv\ttau\ttau_pv\n");
  for (auto i = 0UL; i < spearman.size(); ++i) {
    absl::PrintF("%lu\t%.16G\t%.16G\t%.16G\t%.16G\n", i, spearman[i].first, spearman[i].second,
                 kendall[i].first, kendall[i].second);
  }

  return 0;
}
