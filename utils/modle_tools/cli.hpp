#pragma once
#include <cstdint>
#include <string>

#include "CLI/CLI.hpp"
#include "absl/container/btree_set.h"
#include "absl/strings/match.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"

namespace modle::utils {

struct config {
  std::string path_to_input_matrix;
  std::string out_dir;
  std::string chr_sizes;
  bool convert_to_hic{true};
  bool convert_to_tsv{true};
  bool compress{true};
  bool add_noise{false};
  bool compute_spearman{true};
  bool compute_kendall{true};
  bool force{false};
  bool keep_tmp_files{false};
  uint64_t sliding_window_size{3000};
  uint64_t sliding_window_overlap{2400};
  std::string base_name;
  std::string path_to_juicer_tools;
  std::string tmp_dir{std::filesystem::temp_directory_path()};
  int exit_code{0};
  uint64_t seed{0};
  double noise_stdev{3};
  uint64_t noise_range{100'000};
  uint64_t juicer_tools_mem{2 * 1024 * 1024 * 1024ULL};  // 2GiB
};

class Cli {
 public:
  enum subcommand { convert, eval, help };

 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  config _config{};
  CLI::App _cli{};
  subcommand _subcommand{subcommand::help};
  const absl::btree_set<std::string_view> _allowed_genome_ids{
      // See https://github.com/aidenlab/juicer/wiki/Pre#usage
      {"hg18"},      {"hg19"},      {"hg38"},    {"dMel"},    {"mm9"},     {"mm10"},
      {"anasPlat1"}, {"bTaurus3"},  {"canFam3"}, {"equCab2"}, {"galGal4"}, {"Pf3D7"},
      {"sacCer3"},   {"sCerS288c"}, {"susScr3"}, {"TAIR10"}};

  void MakeCli() {
    // clang-format off
    this->_cli.name(this->_exec_name);
    this->_cli.description("Modle's helper tool. This tool allows to post-process the contact matrix produced by Modle in various ways.");
    this->_cli.require_subcommand(1);
    this->_cli.add_option("-i,--input", this->_config.path_to_input_matrix, "Path to Modle's contact matrix")->check(CLI::ExistingFile)->required();
    this->_cli.add_option("-o,--output-dir", this->_config.out_dir, "Path where to save the output files.")->required();
    this->_cli.add_option("--tmp-dir", this->_config.tmp_dir, "Path where to store temporary files.")->capture_default_str();
    this->_cli.add_flag("--keep-temporary-files", this->_config.keep_tmp_files, "Do not delete temporary files.")->capture_default_str();
    this->_cli.add_flag("-f,--force", this->_config.force, "Overwrite existing file(s).")->capture_default_str();

    auto convert = this->_cli.add_subcommand("convert", "Convert Modle's contact matrix into several popular formats.")->fallthrough();
    auto eval = this->_cli.add_subcommand("evaluate", "Compare Model's output with other contact matrices using various correlation tests.")->fallthrough();
    eval->alias("eval");

    convert->add_flag("--hic,!--no-hic", this->_config.convert_to_hic, "Convert contact matrix to .hic format (requires Juicer Tools).")->capture_default_str();
    convert->add_flag("--tsv,!--no-tsv", this->_config.convert_to_tsv, "Convert contact matrix to TSV format.")->capture_default_str();
    convert->add_flag("--compress", this->_config.compress, "Compress output using bzip2.")->capture_default_str();
    convert->add_flag("--add-noise,--make-realistic", this->_config.add_noise, "Add noise to make the contact matrix more similar to the matrices produced by Hi-C experiments.")->capture_default_str();
    convert->add_option("--noise-stdev", this->_config.noise_stdev, "Standard deviation to use when sampling noise to add to contact matrices.")->check(CLI::PositiveNumber)->needs(convert->get_option("--add-noise"))->capture_default_str();
    convert->add_option("--noise-range", this->_config.noise_range, "Range to use when adding noise to contact matrices. 99.9%% of the sampled numbers will fall within this range.")->check(CLI::PositiveNumber)->needs(convert->get_option("--add-noise"))->excludes(convert->get_option("--noise-stdev"))->capture_default_str();
    convert->add_option("-j,--path-to-juicer-tools", this->_config.path_to_juicer_tools, "Path to Juicer tools jar. When this is not specified, we will look in PATH for an executable called juicer_tools, juicer_tools.exe or juicer_tools.sh.")->check(CLI::ExistingFile);
    convert->add_option("--juicer-tools-max-mem", this->_config.juicer_tools_mem, "Maximum memory allocation pool for Juicer Tools (JVM).")->needs(convert->get_option("--path-to-juicer-tools"))->transform(CLI::AsSizeValue(false))->capture_default_str();
    convert->add_option("-c,--chr-sizes", this->_config.chr_sizes, absl::StrFormat("Path to file containing chromosome size(s). Can also be one of the following genome IDs: %s.", absl::StrJoin(this->_allowed_genome_ids, ", ")));
    convert->add_option("--seed", this->_config.seed, "Seed value to use when adding noise to the contact matrix.")->capture_default_str();

    eval->add_flag("--spearman", this->_config.compute_spearman, "Compute Spearman rank correlation.");
    eval->add_flag("--kendall", this->_config.compute_kendall, "Compute Kendall rank correlation.");
    eval->add_option("-w,--sliding-window-size", this->_config.sliding_window_size, "Sliding window size.")->check(CLI::PositiveNumber)->capture_default_str();
    eval->add_option("--sliding-window-overlap", this->_config.sliding_window_overlap, "Overlap between consecutive sliding-windows.")->check(CLI::PositiveNumber)->capture_default_str();
    // clang-format on
  }

  [[nodiscard]] bool validate() {
    std::string errors;
    auto& c = this->_config;
    if (!std::filesystem::is_directory(c.out_dir) && std::filesystem::exists(c.out_dir)) {
      absl::StrAppendFormat(
          &errors, "--output-dir should point to a directory or a non-existing path. Is '%s'",
          c.out_dir);
    }

    if (c.path_to_input_matrix.rfind(".tsv.bz2") == std::string::npos) {
      absl::StrAppend(&errors, "The matrix passed to --input should be in .tsv.bz2 format.");
    }
    c.base_name = std::filesystem::path(
                      c.path_to_input_matrix.substr(0, c.path_to_input_matrix.rfind(".tsv.bz2")))
                      .filename();

    if (auto s = "/modle_tools/"; !absl::EndsWithIgnoreCase(c.tmp_dir, s)) c.tmp_dir += s;
    if (!std::filesystem::is_directory(c.tmp_dir) && std::filesystem::exists(c.tmp_dir)) {
      absl::StrAppendFormat(&errors,
                            "--tmp-dir should point to a directory or a non-existing path. Is '%s'",
                            c.tmp_dir);
    }

    if (!c.force) {
      auto base_name = absl::StrFormat("%s/%s", c.out_dir,
                                       std::string_view(c.path_to_input_matrix.data(),
                                                        c.path_to_input_matrix.rfind(".tsv.bz2")));
      std::vector<std::string> collisions;
      if (this->_cli.get_subcommand("convert")->parsed()) {
        if (this->_config.convert_to_hic) {
          if (auto file =
                  absl::StrFormat("%s%s.hic", this->_config.add_noise ? "_w_noise" : "", base_name);
              std::filesystem::exists(file)) {
            collisions.emplace_back(std::move(file));
          }
        }
        if (this->_config.convert_to_tsv) {
          if (auto file = absl::StrFormat("%s%s_symmetric.tsv%s", base_name,
                                          this->_config.add_noise ? "_w_noise" : "",
                                          this->_config.compress ? ".bz2" : "");
              std::filesystem::exists(file)) {
            collisions.emplace_back(std::move(file));
          }
        }
      }
      if (this->_cli.get_subcommand("eval")->parsed()) {
        for (std::string_view suffix : {"rho", "tau", "rho_pv", "tau_pv"}) {
          for (std::string_view ext : {"tsv.bz2", "wig"}) {
            if (auto file = absl::StrFormat("%s_%s.%s", c.base_name, suffix, ext);
                std::filesystem::exists(file)) {
              collisions.emplace_back(std::move(file));
            }
          }
        }
      }
      if (!collisions.empty()) {
        absl::StrAppendFormat(
            &errors,
            "Detected %lu file name collisions: refusing to proceed. Pass --force to "
            "overwrite existing file(s).\nColliding file(s):\n - %s",
            collisions.size(), absl::StrJoin(collisions, "\n - "));
      }
    }

    if (this->_cli.get_subcommand("convert")->parsed()) {
      if (std::filesystem::exists(c.chr_sizes)) {
        if (std::filesystem::is_directory(c.chr_sizes)) {
          absl::StrAppendFormat(
              &errors,
              "--chr-sizes should be the path to a file or one of %s, but is a directory\n",
              absl::StrJoin(this->_allowed_genome_ids, ", "));
        }
      } else {
        if (!this->_allowed_genome_ids.contains(c.chr_sizes)) {
          absl::StrAppendFormat(
              &errors, "--chr-sizes='%s' should be the path to an existing file or one of %s.\n",
              c.chr_sizes, absl::StrJoin(this->_allowed_genome_ids, ", "));
        }
      }
    }

    if (this->_cli.get_subcommand("eval")->parsed() &&
        c.sliding_window_size <= c.sliding_window_overlap) {
      absl::StrAppendFormat(&errors,
                            "--sliding-window-size should be > --sliding-window-overlap: "
                            "--sliding-window-size=%lu --sliding-window-overlap=%lu.",
                            c.sliding_window_size, c.sliding_window_overlap);
    }
    if (!errors.empty()) {
      absl::FPrintF(stderr, "The following issues have been detected:\n%s\n", errors);
    }
    return errors.empty();
  }

 public:
  Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(argv[0]) { this->MakeCli(); }
  [[nodiscard]] bool is_ok() const {
    return this->_config.exit_code && this->_subcommand != subcommand::help;
  }
  [[nodiscard]] subcommand get_subcommand() const { return this->_subcommand; }
  [[nodiscard]] config parse_arguments() {
    this->_exec_name = this->_argv[0];
    try {
      this->_cli.parse(this->_argc, this->_argv);
      if (this->_cli.get_subcommand("convert")->parsed())
        this->_subcommand = subcommand::convert;
      else if (this->_cli.get_subcommand("evaluate")->parsed())
        this->_subcommand = subcommand::eval;
      else
        this->_subcommand = subcommand::help;
    } catch (const CLI::ParseError& e) {
      //  This takes care of formatting and printing error messages (if any)
      this->_config.exit_code = this->_cli.exit(e);
      return this->_config;
    }
    this->_config.exit_code = this->validate();
    return this->_config;
  }
};

}  // namespace modle::utils