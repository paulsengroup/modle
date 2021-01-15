#pragma once
#include <absl/container/btree_set.h>
#include <absl/strings/match.h>
#include <absl/strings/str_format.h>
#include <absl/strings/str_join.h>
#include <absl/strings/strip.h>
#include <fmt/printf.h>

#include <CLI/CLI.hpp>
#include <boost/process.hpp>
#include <cstdint>
#include <string>

#include "modle/cooler.hpp"
#include "modle_tools/config.hpp"

namespace modle::tools {

class Cli {
 public:
  enum subcommand : uint8_t { convert, eval, stats, help };

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

  inline void MakeCli() {
    // clang-format off
    this->_cli.name(this->_exec_name);
    this->_cli.description("Modle's helper tool. This tool allows to post-process the contact matrix produced by Modle in various ways.");
    this->_cli.require_subcommand(1);
    this->_cli.add_option("-i,--input", this->_config.path_to_input_matrices, "Path to a contact matrix in Cooler format")->check(CLI::ExistingFile)->required();
    this->_cli.add_option("--tmp-dir", this->_config.tmp_dir, "Path where to store temporary files.")->capture_default_str();
    this->_cli.add_flag("--keep-temporary-files", this->_config.keep_tmp_files, "Do not delete temporary files.")->capture_default_str();
    this->_cli.add_flag("-f,--force", this->_config.force, "Overwrite existing file(s).")->capture_default_str();
    this->_cli.add_option("-j,--path-to-juicer-tools", this->_config.path_to_juicer_tools, "Path to Juicer tools jar. When this is not specified, we will look in PATH for an executable called juicer_tools, juicer_tools.exe or juicer_tools.sh.")->check(CLI::ExistingFile);
    this->_cli.add_option("-t,--threads", this->_config.nthreads, "CPU threads to allocate.")->check(CLI::PositiveNumber);

    auto *convert_sc = this->_cli.add_subcommand("convert", "Convert Modle's contact matrix into several popular formats.")->fallthrough();
    auto *eval_sc = this->_cli.add_subcommand("evaluate", "Compare Model's output with other contact matrices using various correlation tests.")->fallthrough();
    auto *stats_sc = this->_cli.add_subcommand("statistics", "Compute several useful statistics for a given Cooler file.")->fallthrough();
    eval_sc->alias("eval");
    stats_sc->alias("stats");

    convert_sc->add_option("-o,--output-base-name", this->_config.output_base_name, "Base file name (including directories) to use for output.")->required();
    convert_sc->add_flag("--output-format", this->_config.output_format, "Output file format. Accepted formats are: HIC and TSV.")->capture_default_str();
    convert_sc->add_flag("--add-noise,--make-realistic", this->_config.add_noise, "Add noise to make the contact matrix more similar to the matrices produced by Hi-C experiments.")->capture_default_str();
    convert_sc->add_option("--noise-stdev", this->_config.noise_stdev, "Standard deviation to use when sampling noise to add to contact matrices.")->check(CLI::PositiveNumber)->needs(convert_sc->get_option("--add-noise"));
    convert_sc->add_option("--noise-range", this->_config.noise_range, "Range to use when adding noise to contact matrices. By default 99.9% of the sampled numbers will fall within this range.")->check(CLI::PositiveNumber)->needs(convert_sc->get_option("--add-noise"))->excludes(convert_sc->get_option("--noise-stdev"))->capture_default_str();
    convert_sc->add_option("--juicer-tools-max-mem", this->_config.juicer_tools_mem, "Maximum memory allocation pool for Juicer Tools (JVM).")->needs(this->_cli.get_option("--path-to-juicer-tools"))->transform(CLI::AsSizeValue(false))->capture_default_str();
    convert_sc->add_option("-c,--chr-sizes", this->_config.chr_sizes, fmt::format("Path to file containing chromosome size(s). Can also be one of the following genome IDs: {}.", absl::StrJoin(this->_allowed_genome_ids, ", ")));
    convert_sc->add_option("--seed", this->_config.seed, "Seed value to use when adding noise to the contact matrix.")->capture_default_str();

    eval_sc->add_option("-o,--output-base-name", this->_config.output_base_name, "Base file name (including directories) to use for output.")->required();
    eval_sc->add_option("--reference-matrix", this->_config.path_to_reference_matrix, "Path to contact matrix to use as reference when computing the correlation. Formats accepted: Cooler.")->required();
    eval_sc->add_option("-n,--chr-name", this->_config.chr_name_hic, "Name of the chromosome whose contacts should be extracted from an .hic file. Required only if different from the name stored in the header of Modle's contact matrix.");
    eval_sc->add_option("--chr-offset", this->_config.chr_offset_hic, "Offset to apply to coordinates read from Modle's contact matrix.")->check(CLI::NonNegativeNumber)->transform(CLI::AsSizeValue(true));
    eval_sc->add_flag("--pearson", this->_config.compute_pearson, "Compute Pearson correlation.");
    eval_sc->add_flag("--spearman", this->_config.compute_spearman, "Compute Spearman rank correlation.");
    eval_sc->add_option("-w,--diagonal-width", this->_config.diagonal_width, "Diagonal width to use when computing correlation coefficients.")->check(CLI::NonNegativeNumber)->required();
    eval_sc->add_option("--sliding-window-size", this->_config.sliding_window_size, "Sliding window size. By default this is set to the diagonal width of ModLE's contact matrix.")->check(CLI::NonNegativeNumber);
    eval_sc->add_option("--sliding-window-overlap", this->_config.sliding_window_overlap, "Overlap between consecutive sliding-windows.")->check(CLI::NonNegativeNumber)->capture_default_str();

    stats_sc->add_option("--bin-size", this->_config.bin_size, "Bin size to use when calculating the statistics. Required in case of MCool files.")->check(CLI::PositiveNumber);
    stats_sc->add_option("-w,--diagonal-width", this->_config.diagonal_width, "Diagonal width.")->check(CLI::PositiveNumber)->required();
    stats_sc->add_flag("--dump-depleted-matrices", this->_config.dump_depleted_matrices, "Dump contact matrices used to calculate the normalized average contact density.");
    stats_sc->add_option("--path-to-histograms", this->_config.output_path_for_histograms, "Path where to output contact histograms.");
    stats_sc->add_option("--exclude-chromosomes", this->_config.chromosomes_excluded_vect, "Comma-separated list of chromosomes to skip when calculating chromosome statistics.")->delimiter(',');
    stats_sc->add_option("--depletion-multiplier", this->_config.depletion_multiplier, "Multiplier used to control the magnitude of the depletion.")->check(CLI::NonNegativeNumber)->capture_default_str();
    // clang-format on
  }

  [[nodiscard]] inline bool validate() {  // TODO: Refactor this function
    std::string errors;
    auto& c = this->_config;
    if (!std::filesystem::is_directory(c.output_base_name) &&
        std::filesystem::exists(c.output_base_name)) {
      absl::StrAppendFormat(
          &errors, "--output-dir should point to a directory or a non-existing path. Is '%s'",
          c.output_base_name);
    }

    if (const auto* s = "/modle_tools/"; !absl::EndsWithIgnoreCase(c.tmp_dir, s)) {
      c.tmp_dir += s;
    }
    if (!std::filesystem::is_directory(c.tmp_dir) && std::filesystem::exists(c.tmp_dir)) {
      absl::StrAppendFormat(
          &errors, "--tmp-dir should point to a directory or a non-existing path. Is '%s'\n",
          c.tmp_dir);
    }

    if (!c.force) {
      std::vector<std::string> collisions;
      if (this->_cli.get_subcommand("convert")->parsed()) {
        if (this->_config.output_format == config::output_format::hic) {
          if (auto file =
                  fmt::format("{}{}.hic", c.add_noise ? "_w_noise" : "", c.output_base_name);
              std::filesystem::exists(file)) {
            collisions.push_back(file);
          }
        }
        if (this->_config.output_format == config::output_format::tsv) {
          if (auto file = fmt::format("{}{}_symmetric.tsv.bz2", c.output_base_name,
                                      c.add_noise ? "_w_noise" : "");
              std::filesystem::exists(file)) {
            collisions.push_back(file);
          }
        }
      }
      if (this->_cli.get_subcommand("eval")->parsed()) {
        if (c.path_to_input_matrices.size() > 1) {
          absl::StrAppendFormat(&errors,
                                "You have specified %d input files through the -i option, but "
                                "only one file is allowed when using the eval subcommand.\n",
                                c.path_to_input_matrices.size());
        }

        for (std::string_view suffix :
             {"rho", "tau", "rho_pv", "tau_pv"}) {  // TODO Update this section
          for (std::string_view ext : {"tsv.bz2", "wig"}) {
            if (auto file = fmt::format("{}_{}.{}", c.output_base_name, suffix, ext);
                std::filesystem::exists(file)) {
              collisions.emplace_back(file);
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

    if (!c.chr_sizes.empty()) {
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

    if (this->_cli.get_subcommand("eval")->parsed() && c.sliding_window_size != 0 &&
        c.sliding_window_size <= c.sliding_window_overlap) {
      absl::StrAppendFormat(&errors,
                            "--sliding-window-size should be > --sliding-window-overlap: "
                            "--sliding-window-size=%lu --sliding-window-overlap=%lu.",
                            c.sliding_window_size, c.sliding_window_overlap);
    }
    if (this->_config.output_format == config::output_format::hic ||
        absl::StartsWith("http", this->_config.path_to_reference_matrix) ||
        absl::EndsWithIgnoreCase(c.path_to_reference_matrix, ".hic")) {
      if (c.path_to_juicer_tools.empty()) {
        for (std::string file : {"juicer_tools", "juicer_tools.sh", "juicer_tools.exe"}) {
          if (auto p = boost::process::search_path(file).string(); !p.empty()) {
            c.path_to_juicer_tools = std::move(p);
            break;
          }
        }
      }
      if (c.output_format == config::output_format::hic && c.path_to_juicer_tools.empty()) {
        absl::StrAppendFormat(
            &errors,
            "--path-to-juicer-tools was not specified and ModLE tools was unable to find Juicer "
            "tools in your path");
      }
    }

    if (this->_cli.get_subcommand("stats")->parsed()) {  // Checks whether the bin size is available
      // TODO catch error and print a friendly error message
      cooler::Cooler f(c.path_to_input_matrices[0], cooler::Cooler::READ_ONLY, c.bin_size);
    }

    this->validate_stats_args(errors);

    if (!errors.empty()) {
      fmt::print(stderr, "The following issues have been detected:\n{}\n", errors);
    }
    return errors.empty();
  }

  inline void validate_stats_args(std::string& errors) {
    if (!this->_cli.get_subcommand("stats")->parsed()) {
      return;
    }
    const auto& c = this->_config;
    if (c.path_to_input_matrices.size() != 1) {
      absl::StrAppendFormat(&errors,
                            "modle_tools stats expects one contact matrix, but %d were found.\n",
                            c.path_to_input_matrices.size());
    }

    try {
      cooler::Cooler f(c.path_to_input_matrices[0], cooler::Cooler::READ_ONLY, c.bin_size);
    } catch (const std::runtime_error& e) {
      if (absl::EndsWith(e.what(),
                         "A bin size other than 0 is required when calling "
                         "Cooler::validate_multires_cool_flavor()")) {
        absl::StrAppendFormat(&errors,
                              "File '%s' appears to be a multi-resolution Cooler. --bin-size is a "
                              "mandatory argument when processing .mcool files.\n",
                              c.path_to_input_matrices[0]);
      }
    }

    if (!c.force) {
      const auto ext = std::filesystem::path(c.path_to_input_matrices[0]).extension().string();
      const auto path =
          absl::StrCat(absl::StripSuffix(c.path_to_input_matrices[0], ext), "_depl.cool");
      if (std::filesystem::exists(path)) {
        absl::StrAppendFormat(&errors, "File '%s' already exists. Pass --force to overwrite", path);
      }
      if (std::filesystem::exists(c.output_path_for_histograms)) {
        absl::StrAppendFormat(&errors, "File '%s' already exists. Pass --force to overwrite",
                              c.output_path_for_histograms);
      }
    }
  }

  inline void post_process_cli_args() {
    this->_config.chromosomes_excluded =
        absl::flat_hash_set<std::string>(this->_config.chromosomes_excluded_vect.begin(),
                                         this->_config.chromosomes_excluded_vect.end());
    std::vector<std::string> v;
    std::swap(this->_config.chromosomes_excluded_vect, v);
  }

 public:
  inline Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) {
    this->MakeCli();
  }
  [[nodiscard]] inline bool is_ok() const {
    return this->_config.exit_code && this->_subcommand != subcommand::help;
  }
  [[nodiscard]] inline subcommand get_subcommand() const { return this->_subcommand; }
  [[nodiscard]] inline config parse_arguments() {
    this->_exec_name = *this->_argv;
    try {
      this->_cli.parse(this->_argc, this->_argv);
      if (this->_cli.get_subcommand("convert")->parsed()) {
        this->_subcommand = subcommand::convert;
      } else if (this->_cli.get_subcommand("evaluate")->parsed()) {
        this->_subcommand = subcommand::eval;
      } else if (this->_cli.get_subcommand("statistics")->parsed()) {
        this->_subcommand = subcommand::stats;
      } else {
        this->_subcommand = subcommand::help;
      }
    } catch (const CLI::ParseError& e) {
      //  This takes care of formatting and printing error messages (if any)
      this->_config.exit_code = this->_cli.exit(e);
      return this->_config;
    }
    this->_config.exit_code = this->validate();
    this->post_process_cli_args();
    return this->_config;
  }
};

}  // namespace modle::tools
