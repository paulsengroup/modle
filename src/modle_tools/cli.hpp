#pragma once

#include <absl/container/btree_set.h>      // for btree_set
#include <absl/container/flat_hash_set.h>  // for BitMask, flat_hash_set, raw_hash_set
#include <absl/strings/match.h>            // for EndsWith, EndsWithIgnoreCase
#include <absl/strings/str_cat.h>          // for StrCat
#include <absl/strings/str_format.h>       // for StrAppendFormat
#include <absl/strings/str_join.h>         // for StrJoin
#include <absl/strings/strip.h>            // for StripSuffix
#include <fmt/format.h>                    // for format, print
#include <fmt/ostream.h>                   // for formatbuf<>::int_type

#include <CLI/CLI.hpp>       // IWYU pragma: keep
#include <algorithm>         // for max
#include <cassert>           // for assert
#include <cstdint>           // for uint32_t, uint_fast8_t
#include <cstdio>            // for stderr
#include <filesystem>        // for path, exists, is_directory, operat...
#include <initializer_list>  // for initializer_list
#include <sstream>           // for size_t, basic_stringbuf<>::int_type
#include <stdexcept>         // for invalid_argument, out_of_range
#include <string>            // for string, allocator, basic_string
#include <string_view>       // for string_view
#include <vector>            // for vector, swap

#include "modle/bed.hpp"           // for Parser
#include "modle/common/utils.hpp"  // for throw_with_trace
#include "modle/cooler.hpp"        // for Cooler, Cooler::READ_ONLY
#include "modle_tools/config.hpp"  // for config

namespace modle::tools {

class Cli {
 public:
  enum subcommand : uint_fast8_t { eval, noisify, stats, help };

 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  int _exit_code{1};
  config _config{};
  CLI::App _cli{};
  subcommand _subcommand{subcommand::help};

  inline void make_eval_subcommand() {
    // clang-format off
    auto *sc = this->_cli.add_subcommand("evaluate", "Compare Model's output with other contact matrices using various correlation tests.")->fallthrough();
    sc->alias("eval");

    sc->add_option("-i,--input", this->_config.path_to_input_matrix, "Path to a contact matrix in Cooler format")->check(CLI::ExistingFile)->required();
    sc->add_option("--chromosome-subrange-file", this->_config.path_to_chrom_subranges, "Path to BED file with subranges of the chromosomes to be processed.")->check(CLI::ExistingFile);
    sc->add_option("--tmp-dir", this->_config.tmp_dir, "Path where to store temporary files.")->capture_default_str();
    sc->add_flag("--keep-temporary-files", this->_config.keep_tmp_files, "Do not delete temporary files.")->capture_default_str();
    sc->add_flag("-f,--force", this->_config.force, "Overwrite existing file(s).")->capture_default_str();
    sc->add_option("-t,--threads", this->_config.nthreads, "CPU threads to allocate.")->check(CLI::PositiveNumber);
    sc->add_option("-o,--output-base-name", this->_config.output_base_name, "Base file name (including directories) to use for output.")->required();
    sc->add_option("--reference-matrix", this->_config.path_to_reference_matrix, "Path to contact matrix to use as reference when computing the correlation. Formats accepted: Cooler.")->required();
    sc->add_flag("--eucl-dist", this->_config.compute_edist, "Compute Euclidean distance.");
    sc->add_flag("--pearson", this->_config.compute_pearson, "Compute Pearson correlation.");
    sc->add_flag("--spearman", this->_config.compute_spearman, "Compute Spearman rank correlation.");
    sc->add_option("-w,--diagonal-width", this->_config.diagonal_width, "Diagonal width to use when computing correlation coefficients.")->check(CLI::NonNegativeNumber)->required();
    sc->add_option("--sliding-window-size", this->_config.sliding_window_size, "Sliding window size. By default this is set to the diagonal width of ModLE's contact matrix.")->check(CLI::NonNegativeNumber);
    sc->add_option("--sliding-window-overlap", this->_config.sliding_window_overlap, "Overlap between consecutive sliding-windows.")->check(CLI::NonNegativeNumber)->capture_default_str();
    sc->add_option("--depletion-multiplier", this->_config.depletion_multiplier, "Multiplier used to control the magnitude of the depletion.")->check(CLI::NonNegativeNumber)->capture_default_str();
    sc->add_flag("--deplete-reference-contacts,!--no-deplete-reference-contacts", this->_config.deplete_contacts_from_reference, "Deplete contacts along the diagonal from the reference matrix.")->capture_default_str();
    // clang-format on
  }

  inline void make_noisify_subcommand() {
    // clang-format off
    auto *sc = this->_cli.add_subcommand("noisify", "Add noise to ModLE's contact matrix in Cooler format.")->fallthrough();
    sc->add_option("-i,--input", this->_config.path_to_input_matrix, "Path to a contact matrix in Cooler format")->check(CLI::ExistingFile)->required();
    sc->add_option("-o,--output-base-name", this->_config.path_to_output_matrix, "Output file name (possibly including directories) to use for output.")->required();
    sc->add_option("-w,--diagonal-width", this->_config.diagonal_width, "Diagonal width of the input contact matrix.")->check(CLI::NonNegativeNumber)->required();
    sc->add_flag("-f,--force", this->_config.force, "Overwrite existing file(s).")->capture_default_str();
    sc->add_option("--mu,--location", this->_config.genextreme_mu, "Location parameter (mu) of the generalized extreme value used to add noise to the contact matrix.")->check(CLI::NonNegativeNumber)->capture_default_str();
    sc->add_option("--sigma,--scale", this->_config.genextreme_sigma, "Scale parameter (sigma) of the generalized extreme value used to add noise to the contact matrix.")->check(CLI::NonNegativeNumber)->capture_default_str();
    sc->add_option("--xi,--shape", this->_config.genextreme_xi, "Shape parameter (xi) of the generalized extreme value used to add noise to the contact matrix.")->check(CLI::NonNegativeNumber)->capture_default_str();
    sc->add_option("--seed", this->_config.seed, "Seed used to initialize the PRNG.")->check(CLI::NonNegativeNumber)->capture_default_str();
    // clang-format on
  }

  inline void make_stats_subcommand() {
    // clang-format off
    auto *sc = this->_cli.add_subcommand("statistics", "Compute several useful statistics for a given Cooler file.")->fallthrough();
    sc->alias("stats");
    sc->add_option("-i,--input", this->_config.path_to_input_matrix, "Path to a contact matrix in Cooler format")->check(CLI::ExistingFile)->required();
    sc->add_option("--chromosome-subrange-file", this->_config.path_to_chrom_subranges, "Path to BED file with subranges of the chromosomes to be processed.")->check(CLI::ExistingFile);
    sc->add_flag("-f,--force", this->_config.force, "Overwrite existing file(s).")->capture_default_str();
    sc->add_option("--bin-size", this->_config.bin_size, "Bin size to use when calculating the statistics. Required in case of MCool files.")->check(CLI::PositiveNumber);
    sc->add_option("-w,--diagonal-width", this->_config.diagonal_width, "Diagonal width.")->check(CLI::PositiveNumber)->required();
    sc->add_flag("--dump-depleted-matrices", this->_config.dump_depleted_matrices, "Dump contact matrices used to calculate the normalized average contact density.");
    sc->add_option("--path-to-histograms", this->_config.output_path_for_histograms, "Path where to output contact histograms.");
    sc->add_option("--exclude-chromosomes", this->_config.chromosomes_excluded_vect, "Comma-separated list of chromosomes to skip when calculating chromosome statistics.")->delimiter(',');
    sc->add_option("--depletion-multiplier", this->_config.depletion_multiplier, "Multiplier used to control the magnitude of the depletion.")->check(CLI::NonNegativeNumber)->capture_default_str();
    // clang-format on
  }

  inline void MakeCli() {
    // clang-format off
    this->_cli.name(this->_exec_name);
    this->_cli.description("Modle's helper tool. This tool allows to post-process the contact matrix produced by Modle in various ways.");
    this->_cli.require_subcommand(1);
    // clang-format on
    this->make_eval_subcommand();
    this->make_noisify_subcommand();
    this->make_stats_subcommand();
  }

  [[nodiscard]] inline std::string validate_eval_subcommand() {
    assert(this->_cli.get_subcommand("eval")->parsed());  // NOLINT
    std::string errors;
    auto& c = this->_config;
    if (!std::filesystem::is_directory(c.output_base_name) &&
        std::filesystem::exists(c.output_base_name)) {
      absl::StrAppendFormat(
          &errors, "--output-dir should point to a directory or a non-existing path. Is '%s'",
          c.output_base_name);
    }

    if (const auto* s = "/modle_tools/"; !absl::EndsWithIgnoreCase(c.tmp_dir.string(), s)) {
      c.tmp_dir += s;
    }
    if (!std::filesystem::is_directory(c.tmp_dir) && std::filesystem::exists(c.tmp_dir)) {
      absl::StrAppendFormat(
          &errors, "--tmp-dir should point to a directory or a non-existing path. Is '%s'\n",
          c.tmp_dir);
    }

    if (!c.force) {
      std::vector<std::string> collisions;
      for (std::string_view suffix :
           {"rho", "tau", "rho_pv", "tau_pv"}) {  // TODO Update this section
        for (std::string_view ext : {"tsv.bz2", "bwig"}) {
          if (auto file = fmt::format("{}_{}.{}", c.output_base_name, suffix, ext);
              std::filesystem::exists(file)) {
            collisions.emplace_back(file);
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

    if (!c.chrom_sizes.empty()) {
      if (std::filesystem::exists(c.chrom_sizes)) {
        if (std::filesystem::is_directory(c.chrom_sizes)) {
          absl::StrAppendFormat(
              &errors,
              "--chrom-sizes='%s' should be the path to an existing chrom.sizes file, "
              "but is a directory.\n",
              c.chrom_sizes);
        }
      } else {
        absl::StrAppendFormat(
            &errors, "--chrom-sizes='%s' should be the path to an existing file or one of.\n",
            c.chrom_sizes);
      }
    }

    if (c.sliding_window_size != 0 && c.sliding_window_size <= c.sliding_window_overlap) {
      absl::StrAppendFormat(&errors,
                            "--sliding-window-size should be > --sliding-window-overlap: "
                            "--sliding-window-size=%lu --sliding-window-overlap=%lu.",
                            c.sliding_window_size, c.sliding_window_overlap);
    }
    return errors;
  }

  [[nodiscard]] inline std::string validate_stats_subcommand() const {
    std::string errors;
    assert(this->_cli.get_subcommand("stats")->parsed());  // NOLINT
    const auto& c = this->_config;

    if (!c.path_to_chrom_subranges.empty()) {
      auto p = modle::bed::Parser(c.path_to_chrom_subranges);
      if (const auto s = p.validate(); !s.empty()) {
        absl::StrAppendFormat(&errors,
                              "Validation of file '%s' failed with the following error: %s.\n",
                              c.path_to_input_matrix, s);
      }
    }

    try {
      cooler::Cooler f(c.path_to_input_matrix, cooler::Cooler::READ_ONLY, c.bin_size);
    } catch (const std::runtime_error& e) {
      if (absl::EndsWith(e.what(),  // NOLINT
                         "A bin size other than 0 is required when calling "
                         "Cooler::validate_multires_cool_flavor()")) {
        absl::StrAppendFormat(&errors,
                              "File '%s' appears to be a multi-resolution Cooler. --bin-size is a "
                              "mandatory argument when processing .mcool files.\n",
                              c.path_to_input_matrix);
      } else {
        absl::StrAppendFormat(&errors,
                              "Validation of file '%s' failed with the following error: %s.\n",
                              c.path_to_input_matrix, e.what());
      }
    }

    if (!c.force) {
      const auto ext = c.path_to_input_matrix.extension().string();
      const auto path =
          absl::StrCat(absl::StripSuffix(c.path_to_input_matrix.string(), ext), "_depl.cool");
      if (std::filesystem::exists(path)) {
        absl::StrAppendFormat(&errors, "File '%s' already exists. Pass --force to overwrite", path);
      }
      if (std::filesystem::exists(c.output_path_for_histograms)) {
        absl::StrAppendFormat(&errors, "File '%s' already exists. Pass --force to overwrite",
                              c.output_path_for_histograms);
      }
    }
    return errors;
  }

  [[nodiscard]] inline bool validate() {  // TODO: Refactor this function
    std::string errors;
    if (this->_cli.get_subcommand("eval")->parsed()) {
      errors = this->validate_eval_subcommand();
    } else if (this->_cli.get_subcommand("noisify")->parsed()) {
      // errors = this->validate_noisify_subcommand();
      (void)0;
    } else if (this->_cli.get_subcommand("stats")->parsed()) {
      errors = this->validate_stats_subcommand();
    } else {
      modle::utils::throw_with_trace(std::logic_error("Unreachable code"));
    }

    if (!errors.empty()) {
      fmt::print(stderr, "The following issues have been detected:\n{}\n", errors);
    }
    return errors.empty();
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
    return (this->_exit_code != 0) && this->_subcommand != subcommand::help;
  }
  [[nodiscard]] inline subcommand get_subcommand() const { return this->_subcommand; }
  [[nodiscard]] inline config parse_arguments() {
    this->_exec_name = *this->_argv;
    try {
      this->_cli.parse(this->_argc, this->_argv);
      if (this->_cli.get_subcommand("evaluate")->parsed()) {
        this->_subcommand = subcommand::eval;
      } else if (this->_cli.get_subcommand("noisify")->parsed()) {
        this->_subcommand = subcommand::noisify;
      } else if (this->_cli.get_subcommand("statistics")->parsed()) {
        this->_subcommand = subcommand::stats;
      } else {
        this->_subcommand = subcommand::help;
      }
    } catch (const CLI::ParseError& e) {
      //  This takes care of formatting and printing error messages (if any)
      this->_exit_code = this->_cli.exit(e);
      return this->_config;
    }
    this->_exit_code = static_cast<int>(this->validate());
    this->post_process_cli_args();
    return this->_config;
  }

  [[nodiscard]] inline int get_exit_code() const { return this->_exit_code; }
};

}  // namespace modle::tools
