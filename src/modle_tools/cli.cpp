#include "./cli.hpp"  // for Cli

#include <absl/container/btree_set.h>      // for btree_set
#include <absl/container/flat_hash_map.h>  // for flat_hash_map
#include <absl/container/flat_hash_set.h>  // for operator!=, BitMask, flat_hash_set
#include <absl/strings/match.h>            // for EndsWith, EndsWithIgnoreCase
#include <absl/strings/str_cat.h>          // for StrCat
#include <absl/strings/str_format.h>       // for StrAppendFormat
#include <absl/strings/str_join.h>         // for StrJoin
#include <absl/strings/strip.h>            // for StripPrefix, StripSuffix
#include <fmt/format.h>                    // for format, FMT_STRING, print
#include <fmt/ostream.h>                   // for formatbuf<>::int_type

#include <CLI/CLI.hpp>
#include <boost/filesystem/path.hpp>  // for path, exists, is_directory, operat...
#include <cassert>                    // for assert
#include <cstdint>                    // for uint64_t, uint32_t, uint_fast8_t
#include <cstdio>                     // for size_t, stderr
#include <exception>                  // for exception
#include <initializer_list>           // for initializer_list
#include <stdexcept>                  // for invalid_argument, out_of_range
#include <string>                     // for string, basic_string, allocator
#include <string_view>                // for string_view
#include <type_traits>                // for remove_reference<>::type
#include <vector>                     // for vector, swap

#include "modle/bed.hpp"           // for BED::Dialect, Parser, bed_dialects
#include "modle/common/utils.hpp"  // for throw_with_trace
#include "modle/cooler.hpp"        // for Cooler, Cooler::READ_ONLY
#include "modle_tools/config.hpp"  // for config

namespace modle::tools {

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { this->make_cli(); }

void Cli::make_eval_subcommand() {
  auto& sc = *this->_cli
                  .add_subcommand("evaluate",
                                  "Compare MoDLE's output with other contact matrices using "
                                  "various correlation tests.")
                  ->fallthrough()
                  ->preparse_callback([this](size_t i) {
                    assert(this->_config.index() == 0);  // NOLINT empty variant
                    fmt::print("callback_arg={}\n", i);
                    this->_config = eval_config{};
                  });
  sc.alias("eval");

  this->_config = eval_config{};
  auto& c = absl::get<eval_config>(this->_config);

  auto& io = *sc.add_option_group("Input/Output", "");
  auto& corr = *sc.add_option_group("Correlation metrics", "");
  auto& ref = *sc.add_option_group("Reference contact matrix", "");
  auto& gen = *sc.add_option_group("Generic", "");

  // clang-format off
  io.add_option(
     "-i,--input",
     c.path_to_input_matrix,
     "Path to a contact matrix in Cooler format.")
     ->check(CLI::ExistingFile)
     ->required();

  io.add_option(
     "-o,--output-base-name",
     c.output_base_name,
     "Base file name (including directories) to use for output.")
     ->required();

  io.add_option(
     "--reference-matrix",
     c.path_to_reference_matrix,
     "Path to contact matrix to use as reference when computing the correlation.\n"
     "Formats accepted: Cooler.")
     ->required();

  io.add_option(
     "--chromosome-subrange-file",
     c.path_to_chrom_subranges,
     "Path to BED file with subranges of the chromosomes to be processed.")
     ->check(CLI::ExistingFile);

  io.add_option(
     "--tmp-dir",
     c.tmp_dir,
     "Path where to store temporary files.")
     ->capture_default_str();

  io.add_flag(
     "--keep-temporary-files",
     c.keep_tmp_files,
     "Do not delete temporary files.")
     ->capture_default_str();

  io.add_flag(
     "-f,--force",
     c.force,
     "Overwrite existing file(s).")
     ->capture_default_str();

  corr.add_flag(
     "--eucl-dist",
     c.compute_edist,
     "Compute Euclidean distance.");

  corr.add_flag(
     "--pearson",
     c.compute_pearson,
     "Compute Pearson correlation.");

  corr.add_flag(
     "--spearman",
     c.compute_spearman,
     "Compute Spearman rank correlation.");

  ref.add_option("-b,--bin-size",
     c.bin_size,
     "Bin size in base pairs.\n"
     "Only used when the contact matrix passed to --reference-matrix is in .mcool format.")
     ->check(CLI::PositiveNumber)
     ->transform(utils::str_float_to_str_int);

  ref.add_option(
     "-w,--diagonal-width",
     c.diagonal_width,
     "Diagonal width to use when computing correlation coefficients.")
     ->check(CLI::NonNegativeNumber)
     ->transform(utils::str_float_to_str_int)
     ->required();
     
  ref.add_option(
     "--depletion-multiplier",
     c.depletion_multiplier,
     "Multiplier used to control the magnitude of the depletion.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();

  ref.add_flag(
     "--deplete-reference-contacts,!--no-deplete-reference-contacts",
     c.deplete_contacts_from_reference,
     "Deplete contacts along the diagonal from the reference matrix.")
     ->capture_default_str();

  gen.add_option(
     "-t,--threads",
     c.nthreads,
     "CPU threads to allocate.")
     ->check(CLI::PositiveNumber);

  gen.add_option(
     "--sliding-window-size",
     c.sliding_window_size,
     "Sliding window size. By default this is set to the diagonal width of ModLE's contact matrix.")
     ->check(CLI::NonNegativeNumber)
     ->transform(utils::str_float_to_str_int);

  gen.add_option(
     "--sliding-window-overlap",
     c.sliding_window_overlap,
     "Overlap between consecutive sliding-windows.")
     ->check(CLI::NonNegativeNumber)
     ->transform(utils::str_float_to_str_int)
     ->capture_default_str();

  // clang-format on
  this->_config = absl::monostate{};
}

void Cli::make_filter_barriers_subcommand() {
  auto& sc =
      *this->_cli
           .add_subcommand("filter-barriers", "Filter extrusion barriers to be used by ModLE.")
           ->fallthrough()
           ->preparse_callback([this](size_t i) {
             assert(this->_config.index() == 0);  // NOLINT empty variant
             fmt::print("callback_arg={}\n", i);
             this->_config = filter_barrier_config{};
           });

  this->_config = filter_barrier_config{};
  auto& c = absl::get<filter_barrier_config>(this->_config);

  auto& io = *sc.add_option_group("Input/Output", "");
  auto& gen = *sc.add_option_group("Generic", "");

  // clang-format off
  io.add_option(
     "--extrusion-barriers-motif-bed",
     c.path_to_extrusion_barrier_motifs_bed,
     "Path to a BED file containing the list of extrusion barriers to be filtered.")
     ->check(CLI::ExistingFile)
     ->required();

  io.add_option(
     "bedFiles",
     c.path_to_bed_files_for_filtering,
     "One or more path to BED files containing records to be used to select records from the file specified through --extrusion-barriers-motif-bed.")
     ->check(CLI::ExistingFile)
     ->required();

  gen.add_option(
     "-f,--filtering-criterion",
     c.filtering_criterion,
     fmt::format(FMT_STRING("Filtering criterion. Accepted values are: {}"), absl::StrJoin(this->_filtering_criteria, ", ")))
     ->check(
     [this](const auto &str){
     if (this->_filtering_criteria.contains(CLI::detail::to_lower(str))) {
     return std::string{};
     }
     return fmt::format(FMT_STRING("\"{}\" is not a valid filtering criterion. Allowed criteria are: {}."), str, absl::StrJoin(this->_filtering_criteria, ", "));
     })
     ->capture_default_str();

  gen.add_option(
     "--bed-dialect",
     c.bed_dialect,
     fmt::format(FMT_STRING("Specify the BED dialect to use when parsing input files.\n"
     "Example: when specifying BED3 through this dialect, we only validate the first three fields.\n"
     "Additional fields (if any) are copied verbatim and can thus in principle contain arbitrary information.\n"
     "Allowed dialects: {}."), absl::StrJoin(bed::bed_dialects, ", ")))
     ->transform(CLI::CheckedTransformer(bed::str_to_bed_dialect_mappings))
     ->capture_default_str();

  gen.add_flag(
     "--strict-bed-validation,!--no-strict-bed-validation",
     c.strict_bed_validation,
     "Toggle strict BED file format validation on or off.")
     ->capture_default_str();
  // clang-format on
  this->_config = absl::monostate{};
}

void Cli::make_noisify_subcommand() {
  auto& sc =
      *this->_cli
           .add_subcommand("noisify", "Add noise to MoDLE's contact matrix in Cooler format.")
           ->fallthrough()
           ->preparse_callback([this](size_t i) {
             assert(this->_config.index() == 0);  // NOLINT empty variant
             fmt::print("callback_arg={}\n", i);
             this->_config = noisify_config{};
           });

  this->_config = noisify_config{};
  auto& c = absl::get<noisify_config>(this->_config);

  auto& io = *sc.add_option_group("Input/Output", "");
  auto& gen = *sc.add_option_group("Generic", "");
  auto& noise = *sc.add_option_group("Noise properties", "");

  // clang-format off
  io.add_option(
     "-i,--input",
     c.path_to_input_matrix,
     "Path to a contact matrix in Cooler format.")
     ->check(CLI::ExistingFile)
     ->required();

  io.add_option(
     "-o,--output-name",
     c.path_to_output_matrix,
     "Output file name (possibly including directories) to use for output.")
     ->required();

  io.add_flag(
     "-f,--force",
     c.force,
     "Overwrite existing file(s).")
     ->capture_default_str();

  gen.add_option(
     "-w,--diagonal-width",
     c.diagonal_width,
     "Diagonal width of the input contact matrix.")
     ->check(CLI::NonNegativeNumber)
     ->transform(utils::str_float_to_str_int)
     ->required();

  gen.add_option(
      "--bin-size",
      c.bin_size,
      "Bin size of the matrix to noisify. Required and only used when the input matrix is in .mcool format.")
      ->check(CLI::PositiveNumber);

  noise.add_option(
     "--mu,--location",
     c.genextreme_mu,
     "Location parameter (mu) of the generalized extreme value used to.add noise to the contact matrix.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();

  noise.add_option(
     "--sigma,--scale",
     c.genextreme_sigma,
     "Scale parameter (sigma) of the generalized extreme value used to.add noise to the contact matrix.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();

  noise.add_option(
     "--xi,--shape",
     c.genextreme_xi,
     "Shape parameter (xi) of the generalized extreme value used to.add noise to the contact matrix.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();

  noise.add_option(
     "--seed",
     c.seed,
     "Seed used to initialize the PRNG.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();
  // clang-format on
  this->_config = absl::monostate{};
}

void Cli::make_stats_subcommand() {
  auto& sc = *this->_cli
                  .add_subcommand("statistics",
                                  "Compute several useful statistics for a given Cooler file.")
                  ->fallthrough()
                  ->preparse_callback([this](size_t i) {
                    assert(this->_config.index() == 0);  // NOLINT empty variant
                    fmt::print("callback_arg={}\n", i);
                    this->_config = stats_config{};
                  });
  sc.alias("stats");

  this->_config = stats_config{};
  auto& c = absl::get<stats_config>(this->_config);

  auto& io = *sc.add_option_group("Input/Output", "");
  auto& cm = *sc.add_option_group("Contact matrix", "");

  // clang-format off
  io.add_option(
     "-i,--input",
     c.path_to_input_matrix,
     "Path to a contact matrix in Cooler format.")
     ->check(CLI::ExistingFile)
     ->required();

  io.add_option(
     "--chromosome-subrange-file",
     c.path_to_chrom_subranges,
     "Path to BED file with subranges of the chromosomes to be processed.")
     ->check(CLI::ExistingFile);

  io.add_option(
      "--path-to-histograms",
      c.output_path_for_histograms,
      "Path where to output contact histograms.");

  io.add_flag(
      "--dump-depleted-matrices",
      c.dump_depleted_matrices,
      "Dump contact matrices used to calculate the normalized average contact density.");

  io.add_flag(
     "-f,--force",
     c.force,
     "Overwrite existing file(s).")
     ->capture_default_str();

  cm.add_option(
     "--bin-size",
     c.bin_size,
     "Bin size to use when calculating the statistics. Required in case of MCool files.")
     ->check(CLI::PositiveNumber);

  cm.add_option(
     "-w,--diagonal-width",
     c.diagonal_width,
     "Diagonal width.")
     ->check(CLI::PositiveNumber)
     ->transform(utils::str_float_to_str_int)
     ->required();

  cm.add_option(
     "--exclude-chromosomes",
     c.chromosomes_excluded_vect,
     "Comma-separated list of chromosomes to skip when calculating chromosome statistics.")
     ->delimiter(',');

  cm.add_option(
     "--depletion-multiplier",
     c.depletion_multiplier,
     "Multiplier used to control the magnitude of the depletion.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();
  // clang-format on
  this->_config = absl::monostate{};
}
/*
void Cli::make_find_barrier_clusters_subcommand() {
     auto* sc = this->_cli
  .add_subcommand("find-barrier-clusters",
     "Detect clusters of extrusion barriers for a given BED file.")
     ->fallthrough();
     sc.alias("fbcl");

     // clang-format off
  sc.add_option(
     "-i,--input",
     this->_config.path_to_input_matrix,
     "Path to a contact matrix in Cooler format.")
     ->check(CLI::ExistingFile)
     ->required();

  sc.add_option(
     "--chromosome-subrange-file",
     this->_config.path_to_chrom_subranges,
     "Path to BED file with subranges of the chromosomes to be processed.")
     ->check(CLI::ExistingFile);

  sc.add_flag(
     "-f,--force",
     this->_config.force,
     "Overwrite existing file(s).")
     ->capture_default_str();

  sc.add_option(
     "--bin-size",
     this->_config.bin_size,
     "Bin size to use when calculating the statistics. Required in case of MCool files.")
     ->check(CLI::PositiveNumber);

  sc.add_option(
     "-w,--diagonal-width",
     this->_config.diagonal_width,
     "Diagonal width.")
     ->check(CLI::PositiveNumber)
     ->transform(utils::str_float_to_str_int)
     ->required();

  sc.add_flag(
     "--dump-depleted-matrices",
     this->_config.dump_depleted_matrices,
     "Dump contact matrices used to calculate the normalized average contact density.");

  sc.add_option(
     "--path-to-histograms",
     this->_config.output_path_for_histograms,
     "Path where to output contact histograms.");

  sc.add_option(
     "--exclude-chromosomes",
     this->_config.chromosomes_excluded_vect,
     "Comma-separated list of chromosomes to skip when calculating chromosome statistics.")
     ->delimiter(',');

  sc.add_option(
     "--depletion-multiplier",
     this->_config.depletion_multiplier,
     "Multiplier used to control the magnitude of the depletion.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();
     // clang-format on
}
     */

void Cli::make_cli() {
  // clang-format off
     this->_cli.name(this->_exec_name);
     this->_cli.description("Modle's helper tool. This tool allows to post-process the contact matrix produced by Modle in various ways.");
     this->_cli.require_subcommand(1);
  // clang-format on
  this->make_eval_subcommand();
  this->make_filter_barriers_subcommand();
  this->make_noisify_subcommand();
  this->make_stats_subcommand();
}

std::string Cli::validate_eval_subcommand() {
  assert(this->_cli.get_subcommand("eval")->parsed());  // NOLINT
  std::string errors;
  auto& c = absl::get<eval_config>(this->_config);
  if (!boost::filesystem::is_directory(c.output_base_name) &&
      boost::filesystem::exists(c.output_base_name)) {
    absl::StrAppendFormat(
        &errors, "--output-dir should point to a directory or a non-existing path. Is \"%s\"",
        c.output_base_name.string());
  }

  if (c.path_to_reference_matrix.extension() == ".mcool" && c.bin_size == 0) {
    absl::StrAppend(&errors,
                    "--bin-size is required when the contact matrix passed with "
                    "--reference-matrix is in .mcool format.\n");
  }

  if (const auto* s = "/modle_tools/"; !absl::EndsWithIgnoreCase(c.tmp_dir.string(), s)) {
    c.tmp_dir += s;
  }
  if (!boost::filesystem::is_directory(c.tmp_dir) && boost::filesystem::exists(c.tmp_dir)) {
    absl::StrAppendFormat(
        &errors, "--tmp-dir should point to a directory or a non-existing path. Is \"%s\"\n",
        c.tmp_dir.string());
  }

  if (!c.force) {
    std::vector<std::string> collisions;
    for (std::string_view suffix :
         {"rho", "tau", "rho_pv", "tau_pv"}) {  // TODO Update this section
      for (std::string_view ext : {"tsv.bz2", "bwig"}) {
        if (auto file = fmt::format("{}_{}.{}", c.output_base_name, suffix, ext);
            boost::filesystem::exists(file)) {
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

  if (c.sliding_window_size != 0 && c.sliding_window_size <= c.sliding_window_overlap) {
    absl::StrAppendFormat(&errors,
                          "--sliding-window-size should be > --sliding-window-overlap: "
                          "--sliding-window-size=%lu --sliding-window-overlap=%lu.",
                          c.sliding_window_size, c.sliding_window_overlap);
  }
  return errors;
}

std::string Cli::validate_filter_barriers_subcommand() const {
  std::string errors;
  assert(this->_cli.get_subcommand("filter-barriers")->parsed());  // NOLINT
  const auto& c = absl::get<filter_barrier_config>(this->_config);

  auto valitade_bed = [&](const auto& path) {
    if (boost::filesystem::is_regular_file(path)) {
      const auto status = bed::Parser(path, c.bed_dialect, c.strict_bed_validation).validate();
      if (!status.empty()) {
        absl::StrAppendFormat(&errors, "Validation failed for file \"%s\": %s\n", path.string(),
                              absl::StripPrefix(status, "An error occurred while reading file"));
      }
    }
  };

  valitade_bed(c.path_to_extrusion_barrier_motifs_bed);
  assert(!c.path_to_extrusion_barrier_motifs_bed.empty());  // NOLINT
  for (const auto& path : c.path_to_extrusion_barrier_motifs_bed) {
    valitade_bed(path);
  }

  return errors;
}

std::string Cli::validate_stats_subcommand() const {
  std::string errors;
  assert(this->_cli.get_subcommand("stats")->parsed());  // NOLINT
  const auto& c = absl::get<stats_config>(this->_config);

  if (!c.path_to_chrom_subranges.empty()) {
    auto p = modle::bed::Parser(c.path_to_chrom_subranges);
    if (const auto s = p.validate(); !s.empty()) {
      absl::StrAppendFormat(&errors,
                            "Validation of file \"%s\" failed with the following error: %s.\n",
                            c.path_to_input_matrix.string(), s);
    }
  }

  try {
    cooler::Cooler f(c.path_to_input_matrix, cooler::Cooler::READ_ONLY, c.bin_size);
  } catch (const std::runtime_error& e) {
    if (absl::EndsWith(e.what(),  // NOLINT
                       "A bin size other than 0 is required when calling "
                       "Cooler::validate_multires_cool_flavor()")) {
      absl::StrAppendFormat(&errors,
                            "File \"%s\" appears to be a multi-resolution Cooler. --bin-size is a "
                            "mandatory argument when processing .mcool files.\n",
                            c.path_to_input_matrix.string());
    } else {
      absl::StrAppendFormat(&errors,
                            "Validation of file \"%s\" failed with the following error: %s.\n",
                            c.path_to_input_matrix.string(), e.what());
    }
  }

  if (!c.force) {
    const auto ext = c.path_to_input_matrix.extension().string();
    const auto path =
        absl::StrCat(absl::StripSuffix(c.path_to_input_matrix.string(), ext), "_depl.cool");
    if (boost::filesystem::exists(path)) {
      absl::StrAppendFormat(&errors, "File \"%s\" already exists. Pass --force to overwrite", path);
    }
    if (boost::filesystem::exists(c.output_path_for_histograms)) {
      absl::StrAppendFormat(&errors, "File \"%s\" already exists. Pass --force to overwrite",
                            c.output_path_for_histograms);
    }
  }
  return errors;
}

std::string Cli::validate_noisify_subcommand() const {
  std::string errors;
  assert(this->_cli.get_subcommand("noisify")->parsed());  // NOLINT
  const auto& c = absl::get<noisify_config>(this->_config);

  assert(boost::filesystem::exists(c.path_to_input_matrix));  // NOLINT
  try {
    cooler::Cooler f(c.path_to_input_matrix, cooler::Cooler::READ_ONLY, c.bin_size);
  } catch (const std::runtime_error& e) {
    if (absl::EndsWith(e.what(),  // NOLINT
                       "A bin size other than 0 is required when calling "
                       "Cooler::validate_multires_cool_flavor()")) {
      absl::StrAppendFormat(
          &errors,
          "File \"%s\" appears to be a multi-resolution Cooler. --bin-size is a "
          "mandatory argument and must be different than 0 when processing .mcool files.\n",
          c.path_to_input_matrix.string());
    } else {
      absl::StrAppendFormat(&errors,
                            "Validation of file \"%s\" failed with the following error: %s.\n",
                            c.path_to_input_matrix.string(), e.what());
    }
  }

  if (!c.force && boost::filesystem::exists(c.path_to_output_matrix)) {
    const auto path_type = boost::filesystem::status(c.path_to_output_matrix).type();
    if (path_type == boost::filesystem::directory_file) {
      absl::StrAppendFormat(&errors,
                            "File \"%s\" already exists and is actually a directory. Please remove "
                            "the directory and try again.",
                            c.path_to_output_matrix.string());
    } else {
      absl::StrAppendFormat(&errors, "File \"%s\" already exists. Pass --force to overwrite",
                            c.path_to_output_matrix.string());
    }
  }
  return errors;
}

bool Cli::validate() {  // TODO: Refactor this function
  std::string errors;
  if (this->_cli.get_subcommand("eval")->parsed()) {
    errors = this->validate_eval_subcommand();
  } else if (this->_cli.get_subcommand("filter-barriers")->parsed()) {
    errors = this->validate_filter_barriers_subcommand();
  } else if (this->_cli.get_subcommand("noisify")->parsed()) {
    errors = this->validate_noisify_subcommand();
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

bool Cli::is_ok() const { return (this->_exit_code != 0) && this->_subcommand != subcommand::help; }
Cli::subcommand Cli::get_subcommand() const { return this->_subcommand; }
modle::tools::config Cli::parse_arguments() {
  this->_exec_name = *this->_argv;
  try {
    this->_cli.parse(this->_argc, this->_argv);
    if (this->_cli.get_subcommand("evaluate")->parsed()) {
      this->_subcommand = subcommand::eval;
    } else if (this->_cli.get_subcommand("filter-barriers")->parsed()) {
      this->_subcommand = subcommand::filter_barriers;
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
  } catch (const std::exception& e) {
    this->_exit_code = 1;
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "An unexpected error has occurred while parsing CLI arguments: {}. If you see this "
            "message, please file an issue on GitHub"),
        e.what()));

  } catch (...) {
    this->_exit_code = 1;
    throw std::runtime_error(
        "An unknown error occurred while parsing CLI arguments! If you see this message, please "
        "file an issue on GitHub");
  }

  try {
    this->_exit_code = static_cast<int>(this->validate());
    return this->_config;
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("The following error occourred while validating CLI arguments: {}"), e.what()));
  }
}

int Cli::get_exit_code() const { return this->_exit_code; }

}  // namespace modle::tools
