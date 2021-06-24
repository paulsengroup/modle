#include "./cli.hpp"

#include <absl/strings/match.h>        // for EndsWith
#include <absl/strings/str_cat.h>      // for StrAppend, StrCat
#include <absl/strings/string_view.h>  // for string_view
#include <absl/strings/strip.h>        // for StripSuffix
#include <fmt/format.h>                // for format, FMT_STRING, print
#include <fmt/ostream.h>               // for formatbuf<>::int_type

#include <CLI/App.hpp>         // for Option_group, App
#include <CLI/Config.hpp>      // IWYU pragma: keep for ConfigBase
#include <CLI/Error.hpp>       // for ParseError
#include <CLI/Formatter.hpp>   // IWYU pragma: keep for Formatter
#include <CLI/Option.hpp>      // for Option
#include <CLI/Validators.hpp>  // for PositiveNumber, NonNegativeNumber, Range, Existing...
#include <cmath>               // for round
#include <cstdint>             // for uint32_t, uint64_t
#include <filesystem>          // for path, exists, operator<<, is_empty, is_directory
#include <limits>              // for numeric_limits
#include <sstream>             // for basic_stringbuf<>::int_type, basic_stringbuf<>::po...
#include <stdexcept>           // for invalid_argument, out_of_range
#include <string>              // for allocator, string, basic_string

#include "modle/common/common.hpp"  // for bp_t
#include "modle/common/config.hpp"  // for Config

namespace modle {

void Cli::make_cli() {
  this->_cli.description("in-silico modeling of DNA loop extrusion.");
  auto* io = this->_cli.add_option_group("Input/Output", "");
  auto* gen = this->_cli.add_option_group("Generic", "");
  auto* prob = this->_cli.add_option_group("Probabilities", "");
  auto* rand = this->_cli.add_option_group("Random", "");
  auto* hidden = this->_cli.add_option_group("", "");
  auto* extr_barr = this->_cli.add_option_group("Extrusion Barriers", "");
  auto* extr_barr_mandatory = extr_barr->add_option_group("Mandatory", "")->require_option(1);

  auto remove_trailing_zero_from_floats = [](const std::string& s) {
    return std::string{absl::StripSuffix(s, ".0")};
  };

  // clang-format off
    io->add_option(
            "-c,--chromosome-size-file",
            this->_config.path_to_chrom_sizes,
            "Path to file with chromosome size (should be in a TSV-like format with two fields per line).")
            ->check(CLI::ExistingFile)->required();

    io->add_option(
            "--chromosome-subrange-file",
            this->_config.path_to_chrom_subranges,
            "Path to BED file with subranges of the chromosomes to simulate.")
            ->check(CLI::ExistingFile);

    io->add_option(
            "-o,--output-prefix",
            this->_config.path_to_output_prefix,
            "Output prefix. Can be a full or relative path without file extension.")
            ->required();

    io->add_flag(
            "--force",
            this->_config.force,
            "Overwrite existing output files.")
            ->capture_default_str();

    io->add_flag(
            "--write-contacts-for-excluded-chroms",
            this->_config.write_contacts_for_ko_chroms,
            "Write contacts for all chromosomes, even those where loop extrusion was not simulated. In the latter case ModLE will only write to the chrom, bins and indexes datasets.")
            ->capture_default_str();

    io->add_option(
            "--feature-beds",
            this->_config.path_to_feature_bed_files,
            "Path to one or more BED files containing features used to compute the total number of contacts between pairs of features. Pairs of features with a non-zero number of contacts will be written to a BEDPE file. When a single BED file is specified, the output will only contain within-feature contacts.")
            ->check(CLI::ExistingFile);

    // TODO Add flag to write contacts in scool format when --feature-beds is passed. Cell 0 should store contacts for the entire genome,
    // while additional cells should store contacts for the smaller simulations used to generate contacts for pairs of features.

    gen->add_option(
            "-b,--bin-size",
            this->_config.bin_size,
            "Bin size in base pairs.")
            ->check(CLI::PositiveNumber)
            ->transform(remove_trailing_zero_from_floats)
            ->capture_default_str();

    gen->add_option(
            "--ncells",
            this->_config.num_cells,
            "Number of cells to simulate.")
            ->check(CLI::PositiveNumber)
            ->transform(remove_trailing_zero_from_floats)
            ->capture_default_str();

    gen->add_option(
            "-t,--threads",
            this->_config.nthreads,
            "Max size of the thread pool to use to run the simulation. By default ModLE will attempt to use threads available.")
            ->check(CLI::PositiveNumber)
            ->transform(remove_trailing_zero_from_floats)
            ->capture_default_str();

    gen->add_option(
            "-w,--diagonal-width",
            this->_config.diagonal_width,
            "Width of the matrix that will be used to store contacts, which visually corresponds to the width of the diagonal in a typical square contact matrix.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zero_from_floats))
            ->capture_default_str();

    gen->add_option(
            "--number-of-iterations",
            this->_config.simulation_iterations,
            "Number of simulation iterations.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zero_from_floats))
            ->capture_default_str();
    
    gen->add_option(
            "--target-contact-density",
            this->_config.target_contact_density,
            "Target per-chromosome average contact number after which the simulation is halted.")
            ->check(CLI::PositiveNumber);

    gen->add_option(
            "--lefs-per-mbp",
            this->_config.number_of_lefs_per_mbp,
            "Number of loop extrusion factors (LEFs) per Mbp of DNA to be simulated.")
            ->check(CLI::NonNegativeNumber)
            ->required();

    gen->add_option(
            "--avg-lef-lifetime",
            this->_config.average_lef_lifetime,
            "Average loop extrusion factor lifetime in base pairs, assuming the extruder is free to move along the DNA without being obstructed.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zero_from_floats))
            ->capture_default_str();

    gen->add_option(
            "--hard-stall-multiplier",
            this->_config.hard_stall_multiplier,
            "Coefficient to control the increase in the stability of the LEF-DNA bind when a LEF is stalled in both direction by two extrusion barriers in convergent orientation. Setting this to 1 makes hard stalls strength identical to that of normal stalls.")
            ->check(CLI::PositiveNumber)
            ->capture_default_str();

    gen->add_option(
            "--soft-stall-multiplier",
            this->_config.soft_stall_multiplier,
            "Coefficient to control the increase in the stability of the LEF-DNA bind when a LEF is stalled by an extrusion barrier in the non-blocking orientation. Setting this to 1 makes soft stalls strength identical to that of normal stalls.")
            ->check(CLI::PositiveNumber)
            ->capture_default_str();

    gen->add_option(
            "--fwd-extrusion-speed",
            this->_config.fwd_extrusion_speed,
            "Extrusion speed in forward direction (i.e. distance in bp traveled by a LEF in forward direction at each simulation iteration). Defaults to half of the distance specified through --bin-size.")
            ->transform((remove_trailing_zero_from_floats))
            ->check(CLI::NonNegativeNumber);

    gen->add_option(
            "--rev-extrusion-speed",
            this->_config.rev_extrusion_speed,
            "Extrusion speed in reverse direction (i.e. distance in bp traveled by a LEF in reverse direction at each simulation iteration). Defaults to half of the distance specified through --bin-size.")
            ->transform((remove_trailing_zero_from_floats))
            ->check(CLI::NonNegativeNumber);

    gen->add_option(
            "--fwd-extrusion-speed-std",
            this->_config.fwd_extrusion_speed_std,
            "Standard deviation of the normal distribution used to sample the distance to advance a given LEF in forward direction at every iteration.\nWhen a number other than 0 is passed to this parameter, the distribution has mean equal to the speed specified through --fwd-extrusion-speed.\nSpecifying a std equal to 0 will cause all LEFs to move exactly --fwd-extrusion-speed bp at every iteration.\nWhen the specified std is less than 1, then it will be interpreted as a percentage (i.e. the actual std will be equal to pct * fwd LEF extr. speed.")
            ->check(CLI::NonNegativeNumber)
            ->capture_default_str();

    gen->add_option(
            "--rev-extrusion-speed-std",
            this->_config.rev_extrusion_speed_std,
            "Standard deviation of the normal distribution used to sample the distance to advance a given LEF in reverse direction at every iteration.\nWhen a number other than 0 is passed to this parameter, the distribution has mean equal to the speed specified through --rev-extrusion-speed.\nSpecifying a std equal to 0 will cause all LEFs to move exactly --rev-extrusion-speed bp at every iteration.\nWhen the specified std is less than 1, then it will be interpreted as a percentage (i.e. the actual std will be equal to pct * rev LEF extr. speed.")
            ->check(CLI::NonNegativeNumber)
            ->capture_default_str();

    gen->add_flag(
            "--skip-burn-in",
            this->_config.skip_burnin,
            "Skip burn-in phase and start counting contacts from the first extrusion round.")
            ->capture_default_str();

    gen->add_option(
            "--seed",
            this->_config.seed,
            "Seed to use for random number generation.")
            ->check(CLI::NonNegativeNumber)
            ->transform((remove_trailing_zero_from_floats))
            ->capture_default_str();

    gen->add_option(
            "--deletion-size",
            this->_config.deletion_size,
            "Size of deletion in bp. Used to enable/disable extrusion barriers and compute the total number of contacts between pairs of features. Specify 0 to compute all possible combinations.")
            ->check(CLI::NonNegativeNumber)
            ->transform((remove_trailing_zero_from_floats))
            ->capture_default_str();

    prob->add_option(
            "--probability-of-lef-bypass",
            this->_config.probability_of_extrusion_unit_bypass,
            "Probability that a loop extruding factor (LEF) will not block when meeting another LEF.")
            ->check(CLI::Range(0.0, 1.0))
            ->capture_default_str();

    rand->add_option(
            "--lef-fraction-for-contact-sampling",
            this->_config.lef_fraction_contact_sampling,
            "Fraction of LEFs to use when sampling interactions at every iteration.")
            ->check(CLI::Range(0.0, 1.0))
            ->capture_default_str();

    extr_barr_mandatory->add_option(
            "--extrusion-barrier-file",
            this->_config.path_to_extr_barriers,
            "Path to BED file containing the extrusion barriers to be used in the simulation. Barriers corresponding to chromosomes that are not being simulated will be silently discarded.")
            ->check(CLI::ExistingFile);

    extr_barr->add_option(
            "--probability-of-barrier-block",
            this->_config.probability_of_extrusion_barrier_block,
            "Probability of block of an extrusion unit by an extrusion barrier. Used to indirectly set --ctcf-occupied-probability-of-transition-to-self.")
            ->check(CLI::Range(0.0, 1.0));

    extr_barr->add_option(
            "--ctcf-occupied-probability-of-transition-to-self",
            this->_config.ctcf_occupied_self_prob,
            "Transition probability from CTCF in occupied state to CTCF in occupied state.")
            ->check(CLI::Range(0.0, 1.0));

    extr_barr->add_option(
            "--ctcf-not-occupied-probability-of-transition-to-self",
            this->_config.ctcf_not_occupied_self_prob,
            "Transition probability from CTCF in not-occupied state to CTCF in not-occupied state.")
            ->check(CLI::Range(0.0, 1.0))
            ->capture_default_str();

    extr_barr->add_flag(
            "--exclude-chrom-wo-barriers,!--keep-chrom-without-barriers",
            this->_config.exclude_chrom_wo_extr_barriers,
            fmt::format(FMT_STRING("Do not simulate loop extrusion on chromosomes without any extrusion barrier. Default: {}"), this->_config.exclude_chrom_wo_extr_barriers ? "--exclude-chrom-wo-barriers" : "--keep-chrom-without-barriers"))
            ->capture_default_str();

    hidden->add_flag("--skip-output", this->_config.skip_output, "Don't write output files and plots. Useful for profiling");

    gen->get_option("--target-contact-density")->excludes(gen->get_option("--number-of-iterations"));
    extr_barr->get_option("--probability-of-barrier-block")->excludes("--ctcf-occupied-probability-of-transition-to-self");
    gen->get_option("--deletion-size")->needs(io->get_option("--feature-beds"));
  // clang-format on
}

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { this->make_cli(); }
const Config& Cli::parse_arguments() {
  this->_cli.name(this->_exec_name);
  this->_cli.parse(this->_argc, this->_argv);
  this->validate_and_transform_args();
  this->_config.argc = _argc;
  this->_config.argv = _argv;
  return this->_config;
}

std::string Cli::process_paths_and_check_for_collisions(modle::Config& c) {
  std::string collisions;

  c.path_to_output_file_cool = c.path_to_output_prefix;
  c.path_to_output_file_bedpe =
      !c.path_to_feature_bed_files.empty() ? c.path_to_output_prefix : std::filesystem::path{};
  c.path_to_log_file = c.path_to_output_prefix;

  c.path_to_output_file_cool += ".cool";
  c.path_to_output_file_bedpe += !c.path_to_output_file_bedpe.empty() ? ".bedpe" : "";
  c.path_to_log_file += ".log";

  if (c.force || c.skip_output) {
    return "";
  }

  auto check_for_path_collisions = [](const std::filesystem::path& path) -> std::string {
    if (std::filesystem::exists(path)) {
      if (std::filesystem::is_directory(path)) {
        return fmt::format(
            FMT_STRING("Refusing to run the simulation because output file {} already "
                       "exist (and is actually a {}directory). {}.\n"),
            path, std::filesystem::is_empty(path) ? "" : "non-empty ",
            std::filesystem::is_empty(path)
                ? " Pass --force to overwrite"
                : "You should specify a different output path, or manually remove the "
                  "existing directory");
      }
      return fmt::format(
          FMT_STRING("Refusing to run the simulation because output file {} already exist. Pass "
                     "--force to overwrite.\n"),
          path);
    }
    if (std::filesystem::is_directory(path) && !std::filesystem::is_empty(path)) {
      return fmt::format(
          FMT_STRING("Refusing to run the simulation because output file {} is a "
                     "non-empty directory. You should specify a different output path, or "
                     "manually remove the existing directory.\n"),
          path);
    }
    return {};
  };

  if (std::filesystem::exists(c.path_to_log_file)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_log_file));
  }
  if (std::filesystem::exists(c.path_to_output_file_cool)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_output_file_cool));
  }
  if (!c.path_to_output_file_bedpe.empty() &&
      std::filesystem::exists(c.path_to_output_file_bedpe)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_output_file_bedpe));
  }

  return collisions;
}

int Cli::exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }

void Cli::validate_and_transform_args() {
  if (auto& speed = this->_config.rev_extrusion_speed; speed == std::numeric_limits<bp_t>::max()) {
    speed =
        static_cast<bp_t>(std::round(static_cast<double>(this->_config.bin_size) / 2.0));  // NOLINT
  }
  if (auto& speed = this->_config.fwd_extrusion_speed; speed == std::numeric_limits<bp_t>::max()) {
    speed =
        static_cast<bp_t>(std::round(static_cast<double>(this->_config.bin_size) / 2.0));  // NOLINT
  }

  if (auto& stddev = this->_config.fwd_extrusion_speed_std; stddev > 0 && stddev < 1) {
    stddev *= static_cast<double>(this->_config.fwd_extrusion_speed);
  }
  if (auto& stddev = this->_config.rev_extrusion_speed_std; stddev > 0 && stddev < 1) {
    stddev *= static_cast<double>(this->_config.rev_extrusion_speed);
  }

  {
    const auto pno = 1.0 - this->_config.ctcf_not_occupied_self_prob;
    if (this->_config.ctcf_occupied_self_prob == 0) {
      const auto occ = this->_config.probability_of_extrusion_barrier_block;
      const auto pon = (pno - (occ * pno)) / occ;
      this->_config.ctcf_occupied_self_prob = 1.0 - pon;
    } else {
      const auto pon = 1.0 - this->_config.ctcf_occupied_self_prob;
      const auto occ = pno / (pno + pon);
      this->_config.probability_of_extrusion_barrier_block = occ;
    }
  }
}

}  // namespace modle
