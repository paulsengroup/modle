#pragma once

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
#include <cstdint>             // for uint32_t, uint64_t
#include <cstdio>              // for stderr
#include <filesystem>          // for path, exists, operator<<, is_empty, is_directory
#include <sstream>             // for basic_stringbuf<>::int_type, basic_stringbuf<>::po...
#include <stdexcept>           // for invalid_argument, out_of_range
#include <string>              // for allocator, string, basic_string

#include "modle/common.hpp"  // for Bp
#include "modle/config.hpp"  // for Config

namespace modle {

class Cli {
 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  Config _config;
  CLI::App _cli{};

  inline void MakeCli() {
    this->_cli.description("in-silico modeling of DNA loop extrusion.");
    auto* io = this->_cli.add_option_group("Input/Output", "");
    auto* gen = this->_cli.add_option_group("Generic", "");
    auto* prob = this->_cli.add_option_group("Probabilities", "");
    auto* rand = this->_cli.add_option_group("Random", "");
    auto* hidden = this->_cli.add_option_group("", "");
    auto* extr_barr = this->_cli.add_option_group("Extrusion Barriers", "");
    auto* extr_barr_mandatory = extr_barr->add_option_group("Mandatory", "")->require_option(1);

    auto remove_trailing_zeros_from_floats = [](const std::string& s) {
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
            "-o,--output-file",
            this->_config.path_to_output_file,
            "Path to output file (cooler format).")
            ->transform([](const std::string& s) {
              if (absl::EndsWith(s, ".cool")) {
                return s;
              }
              return absl::StrCat(s, ".cool");})
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

    io->add_flag(
            "--write-raw-contacts",
            this->_config.write_raw_contacts,
            "Output raw contact frequencies to a .cool file.")
            ->capture_default_str();

    io->add_flag(
            "--write-contacts-w-noise",
            this->_config.write_contacts_with_noise,
            "Output contact frequencies to a .cool file after adding random noise. When this flag is true, contact frequencies with noise are stored in a .cool file suffixed with \"_w_noise\"")
            ->capture_default_str();

    gen->add_option(
            "-b,--bin-size",
            this->_config.bin_size,
            "Bin size in base pairs.")
            ->check(CLI::PositiveNumber)
            ->transform(remove_trailing_zeros_from_floats)
            ->capture_default_str();

    gen->add_option(
            "--ncells",
            this->_config.ncells,
            "Number of cells to simulate.")
            ->check(CLI::PositiveNumber)
            ->transform(remove_trailing_zeros_from_floats)
            ->capture_default_str();

    gen->add_option(
            "-t,--threads",
            this->_config.nthreads,
            "Max size of the thread pool to use to run the simulation. By default ModLE will attempt to use threads available.")
            ->check(CLI::PositiveNumber)
            ->transform(remove_trailing_zeros_from_floats)
            ->capture_default_str();

    gen->add_option(
            "-w,--diagonal-width",
            this->_config.diagonal_width,
            "Width of the matrix that will be used to store contacts, which visually corresponds to the width of the diagonal in a typical square contact matrix.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zeros_from_floats))
            ->capture_default_str();

    gen->add_option(
            "--number-of-iterations",
            this->_config.simulation_iterations,
            "Number of simulation iterations.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zeros_from_floats))
            ->capture_default_str();
    
    gen->add_option(
            "--target-contact-density",
            this->_config.target_contact_density,
            "Target per-chromosome average contact number after which the simulation is halted.")
            ->check(CLI::PositiveNumber)
            ->excludes(gen->get_option("--number-of-iterations"));

    gen->add_option(
            "--lefs-per-mbp",
            this->_config.number_of_lefs_per_mbp,
            "Number of loop extrusion factors (LEFs) per Mbp to be simulated.")
            ->check(CLI::NonNegativeNumber)
            ->required();

    gen->add_option(
            "--avg-lef-lifetime",
            this->_config.average_lef_lifetime,
            "Average loop extrusion factor lifetime in base pairs, assuming the extruder is free to move along the DNA without being obstructed.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zeros_from_floats))
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
            ->check(CLI::NonNegativeNumber);

    gen->add_option(
            "--rev-extrusion-speed",
            this->_config.rev_extrusion_speed,
            "Extrusion speed in reverse direction (i.e. distance in bp traveled by a LEF in reverse direction at each simulation iteration). Defaults to half of the distance specified through --bin-size.")
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
            "--allow-lef-lifetime-extension,!--disable-lef-lifetime-extension",
            this->_config.allow_lef_lifetime_extension,
            "Toggle on and off the ability to extend LEF lifetimes in case of hard stalls.")
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
            ->transform((remove_trailing_zeros_from_floats))
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

    rand->add_option(
            "--contact-noise-mean",
            this->_config.random_noise_mean,
            "Mean parameter (in base pairs) of the normal distribution that will be used to introduce noise in the contact matrix when --write-contacts-w-noise is specified.")
            ->check(CLI::PositiveNumber)
            ->capture_default_str();

    rand->add_option(
            "--contact-noise-std",
            this->_config.random_noise_std,
            "Standard deviation parameter (in base pairs) of the normal distribution that will be used to introduce noise in the contact matrix when --write-contacts-w-noise is specified.")
            ->check(CLI::PositiveNumber)
            ->capture_default_str();

    extr_barr_mandatory->add_option(
            "--extrusion-barrier-file",
            this->_config.path_to_extr_barriers,
            "Path to BED file containing the extrusion barriers to be used in the simulation. Barriers corresponding to chromosomes that are not being simulated will be silently discarded.")
            ->check(CLI::ExistingFile);

    extr_barr->add_option(
            "--ctcf-occupied-probability-of-transition-to-self",
            this->_config.ctcf_occupied_self_prob,
            "Transition probability from CTCF in occupied state to CTCF in occupied state.")
            ->check(CLI::Range(0.0, 1.0))
            ->capture_default_str();

    extr_barr->add_option(
            "--ctcf-not-occupied-probability-of-transition-to-self",
            this->_config.ctcf_not_occupied_self_prob,
            "Transition probability from CTCF in not-occupied state to CTCF in not-occupied state.")
            ->check(CLI::Range(0.0, 1.0))
            ->capture_default_str();

    extr_barr->add_flag(
            "--exclude-chr-wo-barriers,!--keep-chr-without-barriers",
            this->_config.exclude_chr_wo_extr_barriers,
            fmt::format(FMT_STRING("Do not simulate loop extrusion on chromosomes without any extrusion barrier. Default: {}"), this->_config.exclude_chr_wo_extr_barriers ? "--exclude-chr-wo-barriers" : "--keep-chr-without-barriers"))
            ->capture_default_str();

    hidden->add_flag("--skip-output", this->_config.skip_output, "Don't write output files and plots. Useful for profiling");
    // clang-format on
  }

 public:
  inline Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) {
    this->MakeCli();
  }
  [[nodiscard]] const Config& parse_arguments() {
    this->_cli.name(this->_exec_name);
    this->_cli.parse(this->_argc, this->_argv);
    this->validate_and_transform_args();
    this->_config.argc = _argc;
    this->_config.argv = _argv;
    return this->_config;
  }

  [[nodiscard]] static inline std::string process_paths_and_check_for_collisions(modle::Config& c) {
    auto base_out = c.path_to_output_file;
    std::string collisions;
    base_out.replace_extension();
    c.path_to_output_file_w_noise = base_out;
    c.path_to_log_file = base_out;
    c.path_to_output_file_w_noise +=
        absl::StrCat("_w_noise", c.path_to_output_file.extension().string());
    c.path_to_log_file += ".log";

    fmt::print(stderr, "path_to_file_w_noise={}\n", c.path_to_output_file_w_noise);

    if (c.force) {
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

    if (c.write_raw_contacts) {
      if (std::filesystem::exists(c.path_to_output_file)) {
        absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_output_file));
      }
    }

    if (c.write_contacts_with_noise) {
      if (std::filesystem::exists(c.path_to_output_file_w_noise)) {
        absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_output_file_w_noise));
      }
    }
    return collisions;
  }

  [[nodiscard]] inline int exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }

 private:
  inline void validate_and_transform_args() {
    if (auto& speed = this->_config.rev_extrusion_speed;
        speed == std::numeric_limits<bp_t>::max()) {
      speed = static_cast<bp_t>(std::round(static_cast<double>(this->_config.bin_size) / 2.0));
    }
    if (auto& speed = this->_config.fwd_extrusion_speed;
        speed == std::numeric_limits<bp_t>::max()) {
      speed = static_cast<bp_t>(std::round(static_cast<double>(this->_config.bin_size) / 2.0));
    }

    if (auto& stddev = this->_config.fwd_extrusion_speed_std; stddev > 0 && stddev < 1) {
      stddev *= static_cast<double>(this->_config.fwd_extrusion_speed);
    }
    if (auto& stddev = this->_config.rev_extrusion_speed_std; stddev > 0 && stddev < 1) {
      stddev *= static_cast<double>(this->_config.rev_extrusion_speed);
    }
  }
};

}  // namespace modle
