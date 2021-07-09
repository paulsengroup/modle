#include "./cli.hpp"

#include <absl/strings/match.h>        // for EndsWith
#include <absl/strings/str_cat.h>      // for StrAppend, StrCat
#include <absl/strings/string_view.h>  // for string_view
#include <absl/strings/strip.h>        // for StripSuffix
#include <fmt/format.h>                // for format, FMT_STRING, print
#include <fmt/ostream.h>               // for formatbuf<>::int_type

#include <CLI/App.hpp>                // for Option_group, App
#include <CLI/Config.hpp>             // IWYU pragma: keep for ConfigBase
#include <CLI/Error.hpp>              // for ParseError
#include <CLI/Formatter.hpp>          // IWYU pragma: keep for Formatter
#include <CLI/Option.hpp>             // for Option
#include <CLI/Validators.hpp>         // for PositiveNumber, NonNegativeNumber, Range, Existing...
#include <boost/filesystem/path.hpp>  // for path, exists, operator<<, is_empty, is_directory
#include <cmath>                      // for round
#include <cstdint>                    // for uint32_t, uint64_t
#include <limits>                     // for numeric_limits
#include <sstream>                    // for basic_stringbuf<>::int_type, basic_stringbuf<>::po...
#include <stdexcept>                  // for invalid_argument, out_of_range
#include <string>                     // for allocator, string, basic_string

#include "modle/common/common.hpp"  // for bp_t
#include "modle/common/config.hpp"  // for Config

namespace modle {

std::string remove_trailing_zero_from_floats(const std::string& s) {
  return std::string{absl::StripSuffix(s, ".0")};
}

void Cli::make_simulation_subcommand() {
  auto& s = *this->_cli.add_subcommand(
      "simulate",
      "Perform a single genome-wide simulation and output the resulting contacts to a .cool file.");
  auto& io = *s.add_option_group("Input/Output", "");
  auto& gen = *s.add_option_group("Generic", "");
  auto& prob = *s.add_option_group("Probabilities", "");
  auto& rand = *s.add_option_group("Random", "");
  auto& hidden = *s.add_option_group("", "");
  auto& extr_barr = *s.add_option_group("Extrusion Barriers", "");
  auto& extr_barr_mandatory = *extr_barr.add_option_group("Mandatory", "")->require_option(1);

  // clang-format off
    io.add_option(
            "-c,--chrom-sizes",
            this->_config.path_to_chrom_sizes,
            "Path to file with chromosome sizes in chrom.sizes format.")
            ->check(CLI::ExistingFile)->required();

    io.add_option(
            "--chrom-subranges",
            this->_config.path_to_chrom_subranges,
            "Path to BED file with subranges of the chromosomes to simulate.")
            ->check(CLI::ExistingFile);

    io.add_option(
            "-o,--output-prefix",
            this->_config.path_to_output_prefix,
            "Output prefix. Can be a full or relative path including the file name but without file extension.\n"
            "Example: -o /tmp/mymatrix will cause MoDLE to write contacts to a file named \"/tmp/mymatrix.cool\", "
            "while a log file will be saved at \"/tmp/mymatrix.log\".")
            ->required();

    io.add_flag(
            "-f,--force",
            this->_config.force,
            "Force overwrite of output files if they exists.")
            ->capture_default_str();

    io.add_flag(
            "--write-contacts-for-excluded-chroms",
            this->_config.write_contacts_for_ko_chroms,
            "Write contacts for all chromosomes, even those where loop extrusion was not simulated.\n"
            "In the latter case MoDLE will only write to the chrom, bins and indexes datasets to the contact matrix.")
            ->capture_default_str();

    gen.add_option(
            "-b,--bin-size",
            this->_config.bin_size,
            "Bin size in base pairs.")
            ->check(CLI::PositiveNumber)
            ->transform(remove_trailing_zero_from_floats)
            ->capture_default_str();

    gen.add_option(
            "--ncells",
            this->_config.num_cells,
            "Number of cells to simulate.\n"
            "Loop extrusion will be simulated independently on each chromosome across --ncells simulation instances.\n"
            "The total number of contacts for a given chromosome is obtained by aggregating contacts from all the "
            "simulation instances for that chromosome.")
            ->check(CLI::PositiveNumber)
            ->transform(remove_trailing_zero_from_floats)
            ->capture_default_str();

    gen.add_option(
            "-t,--threads",
            this->_config.nthreads,
            "Number of worker threads used to run the simulation.\n"
            "By default MoDLE will try to use all available threads.")
            ->check(CLI::PositiveNumber)
            ->transform(remove_trailing_zero_from_floats)
            ->capture_default_str();

    gen.add_option(
            "-w,--diagonal-width",
            this->_config.diagonal_width,
            "Diagonal width of the in-memory contact matrix that will be used to store contacts during the simulation.\n"
            "This setting affects the maximum distance of a pair of bins whose interactions will be tracked by MoDLE.\n"
            "Setting --diagonal-width to a very large value (i.e. more than few Mbp) will dramatically inflate MoDLE's memory footprint.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zero_from_floats))
            ->capture_default_str();

    gen.add_option(
            "--number-of-iterations",
            this->_config.simulation_iterations,
            "Number of simulation iterations to run on each cell.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zero_from_floats))
            ->capture_default_str();

    gen.add_option(
            "--target-contact-density",
            this->_config.target_contact_density,
            "The average number of contacts to generate for each chromosome before the simulation is halted.")
            ->check(CLI::PositiveNumber);

    gen.add_option(
            "--lefs-per-mbp",
            this->_config.number_of_lefs_per_mbp,
            "Number of loop extrusion factors (LEFs) per Mbp of simulated DNA.")
            ->check(CLI::NonNegativeNumber)
            ->required();

    gen.add_option(
            "--avg-lef-lifetime",
            this->_config.average_lef_lifetime,
            "Average LEF lifetime.\n"
            "This is equal to the average distance an unencumbered LEF is able to move along the DNA before it is released.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zero_from_floats))
            ->capture_default_str();

    gen.add_option(
            "--hard-stall-multiplier", // TODO: Should we rename this?
            this->_config.hard_stall_multiplier,
            "Coefficient to control the DNA-binding stability of LEFs that are stalled on both sides by a pair of "
            "extrusion barriers in convergent orientation.\n"
            "Setting this to 1 makes the DNA-binding stability of hard-stalled LEFs identical to that of unobstructed LEFs.")
            ->check(CLI::NonNegativeNumber)
            ->capture_default_str();

    gen.add_option(
            "--soft-stall-multiplier", // TODO: Should we rename this?
            this->_config.soft_stall_multiplier,
            "Coefficient to control the DNA-binding stability of a LEF that is stalled by one or more extrusion barriers "
            "in a non-blocking orientation.\n"
            "Setting this to 1 makes the DNA-binding stability of soft-stalled LEFs identical to that of unobstructed LEFs.")
            ->check(CLI::NonNegativeNumber)
            ->capture_default_str();

    gen.add_option(
            "--fwd-extrusion-speed",
            this->_config.fwd_extrusion_speed,
            "Average distance in bp covered by a LEF in forward direction during one iteration.\n"
            "By deafult this is set to half of the bin size specified through --bin-size.")
            ->transform((remove_trailing_zero_from_floats))
            ->check(CLI::NonNegativeNumber);

    gen.add_option(
            "--rev-extrusion-speed",
            this->_config.rev_extrusion_speed,
            "Same as --fwd-extrusion-speed but for the reverse direction.")
            ->transform((remove_trailing_zero_from_floats))
            ->check(CLI::NonNegativeNumber);

    gen.add_option(
            "--fwd-extrusion-speed-std",
            this->_config.fwd_extrusion_speed_std,
            "Standard deviation of the normal distribution from which LEF forward moves are sampled.\n"
            "When a number other than 0 is passed to this parameter, the move distribution has mean "
            "equal to the speed specified through --fwd-extrusion-speed.\n"
            "Specifying a --fwd-extrusion-speed-std=0 will cause all LEFs to move exactly --fwd-extrusion-speed bp "
            "at every iteration.\n"
            "When the specified --fwd-extrusion-speed-std is less than 1, then it will be interpreted as a percentage "
            "of the forward extrusion speed (i.e. the standard deviation of the distribution will be equal to "
            "--fwd-extrusion-speed-std * --fwd-extrusion-speed).")
            ->check(CLI::NonNegativeNumber)
            ->capture_default_str();

    gen.add_option(
            "--rev-extrusion-speed-std",
            this->_config.rev_extrusion_speed_std,
            "Same as --fwd-extrusion-speed-std but for the reverse direction.")
            ->check(CLI::NonNegativeNumber)
            ->capture_default_str();

    gen.add_flag(
            "--skip-burn-in",
            this->_config.skip_burnin,
            "Skip the burn-in phase and start counting contacts from the first extrusion round.")
            ->capture_default_str();

    gen.add_option(
            "--seed",
            this->_config.seed,
            "Base seed to use for random number generation.")
            ->check(CLI::NonNegativeNumber)
            ->transform((remove_trailing_zero_from_floats))
            ->capture_default_str();

    prob.add_option(
            "--probability-of-lef-bypass",
            this->_config.probability_of_extrusion_unit_bypass,
            "Probability that two LEFs will bypass each other when meeting.")
            ->check(CLI::Range(0.0, 1.0))
            ->capture_default_str();

    rand.add_option(
            "--lef-fraction-for-contact-sampling",
            this->_config.lef_fraction_contact_sampling,
            "Fraction of LEFs to use when sampling interactions at every iteration.")
            ->check(CLI::Range(0.0, 1.0))
            ->capture_default_str();

    rand.add_flag(
            "--randomize-contacts",
            this->_config.randomize_contacts,
            "Enable randomization when collecting contacts.\n"
            "Contacts are randomized by drawing offsets from a genextreme distrubution.\n"
            "The distribution properties can be controlled through the parameters --mu, --sigma and --xi.")
            ->capture_default_str();

    rand.add_option(
            "--mu,--genextr-location",
            this->_config.genextreme_mu,
            "Location parameter (mu) of the generalized extreme value distribution used to add noise to the contact matrix.")
            ->check(CLI::NonNegativeNumber)
            ->capture_default_str();

    rand.add_option(
            "--sigma,--genextr-scale",
            this->_config.genextreme_sigma,
            "Scale parameter (sigma) of the generalized extreme value distribution used to add noise to the contact matrix.")
            ->check(CLI::NonNegativeNumber)
            ->capture_default_str();

    rand.add_option(
            "--xi,--genextr-shape",
            this->_config.genextreme_xi,
            "Shape parameter (xi) of the generalized extreme value distribution used to add noise to the contact matrix.")
            ->check(CLI::NonNegativeNumber)
            ->capture_default_str();

    extr_barr_mandatory.add_option(
            "--extrusion-barrier-file",
            this->_config.path_to_extr_barriers,
            "Path to BED file containing the extrusion barriers to be used in the simulation.\n"
            "Barriers corresponding to chromosomes that are not part of the simulation will be ignored.")
            ->check(CLI::ExistingFile);

    extr_barr.add_option(
            "--probability-of-barrier-block",
            this->_config.probability_of_extrusion_barrier_block,
            "Probability that a LEF will be stalled upon reaching an extrusion barrier.\n"
            "This parameter provides an easier mean to set --ctcf-occupied-probability-of-transition-to-self.")
            ->check(CLI::Range(0.0, 1.0));

    extr_barr.add_option(
            "--ctcf-occupied-probability-of-transition-to-self",
            this->_config.ctcf_occupied_self_prob,
            "Probability that an extrusion barrier that was occupied during the previous iteration will also be "
            "occupied in the current iteration.")
            ->check(CLI::Range(0.0, 1.0));

    extr_barr.add_option(
            "--ctcf-not-occupied-probability-of-transition-to-self",
            this->_config.ctcf_not_occupied_self_prob,
            "Probability that an extrusion barrier that was not occupied during the previous iteration will remain "
            "not occupied in the current iteration.")
            ->check(CLI::Range(0.0, 1.0))
            ->capture_default_str();

    extr_barr.add_flag(
            "--exclude-chrom-wo-barriers,!--keep-chrom-without-barriers",
            this->_config.exclude_chrom_wo_extr_barriers,
            "Do not simulate loop extrusion on chromosomes without any extrusion barrier.")
            ->capture_default_str();

    hidden.add_flag("--skip-output", this->_config.skip_output, "Don't write output files. Useful for profiling");

    gen.get_option("--target-contact-density")->excludes(gen.get_option("--number-of-iterations"));
    extr_barr.get_option("--probability-of-barrier-block")->excludes("--ctcf-occupied-probability-of-transition-to-self");
    rand.get_option("--mu")->needs(rand.get_option("--randomize-contacts"));
    rand.get_option("--sigma")->needs(rand.get_option("--randomize-contacts"));
    rand.get_option("--xi")->needs(rand.get_option("--randomize-contacts"));
  // clang-format on
}

void Cli::make_cli() {
  this->_cli.description("Stochastic modeling of DNA loop extrusion.");

  this->make_simulation_subcommand();
  // this->make_perturbate_subcommand();
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
      !c.path_to_feature_bed_files.empty() ? c.path_to_output_prefix : boost::filesystem::path{};
  c.path_to_log_file = c.path_to_output_prefix;

  c.path_to_output_file_cool += ".cool";
  c.path_to_output_file_bedpe += !c.path_to_output_file_bedpe.empty() ? ".bedpe.gz" : "";
  c.path_to_log_file += ".log";

  if (c.force || c.skip_output) {
    return "";
  }

  auto check_for_path_collisions = [](const boost::filesystem::path& path) -> std::string {
    if (boost::filesystem::exists(path)) {
      if (boost::filesystem::is_directory(path)) {
        return fmt::format(
            FMT_STRING("Refusing to run the simulation because output file {} already "
                       "exist (and is actually a {}directory). {}.\n"),
            path, boost::filesystem::is_empty(path) ? "" : "non-empty ",
            boost::filesystem::is_empty(path)
                ? " Pass --force to overwrite"
                : "You should specify a different output path, or manually remove the "
                  "existing directory");
      }
      return fmt::format(
          FMT_STRING("Refusing to run the simulation because output file {} already exist. Pass "
                     "--force to overwrite.\n"),
          path);
    }
    if (boost::filesystem::is_directory(path) && !boost::filesystem::is_empty(path)) {
      return fmt::format(
          FMT_STRING("Refusing to run the simulation because output file {} is a "
                     "non-empty directory. You should specify a different output path, or "
                     "manually remove the existing directory.\n"),
          path);
    }
    return {};
  };

  if (boost::filesystem::exists(c.path_to_log_file)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_log_file));
  }
  if (boost::filesystem::exists(c.path_to_output_file_cool)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_output_file_cool));
  }
  if (!c.path_to_output_file_bedpe.empty() &&
      boost::filesystem::exists(c.path_to_output_file_bedpe)) {
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
