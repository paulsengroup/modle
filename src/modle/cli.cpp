// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./cli.hpp"

#include <absl/strings/str_cat.h>      // for StrAppend
#include <absl/strings/str_split.h>    // for SplitIterator, Splitter, StrSplit, operator!=
#include <absl/strings/string_view.h>  // for string_view, basic_string_view, operator""sv
#include <absl/types/span.h>           // for MakeSpan
#include <fmt/format.h>                // for format, FMT_STRING, join, print, make_format_...
#include <fmt/os.h>                    // for output_file, ostream
#include <fmt/ostream.h>               // for formatbuf<>::int_type
#include <toml++/toml.h>               // for array::operator[], operator<<, parse, print_to...

#include <CLI/CLI.hpp>                      // for Option_group, App
#include <algorithm>                        // for max
#include <array>                            // for array
#include <boost/filesystem/operations.hpp>  // for exists, is_empty, is_directory
#include <boost/filesystem/path.hpp>        // for path, operator<<
#include <cassert>                          // for assert
#include <cmath>                            // for round
#include <cstdio>                           // for stderr
#include <exception>                        // for exception
#include <limits>                           // for numeric_limits
#include <sstream>                          // for streamsize, stringstream, basic_ostream
#include <stdexcept>                        // for invalid_argument, out_of_range, runtime_error
#include <string>                           // for allocator, string, basic_string
#include <thread>                           // for hardware_concurrency
#include <vector>                           // for vector

#include "modle/common/common.hpp"  // for bp_t, i64, modle_version_long
#include "modle/common/config.hpp"  // for Config
#include "modle/common/utils.hpp"   // for str_float_to_str_int, parse_numeric_or_throw
#include "modle/cooler/cooler.hpp"  // for Cooler

namespace modle {

static std::string is_odd_number(std::string_view s) {
  try {
    const auto n = utils::parse_numeric_or_throw<i64>(s);
    if (n % 2 == 0) {
      return fmt::format(FMT_STRING("{} is not an odd number"), n);
    }
  } catch (const std::exception& e) {
    return fmt::format(FMT_STRING("Failed parsing number: ({}): {}"), s, e.what());
  }
  return "";
}

static void add_common_options(CLI::App& subcommand, modle::Config& c) {
  auto& s = subcommand;

  auto& io = *s.add_option_group("Input/Output", "");
  auto& gen = *s.add_option_group("Generic", "");
  auto& prob = *s.add_option_group("Probabilities", "");
  auto& rand = *s.add_option_group("Random", "");
  auto& extr_barr = *s.add_option_group("Extrusion Barriers", "");
  auto& dbg = *s.add_option_group("Debug/Profiling", "");

  // clang-format off
  io.add_option(
      "-c,--chrom-sizes",
      c.path_to_chrom_sizes,
      "Path to file with chromosome sizes in chrom.sizes format.")
      ->check(CLI::ExistingFile)->required();

  io.add_option(
      "--chrom-subranges",
      c.path_to_chrom_subranges,
      "Path to BED file with subranges of the chromosomes to simulate.")
      ->check(CLI::ExistingFile);

  io.add_flag(
      "-f,--force",
      c.force,
      "Force overwrite of output files if they exists.")
      ->capture_default_str();

  io.add_flag(
      "-q,--quiet",
      c.quiet,
      "Only log fatal errors.")
      ->capture_default_str();

  gen.add_option(
      "-b,--bin-size",
      c.bin_size,
      "Bin size in base pairs.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  gen.add_option(
      "-t,--threads",
      c.nthreads,
      "Number of simulate_worker threads used to run the simulation.\n"
      "By default MoDLE will try to use all available threads.")
      ->check(CLI::PositiveNumber)
      ->transform(CLI::Bound(1U, std::thread::hardware_concurrency()))
      ->capture_default_str();

  gen.add_option(
      "-w,--diagonal-width",
      c.diagonal_width,
      "Diagonal width of the in-memory contact matrix that will be used to store contacts during the simulation.\n"
      "This setting affects the maximum distance of a pair of bins whose interactions will be tracked by MoDLE.\n"
      "Setting --diagonal-width to a very large value (i.e. more than few Mbp) will dramatically inflate MoDLE's memory footprint.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  gen.add_option(
      "--number-of-iterations",
      c.simulation_epochs,
      "Number of simulation iterations to run on each cell.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  gen.add_option(
      "--target-contact-density",
      c.target_contact_density,
      "The average number of contacts to generate for each chromosome before the simulation is halted.")
      ->check(CLI::PositiveNumber);

  gen.add_option(
      "--lefs-per-mbp",
      c.number_of_lefs_per_mbp,
      "Number of loop extrusion factors (LEFs) per Mbp of simulated DNA.")
      ->check(CLI::NonNegativeNumber)
      ->required();

  gen.add_option(
      "--hard-stall-lef-stability-multiplier",
      c.hard_stall_lef_stability_multiplier,
      "Coefficient to control the DNA-binding stability of LEFs that are stalled on both sides by a pair of "
      "extrusion barriers in convergent orientation.\n"
      "Setting this to 1 makes the DNA-binding stability of hard-stalled LEFs identical to that of unobstructed LEFs.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  gen.add_option(
      "--soft-stall-lef-stability-multiplier",
      c.soft_stall_lef_stability_multiplier,
      "Coefficient to control the DNA-binding stability of a LEF that is stalled by one or more extrusion barriers "
      "in a non-blocking orientation.\n"
      "Setting this to 1 makes the DNA-binding stability of soft-stalled LEFs identical to that of unobstructed LEFs.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  gen.add_option(
      "--fwd-extrusion-speed",
      c.fwd_extrusion_speed,
      "Average distance in bp covered by a LEF in forward direction during one iteration.\n"
      "By deafult this is set to half of the bin size specified through --bin-size.")
      ->transform(utils::str_float_to_str_int)
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  gen.add_option(
      "--rev-extrusion-speed",
      c.rev_extrusion_speed,
      "Same as --fwd-extrusion-speed but for the reverse direction.")
      ->transform(utils::str_float_to_str_int)
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  gen.add_option(
      "--fwd-extrusion-speed-std",
      c.fwd_extrusion_speed_std,
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
      c.rev_extrusion_speed_std,
      "Same as --fwd-extrusion-speed-std but for the reverse direction.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  gen.add_flag(
      "--skip-burnin",
      c.skip_burnin,
      "Skip the burn-in phase and start counting contacts from the first extrusion round.")
      ->capture_default_str();

  gen.add_option(
      "--burnin-target-epochs-for-lef-activation",
       c.burnin_target_epochs_for_lef_activation,
      "Number of epochs over which LEFs are progressively activated and bound to DNA.\n"
      "Note: this number is approximate, as LEFs are activate using a Possion process.\n"
      "By default this parameter is computed based on the average LEF lifetime, overall extrusion speed "
      "and burn-in extrusion speed coefficient (controlled by --avg-lef-lifetime, --fwd/rev-extrusion-speed "
      "and --burn-in-extr-speed-coefficient respectively).")
       ->check(CLI::PositiveNumber)
       ->capture_default_str();

  gen.add_option(
      "--burnin-history-length",
      c.burnin_history_length,
      "Number of epochs used to determine whether a simulation instance has reached a stable state.\n"
      "This is used to decide whether to terminate the burn-in phase at a given epoch.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

   gen.add_option(
      "--burnin-smoothing-window-size",
      c.burnin_smoothing_window_size,
      "Window size used to smooth values during the burnin phase.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  gen.add_option(
      "--max-burnin-epochs",
       c.max_burnin_epochs,
       "Maximum number of epochs to spend in burn-in phase.")
       ->check(CLI::PositiveNumber)
       ->capture_default_str();

  gen.add_option(
      "--burnin-extr-speed-coefficient",
      c.burnin_speed_coefficient,
      "Extrusion speed coefficient to apply during the burn-in phase.\n"
      "Setting this to numbers > 1.0 will speed-up the burn-in phase.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  gen.add_option(
      "--seed",
      c.seed,
      "Base seed to use for random number generation.")
      ->check(CLI::NonNegativeNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  gen.add_option(
      "--prob-of-lef-release",
      c.prob_of_lef_release,
      "Probability of LEF release.\n"
      "The actual probability that is used for a specific LEF is affected --hard-stall-multiplier and --hard-stall-multiplier.")
      ->check(CLI::PositiveNumber & CLI::Range(0.0, 1.0))
      ->capture_default_str();

  prob.add_option(
      "--probability-of-lef-bypass",
      c.probability_of_extrusion_unit_bypass,
      "Probability that two LEFs will bypass each other when meeting.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  rand.add_option(
      "--lef-fraction-for-contact-sampling",
      c.lef_fraction_contact_sampling,
      "Fraction of LEFs to use when sampling interactions at every iteration.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  extr_barr.add_option(
      "--extrusion-barrier-file",
      c.path_to_extr_barriers,
      "Path to BED file containing the extrusion barriers to be used in the simulation.\n"
      "Barriers corresponding to chromosomes that are not part of the simulation will be ignored.")
      ->check(CLI::ExistingFile)
      ->required();

  extr_barr.add_option(
      "--extrusion-barrier-occupancy",
      c.extrusion_barrier_occupancy,
      "Probability that an extrusion barrier will be active (i.e. occupied) at any given time."
      "This parameter provides an easier mean to set --ctcf-occupied-probability-of-transition-to-self.")
      ->check(CLI::Range(0.0, 1.0));

  extr_barr.add_option(
      "--lef-bar-major-collision-prob",
      c.lef_bar_major_collision_pblock,
      "Collision probability of a LEF moving towards an extrusion barrier in blocking orientation.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  extr_barr.add_option(
      "--lef-bar-minor-collision-prob",
      c.lef_bar_minor_collision_pblock,
      "Collision probability of a LEF moving towards an extrusion barrier in non-blocking orientation.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  extr_barr.add_option(
      "--ctcf-occupied-probability-of-transition-to-self",
      c.ctcf_occupied_self_prob,
      "Probability that an extrusion barrier that was occupied during the previous iteration will also be "
      "occupied in the current iteration.")
      ->check(CLI::Range(0.0, 1.0));

  extr_barr.add_option(
      "--ctcf-not-occupied-probability-of-transition-to-self",
      c.ctcf_not_occupied_self_prob,
      "Probability that an extrusion barrier that was not occupied during the previous iteration will remain "
      "not occupied in the current iteration.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  extr_barr.add_flag(
      "--exclude-chrom-wo-barriers,!--keep-chrom-without-barriers",
      c.exclude_chrom_wo_extr_barriers,
      "Do not simulate loop extrusion on chromosomes without any extrusion barrier.")
      ->capture_default_str();


  rand.add_flag(
      "--randomize-contacts",
      c.randomize_contacts,
      "Enable randomization when collecting contacts.\n"
      "Contacts are randomized by drawing offsets from a genextreme distrubution.\n"
      "The distribution properties can be controlled through the parameters --mu, --sigma and --xi.")
      ->capture_default_str();

  rand.add_option(
      "--mu,--genextr-location",
      c.genextreme_mu,
      "Location parameter (mu) of the generalized extreme value distribution used to add noise to the contact matrix.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  rand.add_option(
      "--sigma,--genextr-scale",
      c.genextreme_sigma,
      "Scale parameter (sigma) of the generalized extreme value distribution used to add noise to the contact matrix.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  rand.add_option(
      "--xi,--genextr-shape",
      c.genextreme_xi,
      "Shape parameter (xi) of the generalized extreme value distribution used to add noise to the contact matrix.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  dbg.add_flag(
      "--skip-output",
      c.skip_output,
      "Don't write output files. Useful for profiling.")
      ->capture_default_str();

    gen.get_option("--target-contact-density")->excludes(gen.get_option("--number-of-iterations"));
  extr_barr.get_option("--extrusion-barrier-occupancy")->excludes("--ctcf-occupied-probability-of-transition-to-self");


  rand.get_option("--mu")->needs(rand.get_option("--randomize-contacts"));
  rand.get_option("--sigma")->needs(rand.get_option("--randomize-contacts"));
  rand.get_option("--xi")->needs(rand.get_option("--randomize-contacts"));
  // clang-format on
}

void Cli::make_simulation_subcommand() {
  auto& s = *this->_cli
                 .add_subcommand("simulate",
                                 "Perform a single genome-wide simulation and output the resulting "
                                 "contacts to a .cool file.")
                 ->fallthrough()
                 ->configurable();
  s.alias("sim");
  add_common_options(s, this->_config);

  auto& c = this->_config;
  auto& io = *s.get_option_group("Input/Output");
  auto& gen = *s.get_option_group("Generic");
  auto& dbg = *s.add_option_group("Debug/Profiling", "");

  // clang-format off
  io.add_option(
      "-o,--output-prefix",
      c.path_to_output_prefix,
      "Output prefix. Can be a full or relative path including the file name but without file extension.\n"
      "Example: -o /tmp/mymatrix will cause MoDLE to write contacts to a file named \"/tmp/mymatrix.cool\", "
      "while a log file will be saved at \"/tmp/mymatrix.log\".")
      ->required();

  io.add_flag(
      "--write-contacts-for-excluded-chroms",
      c.write_contacts_for_ko_chroms,
      "Write contacts for all chromosomes, even those where loop extrusion was not simulated.\n"
      "In the latter case MoDLE will only write to the chrom, bins and indexes datasets to the contact matrix.")
      ->capture_default_str();

  gen.add_option(
      "--ncells",
      c.num_cells,
      "Number of cells to simulate.\n"
      "Loop extrusion will be simulated independently on each chromosome across --ncells simulation instances.\n"
      "The total number of contacts for a given chromosome is obtained by aggregating contacts from all the "
      "simulation instances for that chromosome.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  dbg.add_flag(
      "--log-model-internal-state",
      c.log_model_internal_state,
      "Produce a detailed log of the model internal state throughout a simulation.\n"
      "Passing this flag often causes a noticeable reduction in MoDLE's throughput.")
      ->capture_default_str();

  s.get_option_group("Debug/Profiling")->get_option("--skip-output")->excludes(dbg.get_option("--log-model-internal-state"));
  // clang-format on
}

void Cli::make_perturbate_subcommand() {
  auto& s = *this->_cli.add_subcommand("perturbate", "TODO.")
                 ->configurable()
                 ->fallthrough();  // TODO add description
  s.alias("pert");
  add_common_options(s, this->_config);

  auto& c = this->_config;
  auto& io = *s.get_option_group("Input/Output");
  auto& gen = *s.get_option_group("Generic");

  // clang-format off
  io.add_option(
      "-o,--output-prefix",
      c.path_to_output_prefix,
      "Output prefix. Can be a full or relative path including the file name but without file extension.\n"
      "Example: -o /tmp/mycontacts will cause MoDLE to write interactions to a file named \"/tmp/mycontacts.bedpe.gz\", "
      "while a log file will be saved at \"/tmp/mycontacts.log\".")
      ->required();

  io.add_option(
      "--reference-contacts",
      c.path_to_reference_contacts,
      "Path to a Cooler file with the contact matrix to use as reference")
      ->check(CLI::ExistingFile)
      ->required();

  io.add_flag("--write-header,!--no-write-header",
      c.write_header,
      "Write header with column names to output file.")
      ->capture_default_str();

  io.add_option(
      "--feature-beds",
      c.path_to_feature_bed_files,
      "Path to one or more BED files containing features used to compute the total number of contacts between pairs of features.\n"
      "Pairs of features with a non-zero number of contacts will be written to a BEDPE file.\n"
      "When a single BED file is specified, the output will only contain within-feature contacts (Not yet implemented).")
      ->check(CLI::ExistingFile);

  gen.add_flag(
      "--generate-reference-matrix",
      c.compute_reference_matrix,
      "Compute and write to disk the reference contact matrix."
      "This is equivalent to running modle simulate followed by modle perturbate without changing parameters.")
      ->capture_default_str();

  gen.add_option(
      "--block-size",
      c.block_size,
      "Size of the block of pixels to use when generating contacts for a pair of features. Must be an odd number.")
      ->check(CLI::Range(1UL, std::numeric_limits<decltype(c.block_size)>::max()) |
              CLI::Validator(is_odd_number, "ODD-NUMBER", ""))
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  gen.add_option(
      "--deletion-size",
      c.deletion_size,
      "Size of deletion in bp. Used to enable/disable extrusion barriers and compute the total number of contacts between pairs of feats1.\n"
      "Specify 0 to compute all possible combinations. Ignored when --mode=\"cluster\".")
      ->check(CLI::NonNegativeNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  gen.add_option(
      "--deletion-list",
      c.path_to_deletion_bed,
      "Path to a BED file containing the list deletions to perform when perturbating extrusion barriers.")
      ->check(CLI::ExistingFile);

  gen.add_flag(
      "--write-tasks,!--no-write-taks",
      c.write_tasks_to_disk,
      "Write tasks to disk.")
      ->capture_default_str();
  // clang-format on

  gen.get_option("--deletion-size")->excludes("--deletion-list");
  gen.get_option("--deletion-list")->excludes("--deletion-size");
}

void Cli::make_replay_subcommand() {
  auto& s = *this->_cli.add_subcommand("replay", "TODO")->fallthrough();
  s.alias("rpl");

  auto& c = this->_config;
  auto& io = *s.add_option_group("Input/Output", "");
  auto& gen = *s.add_option_group("Generic", "");

  // clang-format off
  io.add_option(
      "--config",
      c.path_to_config_file,
      "Path to the config file of a previous run of modle perturbate.")
      ->check(CLI::ExistingFile)
      ->required();

  io.add_option(
      "--task-file",
      c.path_to_task_file,
      "Path to the task file produced by modle perturbate.")
      ->check(CLI::ExistingFile)
      ->required();

  io.add_option(
      "--task-filter-file",
      c.path_to_task_filter_file,
      "TSV containing the list of task IDs to be processed.")
      ->check(CLI::ExistingFile);

  io.add_option(
      "-o,--output-prefix",
      c.path_to_output_prefix,
      "Output prefix. Can be a full or relative path including the file name but without file extension.")
      ->required();

  io.add_flag(
      "-f,--force",
      c.force,
      "Force overwrite of output files if they exists.")
      ->capture_default_str();

  gen.add_option(
      "-t,--threads",
      c.nthreads,
      "Number of simulate_worker threads used to run the simulation.\n"
      "By default MoDLE will try to use all available threads.")
      ->check(CLI::PositiveNumber)
      ->transform(CLI::Bound(1U, std::thread::hardware_concurrency()))
      ->capture_default_str();
  // clang-format on
}

void Cli::make_cli() {
  this->_cli.description("Stochastic modeling of DNA loop extrusion.");
  this->_cli.set_version_flag("-V,--version", std::string{modle_version_long});
  this->_cli.require_subcommand(1);
  this->_cli.set_config("--config", "", "Path to MoDLE's config file.");

  this->make_simulation_subcommand();
  this->make_perturbate_subcommand();
  this->make_replay_subcommand();
}

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { this->make_cli(); }

const Config& Cli::parse_arguments() {
  using namespace std::string_view_literals;
  if (this->_cli.parsed()) {
    return this->_config;
  }

  this->_cli.name(this->_exec_name);
  this->_cli.parse(this->_argc, this->_argv);

  if (this->_cli.get_subcommand("simulate")->parsed()) {
    this->_subcommand = simulate;
  } else if (this->_cli.get_subcommand("perturbate")->parsed()) {
    this->_subcommand = pertubate;
  } else {
    assert(this->_cli.get_subcommand("replay")->parsed());
    this->_subcommand = replay;
    // The code in this branch basically tricks CLI11 into parsing a config file by creating a fake
    // argv like: {"modle", "pert", "--config", "myconfig.toml"}
    const auto config_backup = this->_config;
    constexpr std::string_view subcmd_arg = "pert\0"sv;
    constexpr std::string_view config_arg = "--config\0"sv;
    const std::array<const char*, 4> args{this->_exec_name.c_str(), subcmd_arg.data(),
                                          config_arg.data(),
                                          this->_config.path_to_config_file.c_str()};
    this->_cli.parse(args.size(), args.data());
    this->_config.path_to_config_file = config_backup.path_to_config_file;
    this->_config.path_to_task_file = config_backup.path_to_task_file;
    this->_config.path_to_task_filter_file = config_backup.path_to_task_filter_file;
    this->_config.path_to_output_prefix = config_backup.path_to_output_prefix;
    this->_config.nthreads = config_backup.nthreads;
  }

  this->validate_args();
  this->transform_args();
  this->_config.args = absl::MakeSpan(_argv, static_cast<usize>(_argc));

  this->_config.argv_json = this->to_json();
  return this->_config;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity) TODO: reduce complexity
std::string Cli::detect_path_collisions(modle::Config& c) const {
  std::string collisions;

  if (c.force || c.skip_output) {
    return "";
  }

  auto check_for_path_collisions = [](const boost::filesystem::path& path) -> std::string {
    assert(!path.empty());
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
  if (this->get_subcommand() == simulate && !c.path_to_model_state_log_file.empty() &&
      boost::filesystem::exists(c.path_to_model_state_log_file)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_model_state_log_file));
  }
  if (this->get_subcommand() != subcommand::replay &&
      boost::filesystem::exists(c.path_to_config_file)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_config_file));
  }
  if ((this->get_subcommand() == subcommand::simulate ||
       (this->get_subcommand() == subcommand::pertubate &&
        this->_config.compute_reference_matrix)) &&
      boost::filesystem::exists(c.path_to_output_file_cool)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_output_file_cool));
  }
  if ((this->get_subcommand() == subcommand::pertubate ||
       this->get_subcommand() == subcommand::replay) &&
      !c.path_to_output_file_bedpe.empty() &&
      boost::filesystem::exists(c.path_to_output_file_bedpe)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_output_file_bedpe));
  }
  if (this->get_subcommand() == subcommand::pertubate &&
      boost::filesystem::exists(this->_config.path_to_task_file)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_task_file));
  }

  return collisions;
}

int Cli::exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }

void Cli::validate_args() const {
  const auto& c = this->_config;
  std::vector<std::string> errors;

  if (c.burnin_smoothing_window_size > c.burnin_history_length) {
    assert(this->_cli.get_option("--burnin-smoothing-window-size"));
    assert(this->_cli.get_option("--burnin-history-length"));
    errors.emplace_back(fmt::format(
        FMT_STRING("The value passed to {} should be less or equal than that of {} ({} > {})"),
        "--burnin-smoothing-window-size", "--burnin-history-length", c.burnin_smoothing_window_size,
        c.burnin_history_length));
  }

  if (!c.path_to_reference_contacts.empty()) {
    try {
      modle::cooler::Cooler<>(c.path_to_reference_contacts, cooler::Cooler<>::IO_MODE::READ_ONLY,
                              c.bin_size);
    } catch (const std::runtime_error& e) {
      errors.emplace_back(e.what());
    }
  }

  if (!errors.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error(s) where encountered while validating CLI arguments:\n - {}"),
        fmt::join(errors, "\n - ")));
  }
}

void Cli::transform_args() {
  auto& c = this->_config;

  // Generate output file paths from output prefix
  c.path_to_output_file_cool = c.path_to_output_prefix;
  c.path_to_output_file_bedpe =
      !c.path_to_feature_bed_files.empty() ? c.path_to_output_prefix : boost::filesystem::path{};
  c.path_to_log_file = c.path_to_output_prefix;
  c.path_to_config_file = c.path_to_output_prefix;
  if (this->get_subcommand() == subcommand::pertubate) {
    c.path_to_task_file = c.path_to_output_prefix;
    c.path_to_task_file += "_tasks.tsv.gz";
  }
  if (this->get_subcommand() == subcommand::simulate) {
    c.path_to_model_state_log_file = c.path_to_output_prefix;
    c.path_to_model_state_log_file += "_internal_state_log.tsv.gz";
  }

  c.path_to_output_file_cool += ".cool";
  c.path_to_output_file_bedpe += !c.path_to_output_file_bedpe.empty() ? ".bedpe.gz" : "";
  c.path_to_log_file += ".log";
  c.path_to_config_file += "_config.toml";

  // Compute mean and std for rev and fwd extrusion speed
  if (auto& speed = c.rev_extrusion_speed; speed == (std::numeric_limits<bp_t>::max)()) {
    speed = static_cast<bp_t>(std::round(static_cast<double>(c.bin_size) / 2.0));
  }
  if (auto& speed = c.fwd_extrusion_speed; speed == (std::numeric_limits<bp_t>::max)()) {
    speed = static_cast<bp_t>(std::round(static_cast<double>(c.bin_size) / 2.0));
  }

  if (auto& stddev = c.fwd_extrusion_speed_std; stddev > 0 && stddev < 1) {
    stddev *= static_cast<double>(c.fwd_extrusion_speed);
  }
  if (auto& stddev = c.rev_extrusion_speed_std; stddev > 0 && stddev < 1) {
    stddev *= static_cast<double>(c.rev_extrusion_speed);
  }

  c.rev_extrusion_speed_burnin = static_cast<bp_t>(
      std::round(c.burnin_speed_coefficient * static_cast<double>(c.rev_extrusion_speed)));
  c.fwd_extrusion_speed_burnin = static_cast<bp_t>(
      std::round(c.burnin_speed_coefficient * static_cast<double>(c.fwd_extrusion_speed)));

  // Compute the transition probabilities for the HMM used to model extrusion barriers
  const auto pno = 1.0 - c.ctcf_not_occupied_self_prob;
  if (c.ctcf_occupied_self_prob == 0.0) {
    const auto occ = c.extrusion_barrier_occupancy;
    const auto pon = (pno - (occ * pno)) / occ;
    c.ctcf_occupied_self_prob = 1.0 - pon;
  } else {
    const auto pon = 1.0 - c.ctcf_occupied_self_prob;
    const auto occ = pno / (pno + pon);
    c.extrusion_barrier_occupancy = occ;
  }
}

Cli::subcommand Cli::get_subcommand() const { return this->_subcommand; }

void Cli::print_config(bool print_default_args) const {
  fmt::print(stderr, FMT_STRING("{}\n"), this->_cli.config_to_str(print_default_args, true));
}

void Cli::write_config_file(bool write_default_args) const {
  auto fp = fmt::output_file(this->_config.path_to_config_file.string());
  fp.print(FMT_STRING("{}\n"), this->_cli.config_to_str(write_default_args, true));
}

bool Cli::config_file_parsed() const { return !this->_cli.get_config_ptr()->empty(); }

std::string Cli::to_json() const {
  std::string buff;
  for (const auto& line : absl::StrSplit(this->_cli.config_to_str(true, false), '\n')) {
    if (line.empty()) {
      continue;
    }
    if (line.front() == '[') {  // Begin of the section for the active subcommand
      absl::StrAppend(&buff, line, "\n");
      continue;
    }
    // Given two subcommands named comm1 and comm2, assuming comm1 was parsed while comm2 was not,
    // the TOML produced by CLI11 will have values for comm2 formatted as comm2.myarg1=1,
    // comm2.myarg2="a" etc.
    // All we are doing here is to look for an argument name containing '.'.
    // In this way we can filter out entry corresponding to arguments for inactive subcommands
    const auto arg = line.substr(0, line.find('='));
    assert(!arg.empty());
    if (!absl::StrContains(arg, '.')) {
      absl::StrAppend(&buff, line, "\n");
    }
  }

  try {
    auto tt = toml::parse(buff);
    std::stringstream ss;
    ss << toml::json_formatter{tt};
    return ss.str();
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error occurred while converting MoDLE's config from TOML to JSON: {}"),
        e.what()));
  }
}

}  // namespace modle
