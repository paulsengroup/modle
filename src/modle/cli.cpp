// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
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

#include "modle/common/cli_utils.hpp"
#include "modle/common/common.hpp"             // for bp_t, i64
#include "modle/common/simulation_config.hpp"  // for Config
#include "modle/common/utils.hpp"              // for parse_numeric_or_throw
#include "modle/config/version.hpp"            // modle_version_long
#include "modle/cooler/cooler.hpp"             // for Cooler

namespace modle {

static std::string is_odd_number(std::string_view s) noexcept {
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

// These stream operators are needed to properly serialize enums when writing the config file to
// TOML or JSON
std::ostream& operator<<(std::ostream& os, const Config::StoppingCriterion& criterion) {
  os << Cli::stopping_criterion_map.at(criterion);
  return os;
}

std::ostream& operator<<(std::ostream& os, const Config::ContactSamplingStrategy strategy) {
  os << Cli::contact_sampling_strategy_map.at(strategy);
  return os;
}

static void add_common_options(CLI::App& subcommand, modle::Config& config) {
  auto& s = subcommand;
  auto& c = config;

  auto& io = *s.add_option_group("IO", "Options controlling MoDLE input, output and logs.");

  auto& lefbar =
      *s.add_option_group("Extrusion Barriers and Factors",
                          "Options controlling how extrusion barrier and LEFs are simulated.");

  auto& cgen = *s.add_option_group("Contact generation",
                                   "Options affecting contact sampling and registration.");

  auto& stopping = *s.add_option_group("Stopping criterion",
                                       "Options defining the simulation stopping criterion.");

  auto& misc = *s.add_option_group("Miscellaneous");

  auto& io_adv = *s.add_option_group("Advanced/IO",
                                     "Advanced options controlling MoDLE input, output and logs.");
  auto& lef_adv = *s.add_option_group("Advanced/Extrusion Factors",
                                      "Advanced options controlling how LEFs are simulated.");

  auto& barr_adv = *s.add_option_group("Advanced/Extrusion Barriers",
                                       "Advanced options controlling how extrusion are simulated.");

  auto& cgen_adv =
      *s.add_option_group("Advanced/Contact generation",
                          "Advanced options affecting contact sampling and registration.");

  auto& burnin_adv = *s.add_option_group("Advanced/Burn-in",
                                         "Options defining the simulation stopping criterion.");

  // clang-format off
  io.add_option(
      "-c,--chrom-sizes",
      c.path_to_chrom_sizes,
      "Path to file with chromosome sizes in chrom.sizes format.")
      ->check(CLI::ExistingFile)->required();

  io_adv.add_option(
      "--chrom-subranges",
      c.path_to_chrom_subranges,
      "Path to BED file with subranges of the chromosomes to simulate.")
      ->check(CLI::ExistingFile);

  io.add_option(
      "-b,--extrusion-barrier-file",
      c.path_to_extr_barriers,
      "Path to a file in BED6+ format with the genomic coordinates of extrusion barriers to be simulated.\n"
      "The score field in a BED record should be a number between 0 and 1 and is interpreted as the extrusion barrier occupancy for the extrusion barrier described by the record.\n"
      "Barriers mapping on chromosomes not listed in the chrom.sizes file passed through the --chrom-sizes option are ignored.")
      ->check(CLI::ExistingFile)
      ->required();

  io.add_flag(
      "-f,--force",
      c.force,
      "Overwrite existing files (if any).")
      ->capture_default_str();

  io.add_option(
      "-o,--output-prefix",
      c.path_to_output_prefix,
      "Output prefix.\n"
      "Can be an absolute or relative path including the file name but without the extension.\n"
      "Example: running modle sim -o /tmp/my_simulation ... yields the following files:\n"
      "         - /tmp/my_simulation.cool\n"
      "         - /tmp/my_simulation.log\n"
      "         - /tmp/my_simulation_config.toml")
      ->required();

  io.add_flag(
      "-q,--quiet",
      c.quiet,
      "Suppress console output to stderr.\n"
      "Only fatal errors will be logged to the console.\n"
      "Does not affect entries written to the log file.")
      ->capture_default_str();

  io_adv.add_flag(
      "--log-model-internal-state",
      c.log_model_internal_state,
      fmt::format(FMT_STRING(
                      "Collect detailed statistics regarding the internal state of MoDLE simulation instance(s).\n"
                      "Statistics will be written to a compressed file under the prefix specified through the --output-prefix option.\n"
                      "Example: with --output-prefix=/tmp/myprefix, statistics will be written to file /tmp/myprefix_internal_state.log.gz.\n"
                      "Depending on the input file(s) and parameters, specifying this option may hinder simulation throughput.\n"
                      "Currently the following metrics are collected:\n"
                      " - {}"),
                  fmt::join(absl::StrSplit(Config::model_internal_state_log_header, '\t'), "\n - ")))
      ->capture_default_str();

  io_adv.add_flag(
      "--simulate-chromosomes-wo-barriers,!--skip-chromosomes-wo-barriers",
      c.simulate_chromosomes_wo_barriers,
      "Enable/disable simulation of loop extrusion for chromosomes with 0 extrusion barriers.\n"
      "When --skip-chromosomes-wo-barriers is passed, entries for chromosomes without barriers will "
      "still be written to the output .cool file, but no contacts will be generated for those chromosomes.")
      ->capture_default_str();

  io_adv.add_flag(
      "--skip-output",
      c.skip_output,
      "Do not write output files. Mostly useful for profiling.")
      ->capture_default_str();

  lefbar.add_option(
      "--lef-density,--lefs-per-mbp",
      c.number_of_lefs_per_mbp,
      "Loop extrusion factor (LEF) density expressed as the number of LEF per Mbp of DNA simulated.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  lefbar.add_option(
      "--avg-lef-processivity",
      c.avg_lef_processivity,
      "Average LEF processivity in bp.\n"
      "The average LEF processivity corresponds to the average size of loops extruded by unobstructed LEFs.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  lefbar.add_option(
      "--probability-of-lef-bypass",
      c.probability_of_extrusion_unit_bypass,
      "Probability that two colliding LEFs will avoid collision by bypassing each other.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  lefbar.add_option(
      "--extrusion-barrier-occupancy",
      c.extrusion_barrier_occupancy,
      "Probability that an extrusion barrier is occupied (i.e. blocking) at any given time.\n"
      "This parameter provides an easier mean to set --extrusion-barrier-bound-stp.\n"
      "Passing this parameter will override barrier occupancies read from the BED file specified through --extrusion-barrier-file.")
      ->check(CLI::Range(0.0, 1.0));

  lef_adv.add_option(
      "--hard-stall-lef-stability-multiplier",
      c.hard_stall_lef_stability_multiplier,
      "Coefficient to control the DNA-binding stability of LEFs that are stalled on both sides by a pair of "
      "extrusion barriers in convergent orientation.\n"
      "Setting this to 1 makes the DNA-binding stability of hard-stalled LEFs identical to that of unobstructed LEFs.\n"
      "Has no effects when --lef-bar-major-collision-prob=0.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  lef_adv.add_option(
      "--soft-stall-lef-stability-multiplier",
      c.soft_stall_lef_stability_multiplier,
      "Coefficient to control the DNA-binding stability of LEFs stalled by one or more extrusion barriers "
      "in a non-blocking orientation.\n"
      "Setting this to 1 makes the DNA-binding stability of soft-stalled LEFs identical to that of unobstructed LEFs.\n"
      "Has no effects when --lef-bar-minor-collision-prob=0.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  lef_adv.add_option(
      "--fwd-extrusion-speed",
      c.fwd_extrusion_speed,
      "Average extrusion speed expressed in bp/epoch for forward-moving extrusion units.\n"
      "Extrusion units are assigned a candidate moving distance at the beginning of every epoch.\n"
      "Moves are sampled from a normal distribution centered around --fwd-extrusion-speed and with --fwd-extrusion-speed-std as its standard deviation.\n"
      "Candidate moves represent the maximum distance a given extrusion unit is set to travel during the current epoch.\n"
      "Moving distances can be shortened by collision events taking place during the current epoch.\n"
      "By deafult extrusion speed is set to half the bin size specified through the --resolution option.")
      ->transform(utils::str_float_to_str_int)
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  lef_adv.add_option(
      "--rev-extrusion-speed",
      c.rev_extrusion_speed,
      "Same as --fwd-extrusion-speed but for extrusion units moving in 3'-5' direction.")
      ->transform(utils::str_float_to_str_int)
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  lef_adv.add_option(
      "--fwd-extrusion-speed-std",
      c.fwd_extrusion_speed_std,
      "Standard deviation of the normal distribution used to sample candidate moves for extrusion units moving in 5'-3' direction.\n"
      "See help message for --fwd-extrusion-speed for more details regarding candidate moves.\n"
      "Specifying --fwd-extrusion-speed-std=0 will make candidate moves deterministic.\n"
      "Standard deviations between 0 and 1 are interpreted a percentage of the average extrusion speed.\n"
      "Example: when running modle sim --fwd-extrusion-speed=10000 ...\n"
      "         --fwd-extrusion-speed-std=0.1 and --fwd-extrusion-speed-std=1000 are equivalent.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  lef_adv.add_option(
      "--rev-extrusion-speed-std",
      c.rev_extrusion_speed_std,
      "Same as --fwd-extrusion-speed-std but for extrusion units moving in 3'-5' direction.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  barr_adv.add_option(
      "--lef-bar-major-collision-prob",
      c.lef_bar_major_collision_pblock,
      "Probability of collision between LEFs and extrusion barriers where barriers are pointing towards the extrusion direction.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  barr_adv.add_option(
      "--lef-bar-minor-collision-prob",
      c.lef_bar_minor_collision_pblock,
      "Probability of collision between LEFs and extrusion barriers where barriers are pointing away from the extrusion direction.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  barr_adv.add_option(
      "--extrusion-barrier-bound-stp",
      c.ctcf_occupied_self_prob,
      "Self-transition probability for extrusion barriers in the \"bound\" state.\n"
      "In other words, the probability that an extrusion barrier that is active in the current epoch will remain active during the next epoch.")
      ->check(CLI::Range(0.0, 1.0));

  barr_adv.add_option(
      "--extrusion-barrier-not-bound-stp",
      c.ctcf_not_occupied_self_prob,
      "Self-transition probability for extrusion barriers in the \"not-bound\" state.\n"
      "In other words, the probability that an extrusion barrier that is inactive in the current epoch will remain inactive during the next epoch.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  cgen.add_option(
      "--contact-sampling-strategy",
      c.contact_sampling_strategy,
      fmt::format(FMT_STRING("Strategy to use when sampling contacts. Should be one of {}.\n"
                             "When one of the *-with-noise strategies is specified, contacts are"
                             "randomized by applying a random offset to the location of LEF extrusion units.\n"
                             "Offsets are drawn from a genextreme distrubution.\n"
                             "The distribution parameters can be controlled through the options --mu, --sigma and --xi."),
                  utils::format_collection_to_english_list(Cli::contact_sampling_strategy_map.keys_view(), ", ", " or ")))
      ->transform(CLI::CheckedTransformer(Cli::contact_sampling_strategy_map))
      ->capture_default_str();

  cgen.add_option(
      "--lef-fraction-for-contact-sampling",
      c.lef_fraction_contact_sampling,
      "Fraction of LEFs to use for contact sampling.\n"
      "The actual number of LEFs to use for contact sampling in a given epoch is drawn from a Poisson distribution "
      "with lambda=lf*|lefs|, where --lef-fraction-for-contact-sampling=lf and |lefs| is equal to total number of LEFs "
      "in the current simulation instance.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  cgen.add_option(
      "-r,--resolution",
      c.bin_size,
      "Resolution in base-pairs for the output contact matrix in cooler format.\n"
      "NOTE: MoDLE simulation always take place at 1 bp resolution.\n"
      "      This parameter only affects the resolution of the output contact matrix.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  cgen.add_option(
      "-w,--diagonal-width",
      c.diagonal_width,
      "Width of the subdiagonal window for the in-memory contact matrix.\n"
      "This setting affects the maximum distance of a pair of bins whose interactions will be tracked by MoDLE.\n"
      "As a rule of thumb, --diagonal-width should roughly 10x the average LEF processivity\n."
      "Setting --diagonal-width to very large values (i.e. tens of Mbp) will significantly increase MoDLE's memory requirements.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  cgen_adv.add_option(
      "--num-of-tad-contacts-per-sampling-event",
      c.number_of_tad_contacts_per_sampling_event,
      "Number of TAD contacts to sample each sampling event.\n"
      "Setting --num-of-tad-contacts-per-sampling-event=0 will yield a contact matrix where contacts are concentrated on stripes and dots.")
      ->capture_default_str();

  cgen_adv.add_option(
      "--num-of-loop-contacts-per-sampling-event",
      c.number_of_loop_contacts_per_sampling_event,
      "Number of LEF-mediated contacts to sample each sampling event.\n"
      "Setting --num-of-loop-contacts-per-sampling-event=0 will yield a contact matrix without stripes and dots.")
      ->capture_default_str();

  cgen_adv.add_option(
      "--mu,--genextr-location",
      c.genextreme_mu,
      "Location parameter (mu) of the generalized extreme value distribution used to add noise to molecular contacts.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  cgen_adv.add_option(
      "--sigma,--genextr-scale",
      c.genextreme_sigma,
      "Scale parameter (sigma) of the generalized extreme value distribution used to add noise to molecular contacts.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  cgen_adv.add_option(
      "--xi,--genextr-shape",
      c.genextreme_xi,
      "Shape parameter (xi) of the generalized extreme value distribution used to add noise to molecular contacts.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  stopping.add_option(
      "-s,--stopping-criterion",
      c.stopping_criterion,
      fmt::format(FMT_STRING("Simulation stopping criterion. Should be one of {}."),
                  utils::format_collection_to_english_list(Cli::stopping_criterion_map.keys_view(), ", ", " or ")))
      ->transform(CLI::CheckedTransformer(Cli::stopping_criterion_map))
      ->capture_default_str();

  stopping.add_option(
      "--target-number-of-epochs",
      c.target_simulation_epochs,
      "Target number of epochs to be simulated.\n"
      "Each simulation instance will run exactly --target-number-of-epochs epochs after burn-in phase.\n"
      "Has no effect when --stopping-criterion=\"contact-density\".")
      ->check(CLI::PositiveNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  stopping.add_option(
      "--target-contact-density",
      c.target_contact_density,
      "Average contact density to be reached before the simulation is halted.\n"
      "The target contact density applies independently to each chromosome.\n"
      "Example: modle sim --chrom-sizes hg38.chr1.chrom.sizes \\\n"
      "                   --diagonal-width 3000000            \\\n"
      "                   --resolution 100000                 \\\n"
      "                   --target-contact-density 2\n"
      "         Assuming chr1 is exactly 250 Mbp long, MoDLE will schedule simulation instances to yield a total of 150 000 contacts for chr1:\n"
      "         2 * (250e6 / 100e3) * (3e6 / 100e3) = 150e3.")
      ->check(CLI::NonNegativeNumber);

  misc.add_option(
      "--ncells",
      c.num_cells,
      "Number of simulation instances or cells to simulate.\n"
      "Loop extrusion will be simulated independently for every chromosome across --ncells simulation instances.\n"
      "The final contact matrix is produced by accumulating contacts generated across all simulation instances.\n"
      "To achieve good performance and CPU utilization we recommend setting --ncells equal to the number of available CPU cores.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  misc.add_option(
      "-t,--threads",
      c.nthreads,
      "Number of worker threads used to run simulation instances.\n"
      "By default MoDLE will spawn a number of worker threads equal to the number of logical CPU cores.\n"
      "Sometimes slightly better performance can be obtained by setting --threads equal to the number of physical CPU cores.\n")
      ->check(CLI::PositiveNumber)
      ->transform(CLI::Bound(1U, std::thread::hardware_concurrency()))
      ->capture_default_str();

  misc.add_option(
      "--seed",
      c.seed,
      "Base seed to use for random number generation.")
      ->check(CLI::NonNegativeNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  burnin_adv.add_flag(
      "--skip-burnin",
      c.skip_burnin,
      "Skip the burn-in phase and start collecting contacts from the first extrusion round.")
      ->capture_default_str();

  burnin_adv.add_option(
      "--burnin-target-epochs-for-lef-activation",
       c.burnin_target_epochs_for_lef_activation,
      "Number of epochs over which LEFs are progressively activated and bound to DNA.\n"
      "Note: this number is approximate, as LEFs are activated using a Possion process.\n"
      "By default this parameter is computed based on the average LEF processivity, overall extrusion speed "
      "and burn-in extrusion speed coefficient (--avg-lef-processivity, --fwd/rev-extrusion-speed "
      "and --burn-in-extr-speed-coefficient respectively).")
       ->check(CLI::PositiveNumber)
       ->capture_default_str();

  burnin_adv.add_option(
      "--burnin-history-length",
      c.burnin_history_length,
      "Number of epochs used to determine whether a simulation instance has reached a stable state.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  burnin_adv.add_option(
      "--burnin-smoothing-window-size",
      c.burnin_smoothing_window_size,
      "Window size used to smooth metrics monitored to decide when the burn-in phase should be terminated.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  burnin_adv.add_option(
      "--max-burnin-epochs",
       c.max_burnin_epochs,
       "Upper bound for the burn-in phase duration.\n"
       "This is especially useful when simlating loop extrusion with very few LEFs.\n"
       "When this is the case, --max-burnin-epochs=100000 can be used as a very conservative threshold.")
       ->check(CLI::PositiveNumber)
       ->capture_default_str();

  burnin_adv.add_option(
      "--burnin-extr-speed-coefficient",
      c.burnin_speed_coefficient,
      "Extrusion speed coefficient to apply during the burn-in phase.\n"
      "Setting this to numbers > 1.0 will speed-up the burn-in phase, as the average loop size will stabilize faster.\n"
      "This parameter has many gotchas. For the time being we don't recommend changing this parameter.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  // Address option dependencies/incompatibilities
  io_adv.get_option("--skip-output")->excludes(io_adv.get_option("--log-model-internal-state"));
  stopping.get_option("--target-contact-density")->excludes(stopping.get_option("--target-number-of-epochs"));
  lefbar.get_option("--extrusion-barrier-occupancy")->excludes(barr_adv.get_option("--extrusion-barrier-bound-stp"));
  // clang-format on
}

void Cli::make_simulation_subcommand() {
  auto& s =
      *this->_cli
           .add_subcommand(
               "simulate",
               "Simulate loop extrusion and write resulting molecular contacts in a .cool file.")
           ->fallthrough()
           ->configurable();
  s.alias("sim");
  add_common_options(s, this->_config);
}

void Cli::make_perturbate_subcommand() {
  // TODO add description
  // clang-format off
  auto& s = *this->_cli.add_subcommand("perturbate", "TODO.")
                 ->configurable()
                 ->fallthrough()
                 ->group("");
  // clang-format on

  s.alias("pert");
  add_common_options(s, this->_config);

  auto& c = this->_config;
  auto& io = *s.get_option_group("IO");
  auto& misc = *s.get_option_group("Miscellaneous");

  auto& io_adv = *s.get_option_group("Advanced/IO");

  // Remove unused flags/options
  io_adv.remove_option(io_adv.get_option("--simulate-chromosomes-wo-barriers"));
  misc.remove_option(misc.get_option("--ncells"));

  // Update flag/option descriptions
  // clang-format off
  io.get_option("--output-prefix")
      ->description(
      "Output prefix.\n"
      "Can be a full or relative path including the file name but without extension.\n"
      "Example: -o /tmp/mymatrix will produce the following files:\n"
      "         - /tmp/mymatrix.bedpe.gz\n"
      "         - /tmp/mymatrix.log\n");

  // Add new flags/options
  io.add_option(
      "--reference-contacts",
      c.path_to_reference_contacts,
      "Path to a cooler file with the contact matrix to use as reference")
      ->check(CLI::ExistingFile)
      ->required();

  io_adv.add_flag("--write-header,!--no-write-header",
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

  io_adv.add_flag(
      "--generate-reference-matrix",
      c.compute_reference_matrix,
      "Compute and write to disk the reference contact matrix."
      "This is equivalent to running modle simulate followed by modle perturbate without changing parameters.")
      ->capture_default_str();

  misc.add_option(
      "--block-size",
      c.block_size,
      "Size of the block of pixels to use when generating contacts for a pair of features. Must be an odd number.")
      ->check(CLI::Range(1UL, std::numeric_limits<decltype(c.block_size)>::max()) |
              CLI::Validator(is_odd_number, "ODD-NUMBER", ""))
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  misc.add_option(
      "--deletion-size",
      c.deletion_size,
      "Size of deletion in bp. Used to enable/disable extrusion barriers and compute the total number of contacts between pairs of feats1.\n"
      "Specify 0 to compute all possible combinations. Ignored when --mode=\"cluster\".")
      ->check(CLI::NonNegativeNumber)
      ->transform(utils::str_float_to_str_int)
      ->capture_default_str();

  misc.add_option(
      "--deletion-list",
      c.path_to_deletion_bed,
      "Path to a BED file containing the list deletions to perform when perturbating extrusion barriers.")
      ->check(CLI::ExistingFile);
  // clang-format on

  misc.get_option("--deletion-size")->excludes("--deletion-list");
}

void Cli::make_replay_subcommand() {
  // TODO add description
  // clang-format off
  auto& s = *this->_cli.add_subcommand("replay", "TODO")
                 ->fallthrough()
                 ->group("");
  // clang-format on
  s.alias("rpl");

  auto& c = this->_config;
  auto& io = *s.add_option_group("IO", "Options controlling MoDLE input, output and logs.");
  auto& various = *s.add_option_group("Various");

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
      "Output prefix.\n"
      "Can be a full or relative path including the file name but without extension.")
      ->required();

  io.add_flag(
      "-f,--force",
      c.force,
      "Force overwrite of output files if they exists.")
      ->capture_default_str();

  various.add_option(
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
  this->_cli.description(
      "High-performance stochastic modeling of DNA loop extrusion interactions.");
  this->_cli.set_version_flag("-V,--version", std::string{modle::config::version::str_long()});
  this->_cli.require_subcommand(1);
  this->_cli.set_config("--config", "", "Path to MoDLE's config file (optional).", false);

  this->make_simulation_subcommand();
  this->make_perturbate_subcommand();
  this->make_replay_subcommand();
}

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { this->make_cli(); }

const Config& Cli::parse_arguments() {
  if (this->_cli.parsed()) {
    return this->_config;
  }

  this->_cli.name(this->_exec_name);
  this->_cli.parse(this->_argc, this->_argv);

  if (this->_cli.get_subcommand("simulate")->parsed()) {
    this->_subcommand = simulate;
  } else if (this->_cli.get_subcommand("perturbate")->parsed()) {
    this->_subcommand = perturbate;
  } else {
    assert(this->_cli.get_subcommand("replay")->parsed());
    this->_subcommand = replay;
    // The code in this branch basically tricks CLI11 into parsing a config file by creating a fake
    // argv like: {"modle", "pert", "--config", "myconfig.toml"}
    const auto config_backup = this->_config;
    constexpr std::string_view subcmd_arg = "pert\0";
    constexpr std::string_view config_arg = "--config\0";
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
       (this->get_subcommand() == subcommand::perturbate &&
        this->_config.compute_reference_matrix)) &&
      boost::filesystem::exists(c.path_to_output_file_cool)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_output_file_cool));
  }
  if ((this->get_subcommand() == subcommand::perturbate ||
       this->get_subcommand() == subcommand::replay) &&
      !c.path_to_output_file_bedpe.empty() &&
      boost::filesystem::exists(c.path_to_output_file_bedpe)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_output_file_bedpe));
  }
  if (this->get_subcommand() == subcommand::perturbate &&
      boost::filesystem::exists(this->_config.path_to_task_file)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_task_file));
  }

  return collisions;
}

int Cli::exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }

void Cli::validate_args() const {
  const auto& c = this->_config;
  std::vector<std::string> errors;

  const auto* subcmd = this->get_subcommand_ptr();

  if (c.burnin_smoothing_window_size > c.burnin_history_length) {
    assert(
        subcmd->get_option_group("Advanced/Burn-in")->get_option("--burnin-smoothing-window-size"));
    assert(subcmd->get_option_group("Advanced/Burn-in")->get_option("--burnin-history-length"));
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

  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  if (!(c.contact_sampling_strategy & Config::ContactSamplingStrategy::noisify)) {
    const auto& group = subcmd->get_option_group("Advanced/Contact generation");
    for (const std::string label : {"--mu", "--sigma", "--xi"}) {
      if (!group->get_option(label)->empty()) {
        errors.emplace_back(fmt::format(FMT_STRING("Option {} requires the strategy passed to "
                                                   "--contact-sampling-strategy to be one of "
                                                   "the *-with-noise strategies."),
                                        label));
      }
    }
  }

  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  if (c.contact_sampling_strategy & Config::ContactSamplingStrategy::tad &&
      !subcmd->get_option_group("Advanced/Contact generation")
           ->get_option("--num-of-tad-contacts-per-sampling-event")
           ->empty()) {
    errors.emplace_back(
        "--num-of-tad-contacts-per-sampling-event is not allowed when the sampling strategy passed "
        "through --contact-sampling-strategy does not involve sampling contacts from TADs.");
  }

  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  if (c.contact_sampling_strategy & Config::ContactSamplingStrategy::loop &&
      !subcmd->get_option_group("Advanced/Contact generation")
           ->get_option("--num-of-loop-contacts-per-sampling-event")
           ->empty()) {
    errors.emplace_back(
        "--num-of-loop-contacts-per-sampling-event is not allowed when the sampling strategy "
        "passed through --contact-sampling-strategy does not involve sampling contacts from "
        "loops.");
  }

  if (c.number_of_tad_contacts_per_sampling_event + c.number_of_loop_contacts_per_sampling_event ==
      0) {
    // clang-format off
    assert(subcmd->get_option_group("Advanced/Contact generation")->get_option("--num-of-tad-contacts-per-sampling-event"));
    assert(subcmd->get_option_group("Advanced/Contact generation")->get_option("--num-of-loop-contacts-per-sampling-event"));
    // clang-format on
    errors.emplace_back(
        "--num-of-tad-contacts-per-sampling-event and --num-of-loop-contacts-per-sampling-event "
        "cannot both be set to zero.");
  }

  if (!errors.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error(s) where encountered while validating CLI arguments:\n - {}"),
        fmt::join(errors, "\n - ")));
  }
}

CLI::App* Cli::get_subcommand_ptr() {
  for (const auto& subcmd : {"sim", "pert", "replay"}) {
    if (auto* s = this->_cli.get_subcommand(subcmd); s->parsed()) {
      return s;
    }
  }
  MODLE_UNREACHABLE_CODE;
}

const CLI::App* Cli::get_subcommand_ptr() const {
  for (const auto& subcmd : {"sim", "pert", "replay"}) {
    if (auto* s = this->_cli.get_subcommand(subcmd); s->parsed()) {
      return s;
    }
  }
  MODLE_UNREACHABLE_CODE;
}

void Cli::transform_args() {
  auto& c = this->_config;

  // Generate output file paths from output prefix
  c.path_to_output_file_cool = c.path_to_output_prefix;
  c.path_to_output_file_bedpe =
      !c.path_to_feature_bed_files.empty() ? c.path_to_output_prefix : boost::filesystem::path{};
  c.path_to_log_file = c.path_to_output_prefix;
  c.path_to_config_file = c.path_to_output_prefix;
  if (this->get_subcommand() == subcommand::perturbate) {
    c.path_to_task_file = c.path_to_output_prefix;
    c.path_to_task_file += "_tasks.tsv.gz";
  }
  if (this->get_subcommand() == subcommand::simulate) {
    c.path_to_model_state_log_file = c.path_to_output_prefix;
    c.path_to_model_state_log_file += "_internal_state.log.gz";
  }

  c.path_to_output_file_cool += ".cool";
  c.path_to_output_file_bedpe += !c.path_to_output_file_bedpe.empty() ? ".bedpe.gz" : "";
  c.path_to_log_file += ".log";
  c.path_to_config_file += "_config.toml";

  // Compute mean and std for rev and fwd extrusion speed
  if (auto& speed = c.rev_extrusion_speed; speed == (std::numeric_limits<bp_t>::max)()) {
    speed = (c.bin_size + 1) / 2;
  }
  if (auto& speed = c.fwd_extrusion_speed; speed == (std::numeric_limits<bp_t>::max)()) {
    speed = (c.bin_size + 1) / 2;
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

  // Compute the probability of LEF release from avg. processivity and two-sided extr. speed
  c.prob_of_lef_release = static_cast<double>(c.rev_extrusion_speed + c.fwd_extrusion_speed) /
                          static_cast<double>(c.avg_lef_processivity);
  c.prob_of_lef_release_burnin =
      static_cast<double>(c.rev_extrusion_speed_burnin + c.fwd_extrusion_speed_burnin) /
      static_cast<double>(c.avg_lef_processivity);

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

  const auto* subcmd = this->get_subcommand_ptr();

  if (!subcmd->get_option_group("Extrusion Barriers and Factors")
           ->get_option("--extrusion-barrier-occupancy")
           ->empty()) {
    c.override_extrusion_barrier_occupancy = true;
  }
}

Cli::subcommand Cli::get_subcommand() const { return this->_subcommand; }

void Cli::print_config(bool print_default_args) const {
  fmt::print(stderr, FMT_STRING("{}\n"), this->_cli.config_to_str(print_default_args, true));
}

void Cli::write_config_file(bool write_default_args) const {
  auto fp = fmt::output_file(this->_config.path_to_config_file.string());
  fp.print(FMT_STRING("# Config created by {} on {}\n{}\n"), modle::config::version::str_long(),
           absl::FormatTime(absl::Now(), absl::UTCTimeZone()),
           this->_cli.config_to_str(write_default_args, true));
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
