// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./cli.hpp"

#include <absl/strings/str_cat.h>    // for StrAppend
#include <absl/strings/str_split.h>  // for SplitIterator, Splitter, StrSplit, operator!=
#include <absl/time/clock.h>
#include <absl/types/span.h>  // for MakeSpan
#include <fmt/format.h>       // for format, FMT_STRING, join, print, make_format_...
#include <fmt/os.h>           // for output_file, ostream
#include <fmt/std.h>
#include <spdlog/spdlog.h>
#include <toml++/toml.h>  // for array::operator[], operator<<, parse, print_to...

#include <CLI/CLI.hpp>  // for Option_group, App
#include <algorithm>    // for max
#include <array>        // for array
#include <cassert>      // for assert
#include <cmath>        // for round, pow, log
#include <coolerpp/coolerpp.hpp>
#include <cstdio>      // for stderr
#include <exception>   // for exception
#include <filesystem>  // for path, operator<<
#include <limits>      // for numeric_limits
#include <sstream>     // for streamsize, stringstream, basic_ostream
#include <stdexcept>   // for invalid_argument, out_of_range, runtime_error
#include <string>      // for allocator, string, basic_string
#include <thread>      // for hardware_concurrency
#include <vector>      // for vector

#include "modle/common/cli_utils.hpp"
#include "modle/common/common.hpp"  // for bp_t, i64
#include "modle/common/fmt_helpers.hpp"
#include "modle/common/simulation_config.hpp"  // for Config
#include "modle/common/utils.hpp"              // for parse_numeric_or_throw
#include "modle/config/version.hpp"            // modle_version_long

namespace modle {

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

static std::vector<CLI::App*> add_common_options(CLI::App& subcommand, modle::Config& config) {
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

  auto& adv = *s.add_option_group("Advanced");

  auto& io_adv =
      *adv.add_option_group("IO", "Advanced options controlling MoDLE input, output and logs.");
  auto& lef_adv = *adv.add_option_group("Extrusion Factors",
                                        "Advanced options controlling how LEFs are simulated.");

  auto& barr_adv = *adv.add_option_group(
      "Extrusion Barriers", "Advanced options controlling how extrusion are simulated.");

  auto& cgen_adv = *adv.add_option_group(
      "Contact generation", "Advanced options affecting contact sampling and registration.");

  auto& burnin_adv =
      *adv.add_option_group("Burn-in", "Options defining the simulation stopping criterion.");

  auto& misc_adv = *adv.add_option_group("Miscellaneous");

  auto& deprecated = *s.add_option_group("");

  // clang-format off
  io.add_option(
      "-c,--chrom-sizes",
      c.path_to_chrom_sizes,
      "Path to file with chromosome sizes in chrom.sizes format.")
      ->check(CLI::ExistingFile)->required();

  io.add_option(
      "-b,--extrusion-barrier-file",
      c.path_to_extr_barriers,
      "Path to a file in BED6+ format with the genomic coordinates of extrusion barriers to be\n"
      "simulated. The score field in a BED record should be a number between 0 and 1 and is\n"
      "interpreted as the extrusion barrier occupancy for the extrusion barrier described\n"
      "by the record.\n"
      "Barriers mapping on chromosomes not listed in the chrom.sizes file passed through\n"
      "the --chrom-sizes option are ignored.")
      ->check(CLI::ExistingFile)
      ->required();

  io.add_option(
      "-g,--genomic-intervals",
      c.path_to_genomic_intervals,
      "Path to BED3+ file with the genomic regions to be simulated.\n"
      "Intervals listed in this file should be a subset of the chromosomes defined in the .chrom.sizes files.\n"
      "Intervals referring to chromosomes not listed in the .chrom.sizes file are ignored.")
      ->check(CLI::ExistingFile);

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

  io.add_option(
      "--assembly-name",
      c.assembly_name,
      "Name of the genome assembly to be simulated.\n"
      "This is only used to populate the \"assembly\" attribute in the output .cool file.")
      ->transform([](std::string str) { return str.empty() ? "unknown" : str; }, "", "")
      ->capture_default_str();

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
                      "Statistics will be written to a compressed file under the prefix specified through the\n"
                      "--output-prefix option.\n"
                      "Example: modle sim --output-prefix=/tmp/myprefix\n"
                               "statistics will be written to file /tmp/myprefix_internal_state.log.gz.\n"
                      "Depending on the input file(s) and parameters, specifying this option may hinder\n"
                      "simulation throughput. Currently the following metrics are collected:\n"
                      " - {}"),
                  fmt::join(absl::StrSplit(Config::model_internal_state_log_header, '\t'), "\n - ")))
      ->capture_default_str();

  io_adv.add_flag(
      "--simulate-chromosomes-wo-barriers,!--skip-chromosomes-wo-barriers",
      c.simulate_chromosomes_wo_barriers,
      "Enable/disable simulation of loop extrusion for chromosomes with 0 extrusion barriers.\n"
      "When --skip-chromosomes-wo-barriers is passed, entries for chromosomes without barriers\n"
      "will still be written to the output .cool file, but no contacts will be generated for\n"
      "those chromosomes.")
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
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  lefbar.add_option(
      "--avg-lef-processivity",
      c.avg_lef_processivity,
      "Average LEF processivity in bp.\n"
      "The average LEF processivity corresponds to the average size of loops extruded by\n"
      "unobstructed LEFs.")
      ->transform(utils::cli::AsGenomicDistance)
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
      "Passing this parameter will override barrier occupancies read from the BED file specified\n"
      "through --extrusion-barrier-file.")
      ->check(CLI::Range(0.0, 1.0));

  lefbar.add_flag(
      "--track-1d-lef-position,!--no-track-1d-lef-position",
      c.track_1d_lef_position,
      "Toggle on/off tracking of LEF positions in 1D space.\n"
      "LEF positions are aggregated across all simulated cells and are written to disk in BigWig format\n"
      "under the prefix specified through --output-prefix.\n")
      ->capture_default_str();

  lef_adv.add_option(
      "--hard-stall-lef-stability-multiplier",
      c.hard_stall_lef_stability_multiplier,
      "Coefficient to control the DNA-binding stability of LEFs that are stalled on both sides by a\n"
      "pair of extrusion barriers in convergent orientation.\n"
      "Setting this to 1 makes the DNA-binding stability of hard-stalled LEFs identical to that of\n"
      "unobstructed LEFs. Has no effects when --lef-bar-major-collision-prob=0.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  lef_adv.add_option(
      "--soft-stall-lef-stability-multiplier",
      c.soft_stall_lef_stability_multiplier,
      "Coefficient to control the DNA-binding stability of LEFs stalled by one or more extrusion\n"
      "barriers in a non-blocking orientation. Setting this to 1 makes the DNA-binding stability\n"
      "of soft-stalled LEFs identical to that of unobstructed LEFs.\n"
      "Has no effects when --lef-bar-minor-collision-prob=0.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  lef_adv.add_option(
      "--fwd-extrusion-speed",
      c.fwd_extrusion_speed,
      "Average extrusion speed expressed in bp/epoch for forward-moving extrusion units.\n"
      "Extrusion units are assigned a candidate moving distance at the beginning of every epoch.\n"
      "Moves are sampled from a normal distribution centered around --fwd-extrusion-speed and with\n"
      "--fwd-extrusion-speed-std as its standard deviation. Candidate moves represent the maximum\n"
      "distance a given extrusion unit is set to travel during the current epoch.\n"
      "Moving distances can be shortened by collision events taking place during the current epoch.\n"
      "By deafult extrusion speed is set to half the bin size specified through the --resolution\n"
      "option.")
      ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit |  utils::cli::AsGenomicDistance)
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  lef_adv.add_option(
      "--rev-extrusion-speed",
      c.rev_extrusion_speed,
      "Same as --fwd-extrusion-speed but for extrusion units moving in 3'-5' direction.")
      ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit |  utils::cli::AsGenomicDistance)
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  lef_adv.add_option(
      "--fwd-extrusion-speed-std",
      c.fwd_extrusion_speed_std,
      "Standard deviation of the normal distribution used to sample candidate moves for extrusion\n"
      "units moving in 5'-3' direction. See help message for --fwd-extrusion-speed for more details\n"
      "regarding candidate moves.\n"
      "Specifying --fwd-extrusion-speed-std=0 will make candidate moves deterministic.\n"
      "Standard deviations between 0 and 1 are interpreted a percentage of the average extrusion\n"
      "speed.\n"
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
      "Probability of collision between LEFs and extrusion barriers where barriers are pointing\n"
      "towards the extrusion direction.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  barr_adv.add_option(
      "--lef-bar-minor-collision-prob",
      c.lef_bar_minor_collision_pblock,
      "Probability of collision between LEFs and extrusion barriers where barriers are pointing\n"
      "away from the extrusion direction.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  barr_adv.add_option(
      "--extrusion-barrier-bound-stp",
      c.barrier_occupied_stp,
      "Self-transition probability for extrusion barriers in the \"bound\" state.\n"
      "In other words, the probability that an extrusion barrier that is active in the current epoch\n"
      "will remain active during the next epoch.")
      ->check(CLI::Range(0.0, 1.0));

  barr_adv.add_option(
      "--extrusion-barrier-not-bound-stp",
      c.barrier_not_occupied_stp,
      "Self-transition probability for extrusion barriers in the \"not-bound\" state.\n"
      "In other words, the probability that an extrusion barrier that is inactive in the current\n"
      "epoch will remain inactive during the next epoch.")
      ->check(CLI::Range(0.0, 1.0))
      ->capture_default_str();

  barr_adv.add_flag(
      "--interpret-extrusion-barrier-name-as-not-bound-stp",
      c.interpret_bed_name_field_as_barrier_not_occupied_stp,
      "Interpret the name field of the extrusion barrier BED file passed through --extrusion-barrier-file\n"
      "as the not-bound stp for the barrier defined by the BED record.")
      ->capture_default_str();

  cgen.add_option(
      "--contact-sampling-strategy",
      c.contact_sampling_strategy,
      fmt::format(FMT_STRING("Strategy to use when sampling contacts.\n"
                             "Should be one of:\n"
                             " - {}\n"
                             "When one of the *-with-noise strategies is specified, contacts are randomized by\n"
                             "applying a random offset to the location of LEF extrusion units.\n"
                             "Offsets are drawn from a genextreme distrubution.\n"
                             "The distribution parameters can be controlled through the options --mu, --sigma\n"
                             "and --xi."),
                  fmt::join(Cli::contact_sampling_strategy_map.keys_view(), "\n - ")))
      ->transform(CLI::CheckedTransformer(Cli::contact_sampling_strategy_map))
      ->capture_default_str();

  cgen.add_option(
      "--contact-sampling-interval",
      c.contact_sampling_interval,
      "Average number of base-pairs extruded by one LEF between two subsequent sampling events\n"
      "(assuming no collision have occurred).\n"
      "Specifying shorter intervals will increase the frequency of sampling frequency, producing more\n"
      "contacts when --stopping-criterion=simulation-epochs or reducing the number of epochs simulated\n"
      "in case --stopping-criterion=contact-density.\n"
      "Using sampling intervals that are significantly smaller than the average LEF processivity is usually\n"
      "not recommended, as will cause MoDLE to sample many contacts from a single loop.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit | utils::cli::AsGenomicDistance)
      ->capture_default_str();

  cgen.add_option(
      "-r,--resolution",
      c.bin_size,
      "Resolution in base-pairs for the output contact matrix in cooler format.\n"
      "NOTE: MoDLE simulation always take place at 1 bp resolution.\n"
      "      This parameter only affects the resolution of the output contact matrix.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit | utils::cli::AsGenomicDistance)
      ->capture_default_str();

  cgen.add_option(
      "-w,--diagonal-width",
      c.diagonal_width,
      "Width of the subdiagonal window for the in-memory contact matrix.\n"
      "This setting affects the maximum distance of a pair of bins whose interactions will be\n"
      "tracked by MoDLE.\n"
      "As a rule of thumb, --diagonal-width should roughly 10x the average LEF processivity.\n"
      "Setting --diagonal-width to very large values (i.e. tens of Mbp) will significantly\n"
      "increase MoDLE's memory requirements.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit | utils::cli::AsGenomicDistance)
      ->capture_default_str();

  cgen_adv.add_option(
      "--tad-to-loop-contact-ratio",
      c.tad_to_loop_contact_ratio,
      "Ratio between the number of TAD and loop contacts.\n"
      "Loop contacts are sampled using the position of the extrusion units of a LEF.\n"
      "This type of contact visually correspond to stripes and dots in the output matrix.\n"
      "TAD contacts are sampled from the body of a loop and visually correspond to TADs.\n"
      "Use \"0\" to disable sampling of TAD contacts.\n"
      "Use \"inf\" to disable sampling of loop contacts.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  cgen_adv.add_option(
      "--mu,--genextr-location",
      c.genextreme_mu,
      "Location parameter (mu) of the generalized extreme value distribution used to add noise to\n"
      "molecular contacts.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  cgen_adv.add_option(
      "--sigma,--genextr-scale",
      c.genextreme_sigma,
      "Scale parameter (sigma) of the generalized extreme value distribution used to add noise to\n"
      "molecular contacts.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  cgen_adv.add_option(
      "--xi,--genextr-shape",
      c.genextreme_xi,
      "Shape parameter (xi) of the generalized extreme value distribution used to add noise to\n"
      "molecular contacts.")
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
      "Each simulation instance will run exactly --target-number-of-epochs epochs after burn-in\n"
      "phase. Has no effect when --stopping-criterion=\"contact-density\".")
      ->check(CLI::PositiveNumber)
      ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit)
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
      "         Assuming chr1 is exactly 250 Mbp long, MoDLE will schedule simulation instances to yield\n"
      "         a total of 150 000 contacts for chr1: 2 * (250e6 / 100e3) * (3e6 / 100e3) = 150e3.")
      ->check(CLI::PositiveNumber);

  misc.add_option(
      "--ncells",
      c.num_cells,
      "Number of simulation instances or cells to simulate.\n"
      "Loop extrusion will be simulated independently for every chromosome across --ncells\n"
      "simulation instances. The final contact matrix is produced by accumulating contacts\n"
      "generated across all simulation instances. To achieve good performance and CPU utilization\n"
      "we recommend setting --ncells equal to the number of available CPU cores.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit)
      ->capture_default_str();

  misc.add_option(
      "-t,--threads",
      c.nthreads,
      "Number of worker threads used to run simulation instances.\n"
      "By default MoDLE will spawn a number of worker threads equal to the number of logical CPU\n"
      "cores. On high core count machines a slightly better performance can usually be obtained\n"
      "by setting --threads equal to the number of physical CPU cores.\n")
      ->check(CLI::PositiveNumber)
      ->transform(CLI::Bound(1U, std::thread::hardware_concurrency()))
      ->capture_default_str();

  misc.add_option(
      "--seed",
      c.seed,
      "Base seed to use for random number generation.")
      ->check(CLI::NonNegativeNumber)
      ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit)
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
      "By default this parameter is computed based on the average LEF processivity, overall\n"
      "extrusion speed and burn-in extrusion speed coefficient (--avg-lef-processivity,\n"
      "--fwd/rev-extrusion-speed and --burn-in-extr-speed-coefficient respectively).")
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
      "Window size used to smooth metrics monitored to decide when the burn-in phase should be\n"
      "terminated.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  burnin_adv.add_option(
      "--min-burnin-epochs",
      c.min_burnin_epochs,
      "Lower bound for the burn-in phase duration.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  burnin_adv.add_option(
      "--max-burnin-epochs",
      c.max_burnin_epochs,
      "Upper bound for the burn-in phase duration.\n"
      "This is especially useful when simlating loop extrusion with very few LEFs.\n"
      "When this is the case, --max-burnin-epochs=100000 can be used as a very conservative\n"
      "threshold.")
      ->check(CLI::PositiveNumber)
      ->default_str("inf");

  burnin_adv.add_option(
      "--burnin-extr-speed-coefficient",
      c.burnin_speed_coefficient,
      "Extrusion speed coefficient to apply during the burn-in phase.\n"
      "Setting this to numbers > 1.0 will speed-up the burn-in phase, as the average loop size will\n"
      "stabilize faster.\n"
      "IMPORTANT: Setting this parameter to values other than 1.0 comes with many gotchas.\n"
      "           For the time being, tuning this parameter is not reccommended.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  misc_adv.add_option(
      "--probability-normalization-factor",
      c.probability_normalization_factor,
      "Normalization factor used to scale probabilities based on the overall average extrusion speed.\n"
      "This is needed in order for the transition and collision probabilities to make sense when simulating\n"
      "using different extrusion speeds.\n"
      "Example: let us consider a simulation where:\n"
      "         - --extrusion-barrier-occupancy=0.5\n"
      "         - --extrusion-barrier-not-bound-stp=0.5\n"
      "         - --rev-extrusion-speed=5000\n"
      "         - --fwd-extrusion-speed=5000\n"
      "         In this scenario we expect extrusion barriers to flip between active and inactive state\n"
      "         every simulation epoch, that is every N bp of DNA extruded, where N is the product between\n"
      "         the overall extrusion speed (10kbp/epoch) and the number of LEFs that are being simulated in\n"
      "         one simulation instance.\n"
      "         If now we consider the same simulation carried out at much slower extrusion speed\n"
      "         (e.g. 100bp/epoch), the situation is quite different, as extrusion barriers switch between\n"
      "         active and inactive states once every N * (100 / 10 000) bp of DNA extruded.\n"
      "This parameter is used to normalize transition and collision probabilities to avoid the issue outlined\n"
      "in the example above. Basically specifying --probability-normalization-factor=10kb and\n"
      "--extrusion-barrier-not-bound-stp=0.8 instructs MoDLE to interpret the self-transition probability\n"
      "specified through the CLI at an extrusion speed of 10kbp, and adjust the transition probability\n"
      "based on the actual extrusion speed.\n"
      "Currently, changing this parameter affects the probabilities corresponding to the following CLI options:\n"
      " - --extrusion-barrier-bound-stp\n"
      " - --extrusion-barrier-not-bound-stp\n"
      " - --lef-bar-major-collision-prob\n"
      " - --lef-bar-minor-collision-prob\n"
      " - --probability-of-lef-bypass")
            ->check(CLI::PositiveNumber)
      ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit | utils::cli::AsGenomicDistance)
      ->capture_default_str();

  misc_adv.add_flag(
      "--normalize-probabilities,!--no-normalize-probabilities",
      c.normalize_probabilities,
      "Toggle on/off normalization of transition and collision probabilities using --probability-normalization-factor.")
      ->capture_default_str();

  // Address option dependencies/incompatibilities
  io_adv.get_option("--skip-output")->excludes(io_adv.get_option("--log-model-internal-state"));
  stopping.get_option("--target-contact-density")->excludes(stopping.get_option("--target-number-of-epochs"));
  lefbar.get_option("--extrusion-barrier-occupancy")->excludes(barr_adv.get_option("--extrusion-barrier-bound-stp"));
  barr_adv.get_option("--interpret-extrusion-barrier-name-as-not-bound-stp")->excludes(barr_adv.get_option("--extrusion-barrier-not-bound-stp"));
  // clang-format on

  // Deprecated options
  deprecated.add_option("--chrom-subranges", c.path_to_genomic_intervals)->check(CLI::ExistingFile);
  deprecated.get_option("--chrom-subranges")->excludes(io.get_option("--genomic-intervals"));

  std::array<std::reference_wrapper<CLI::App>, 10> option_groups{
      io, lefbar, cgen, stopping, misc, io_adv, lef_adv, barr_adv, cgen_adv, burnin_adv};
  std::vector<CLI::App*> option_groups_ptrs(option_groups.size());
  std::transform(option_groups.begin(), option_groups.end(), option_groups_ptrs.begin(),
                 [](auto& grp) { return &(grp.get()); });

  return option_groups_ptrs;
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
  auto option_group_ptrs = add_common_options(s, this->_config);
}

void Cli::make_cli() {
  this->_cli.description(
      "High-performance stochastic modeling of DNA loop extrusion interactions.");
  this->_cli.set_version_flag("-V,--version", std::string{modle::config::version::str_long()});
  this->_cli.require_subcommand(1);
  this->_cli.set_config("--config", "", "Path to MoDLE's config file (optional).", false);
  this->_cli.formatter(std::make_shared<utils::cli::Formatter>());
  this->_cli.get_formatter()->column_width(30);

  this->make_simulation_subcommand();
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
  }

  this->validate_args();
  this->transform_args();
  this->_config.args = absl::MakeSpan(_argv, static_cast<usize>(_argc));

  this->_config.args_json = this->to_json();
  return this->_config;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity) TODO: reduce complexity
std::string Cli::detect_path_collisions(modle::Config& c) const {
  std::string collisions;

  if (c.force || c.skip_output) {
    return "";
  }

  auto check_for_path_collisions = [](const std::filesystem::path& path) -> std::string {
    assert(!path.empty());
    if (std::filesystem::exists(path)) {
      if (std::filesystem::is_directory(path)) {
        return fmt::format(
            FMT_STRING("Refusing to continue because output file {} already "
                       "exist (and is actually a {}directory).\n{}.\n"),
            path, std::filesystem::is_empty(path) ? "" : "non-empty ",
            std::filesystem::is_empty(path)
                ? " Pass --force to overwrite"
                : "You should specify a different output path, or manually remove the "
                  "existing directory");
      }
      return fmt::format(FMT_STRING("Refusing to continue because output file {} already exist.\n"
                                    "Pass --force to overwrite.\n"),
                         path);
    }
    if (std::filesystem::is_directory(path) && !std::filesystem::is_empty(path)) {
      return fmt::format(
          FMT_STRING("Refusing to continue because output file {} is a non-empty directory.\n"
                     "You should specify a different output path, or "
                     "manually remove the existing directory.\n"),
          path);
    }
    return {};
  };

  if (std::filesystem::exists(c.path_to_log_file)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_log_file));
  }

  if (c.track_1d_lef_position && std::filesystem::exists(c.path_to_lef_1d_occupancy_bw_file)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_lef_1d_occupancy_bw_file));
  }

  if (this->get_subcommand() == simulate && !c.path_to_model_state_log_file.empty() &&
      std::filesystem::exists(c.path_to_model_state_log_file)) {
    absl::StrAppend(&collisions, check_for_path_collisions(c.path_to_model_state_log_file));
  }

  return collisions;
}

int Cli::exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }

static void detect_deprecated_option_usage(const CLI::App& subcmd,
                                           std::vector<std::string>& warnings_buff) {
  const auto deprecated_opt_grp = subcmd.get_option_group("");
  if (!deprecated_opt_grp->parsed()) {
    return;
  }

  // Map deprecated options to their replacement
  static constexpr utils::ConstMap<std::string_view, std::string_view, 1>
      deprecated_option_mappings{{"--chrom-subranges", "--genomic-intervals"}};

  for (const auto& opt : deprecated_opt_grp->parse_order()) {
    const auto match = deprecated_option_mappings.find(opt->get_name());
    if (match != deprecated_option_mappings.end()) {
      warnings_buff.emplace_back(fmt::format(FMT_STRING("Option {} is deprecated. Use {} instead."),
                                             *match.first, *match.second));
    }
  }
}

void Cli::validate_args() const {
  const auto& c = this->_config;
  std::vector<std::string> errors;

  const auto* subcmd = this->get_subcommand_ptr();

  detect_deprecated_option_usage(*subcmd, this->_warnings);

  if (c.burnin_smoothing_window_size > c.burnin_history_length) {
    assert(subcmd->get_option_group("Advanced")
               ->get_option_group("Burn-in")
               ->get_option("--burnin-smoothing-window-size"));
    assert(subcmd->get_option_group("Advanced")
               ->get_option_group("Burn-in")
               ->get_option("--burnin-history-length"));
    errors.emplace_back(fmt::format(
        FMT_STRING("The value passed to {} should be less or equal than that of {} ({} > {})"),
        "--burnin-smoothing-window-size", "--burnin-history-length", c.burnin_smoothing_window_size,
        c.burnin_history_length));
  }

  const auto& cgen_adv =
      subcmd->get_option_group("Advanced")->get_option_group("Contact generation");

  using CS = Config::ContactSamplingStrategy;
  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  if (!(c.contact_sampling_strategy & CS::noisify)) {
    for (const std::string label : {"--mu", "--sigma", "--xi"}) {
      if (!cgen_adv->get_option(label)->empty()) {
        this->_warnings.emplace_back(fmt::format(
            FMT_STRING("Option {} has no effect. Reason: {} requires the strategy passed to "
                       "--contact-sampling-strategy to be one of the *-with-noise strategies."),
            label, label));
      }
    }
  }

  if (!cgen_adv->get_option("--tad-to-loop-contact-ratio")->empty()) {
    const bool sample_loop_contacts = c.contact_sampling_strategy & CS::loop;
    const bool sample_tad_contacts = c.contact_sampling_strategy & CS::tad;
    assert(sample_loop_contacts || sample_tad_contacts);

    if (sample_loop_contacts && !sample_tad_contacts && c.tad_to_loop_contact_ratio != 0) {
      this->_warnings.emplace_back(
          fmt::format(FMT_STRING("Option --tad-to-loop-contact-ratio={} has no effect. Reason: "
                                 "--tad-to-loop-contact-ratio is implicitly set to 0 when "
                                 "--contact-sampling-strategy={}"),
                      c.tad_to_loop_contact_ratio,
                      Cli::contact_sampling_strategy_map.at(c.contact_sampling_strategy)));
    }
    if (!sample_loop_contacts && sample_tad_contacts &&
        c.tad_to_loop_contact_ratio != std::numeric_limits<double>::infinity()) {
      this->_warnings.emplace_back(
          fmt::format(FMT_STRING("Option --tad-to-loop-contact-ratio={} has no effect. Reason: "
                                 "--tad-to-loop-contact-ratio is implicitly set to inf when "
                                 "--contact-sampling-strategy={}"),
                      c.tad_to_loop_contact_ratio,
                      Cli::contact_sampling_strategy_map.at(c.contact_sampling_strategy)));
    }
  }

  if (c.stopping_criterion == Config::StoppingCriterion::simulation_epochs &&
      subcmd->get_option_group("Stopping criterion")
          ->get_option("--target-number-of-epochs")
          ->empty()) {
    errors.emplace_back(
        "--stopping-criterion=simulation-epochs requires option --target-number-of-epochs to be "
        "specified and to be finite.");
  }

  if (c.stopping_criterion == Config::StoppingCriterion::contact_density &&
      !subcmd->get_option_group("Stopping criterion")
           ->get_option("--target-number-of-epochs")
           ->empty()) {
    errors.emplace_back(
        "--stopping-criterion=contact-density excludes option --target-number-of-epochs.");
  }

  if (!c.normalize_probabilities && !subcmd->get_option_group("Advanced")
                                         ->get_option_group("Miscellaneous")
                                         ->get_option("--probability-normalization-factor")
                                         ->empty()) {
    this->_warnings.emplace_back(
        fmt::format(FMT_STRING("Option --probability-normalization-factor has no effect. Reason: "
                               "CLI option --no-normalize-probabilities was passed by the user.")));
  }

  if (c.min_burnin_epochs > c.max_burnin_epochs) {
    errors.emplace_back(fmt::format(
        FMT_STRING("--min-burnin-epochs={} cannot be greater than --max-burnin-epochs={}."),
        c.min_burnin_epochs, c.max_burnin_epochs));
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

// The following two functions have been copied from the ExtrusionBarrier class.
// This is not ideal, but I think it is better than make the CLI interface depend on the
// ExtrusionBarrier class (which is part of libmodle_internal)
constexpr double compute_stp_active_from_occupancy(double stp_inactive, double occupancy) noexcept {
  if (MODLE_UNLIKELY(occupancy == 0)) {
    return 0.0;
  }

  const auto tp_inactive_to_active = 1.0 - stp_inactive;
  const auto tp_active_to_inactive =
      (tp_inactive_to_active - (occupancy * tp_inactive_to_active)) / occupancy;
  return std::clamp(1.0 - tp_active_to_inactive, 0.0, 1.0);
}

constexpr double compute_occupancy_from_stp(double stp_active, double stp_inactive) noexcept {
  if (MODLE_UNLIKELY(stp_active + stp_inactive == 0)) {
    return 0.0;
  }

  const auto tp_inactive_to_active = 1.0 - stp_inactive;
  const auto tp_active_to_inactive = 1.0 - stp_active;
  const auto occupancy = tp_inactive_to_active / (tp_inactive_to_active + tp_active_to_inactive);
  return std::clamp(occupancy, 0.0, 1.0);
}

/// Generate output file paths from output prefix
static void cli_update_paths(Cli::subcommand subcommand, Config& c) {
  c.path_to_output_file_cool = c.path_to_output_prefix;
  c.path_to_log_file = c.path_to_output_prefix;
  c.path_to_config_file = c.path_to_output_prefix;
  if (subcommand == Cli::subcommand::simulate) {
    c.path_to_model_state_log_file = c.path_to_output_prefix;
    c.path_to_model_state_log_file += "_internal_state.log.zst";
  }

  if (c.track_1d_lef_position) {
    c.path_to_lef_1d_occupancy_bw_file = c.path_to_output_prefix;
    c.path_to_lef_1d_occupancy_bw_file += "_lef_1d_occupancy.bw";
  }

  c.path_to_output_file_cool += ".cool";
  c.path_to_log_file += ".log";
  c.path_to_config_file += "_config.toml";
}

/// Compute mean and std for rev and fwd extrusion speed
static void cli_update_extr_speed(const CLI::App& cli, Config& c) {
  auto* grp = cli.get_subcommand("simulate")
                  ->get_option_group("Advanced")
                  ->get_option_group("Extrusion Factors");
  const auto rev_speed_parsed = !grp->get_option("--rev-extrusion-speed")->empty();
  const auto fwd_speed_parsed = !grp->get_option("--fwd-extrusion-speed")->empty();

  if (auto& speed = c.rev_extrusion_speed; !rev_speed_parsed) {
    speed = c.bin_size * 8 / 10;  // 0.8 * c.bin_size
  }
  if (auto& speed = c.fwd_extrusion_speed; !fwd_speed_parsed) {
    speed = c.bin_size * 8 / 10;
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
}

/// Compute the probability of LEF release from the avg. processivity and overall extr. speed
static constexpr void cli_compute_prob_of_lef_release(Config& c) {
  c.prob_of_lef_release = static_cast<double>(c.rev_extrusion_speed + c.fwd_extrusion_speed) /
                          static_cast<double>(c.avg_lef_processivity);
  c.prob_of_lef_release_burnin =
      static_cast<double>(c.rev_extrusion_speed_burnin + c.fwd_extrusion_speed_burnin) /
      static_cast<double>(c.avg_lef_processivity);
}

/// Compute the self transition probabilities used to model extrusion barriers
static void cli_update_barrier_stp_and_occupancy(const CLI::App& cli, Config& c) {
  auto* grp = cli.get_subcommand("simulate")->get_option_group("Extrusion Barriers and Factors");

  const auto extrusion_barrier_occupancy_parsed =
      !grp->get_option("--extrusion-barrier-occupancy")->empty();

  if (extrusion_barrier_occupancy_parsed) {
    c.barrier_occupied_stp = compute_stp_active_from_occupancy(c.barrier_not_occupied_stp,
                                                               c.extrusion_barrier_occupancy);
  } else {
    c.extrusion_barrier_occupancy =
        compute_occupancy_from_stp(c.barrier_occupied_stp, c.barrier_not_occupied_stp);
  }
}

/// Normalize transition probabilities based on the overall extrusion speed
static void cli_normalize_probabilities(Config& c) {
  if (const auto ratio = static_cast<double>(c.rev_extrusion_speed + c.fwd_extrusion_speed) /
                         static_cast<double>(c.probability_normalization_factor);
      ratio != 1.0) {
    /// Numerically stable alternative to std::pow(base, exp) for positive bases
    auto stable_pow = [](double base, double exp) noexcept {
      assert(base >= 0);
      if (base == 0.0) {
        return 0.0;
      }
      if (base == 1.0) {
        return 1.0;
      }
      return std::exp(std::log(base) * exp);
    };

    // It is important that we recompute the barrier_occupued_stp after correcting
    // barrier_not_occupied_stp!
    c.barrier_not_occupied_stp = stable_pow(c.barrier_not_occupied_stp, ratio);
    c.barrier_occupied_stp = compute_stp_active_from_occupancy(c.barrier_not_occupied_stp,
                                                               c.extrusion_barrier_occupancy);

    if (const auto p = c.probability_of_extrusion_unit_bypass; p != 0.0 && p != 1.0) {
      c.probability_of_extrusion_unit_bypass = std::min(p * ratio, 1.0);
    }

    c.lef_bar_major_collision_pblock = stable_pow(c.lef_bar_major_collision_pblock, ratio);
    c.lef_bar_minor_collision_pblock = stable_pow(c.lef_bar_minor_collision_pblock, ratio);
  }
}

static constexpr void cli_update_tad_to_loop_contact_ratio(Config& c) {
  using CS = Config::ContactSamplingStrategy;

  const auto sample_loop_contacts = bool(c.contact_sampling_strategy & CS::loop);
  const auto sample_tad_contacts = bool(c.contact_sampling_strategy & CS::tad);

  assert(sample_loop_contacts || sample_tad_contacts);
  if (sample_loop_contacts && !sample_tad_contacts) {
    c.tad_to_loop_contact_ratio = 0;
  }
  if (!sample_loop_contacts && sample_tad_contacts) {
    c.tad_to_loop_contact_ratio = std::numeric_limits<double>::infinity();
  }
}

void cli_update_burnin_params(Config& c) {
  const auto lef_activation_bp = 5 * c.avg_lef_processivity;
  c.burnin_target_epochs_for_lef_activation = std::min(
      c.max_burnin_epochs,
      utils::conditional_static_cast<usize>(
          lef_activation_bp / (c.rev_extrusion_speed_burnin + c.fwd_extrusion_speed_burnin)));
}

void Cli::transform_args() {
  cli_update_paths(this->get_subcommand(), this->_config);
  cli_update_extr_speed(this->_cli, this->_config);
  cli_compute_prob_of_lef_release(this->_config);
  cli_update_barrier_stp_and_occupancy(this->_cli, this->_config);
  cli_update_tad_to_loop_contact_ratio(this->_config);
  cli_update_burnin_params(this->_config);

  if (this->_config.normalize_probabilities) {
    cli_normalize_probabilities(this->_config);
  }

  const auto* subcmd = this->get_subcommand_ptr();

  if (!subcmd->get_option_group("Extrusion Barriers and Factors")
           ->get_option("--extrusion-barrier-occupancy")
           ->empty()) {
    this->_config.override_extrusion_barrier_occupancy = true;
  }

  if (this->_config.stopping_criterion == Config::StoppingCriterion::simulation_epochs) {
    this->_config.target_contact_density = -1;
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

void Cli::log_warnings() const {
  for (const auto& warning : this->_warnings) {
    spdlog::warn(FMT_STRING("{}"), warning);
  }
}

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
    // Given two subcommands named comm1 and comm2, assuming comm1 was parsed while comm2 was
    // not, the TOML produced by CLI11 will have values for comm2 formatted as comm2.myarg1=1,
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
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while converting MoDLE's config "
                               "from TOML to JSON: {}"),
                    e.what()));
  }
}

}  // namespace modle
