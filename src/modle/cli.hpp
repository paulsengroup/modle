#pragma once

#include <CLI/CLI.hpp>

#include "modle/config.hpp"

namespace modle {

class Cli {
 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  config _config;
  CLI::App _cli{};

  inline void MakeCli() {
    this->_cli.description("in-silico modeling of DNA loop extrusion.");
    auto* io = this->_cli.add_option_group("Input/Output", "");
    auto* gen = this->_cli.add_option_group("Generic", "");
    auto* prob = this->_cli.add_option_group("Probabilities", "");
    auto* rand = this->_cli.add_option_group("Random", "");
    auto* burn = this->_cli.add_option_group("Burn-in", "");
    auto* hidden = this->_cli.add_option_group("", "");
    auto* extr_barr = this->_cli.add_option_group("Extrusion Barriers", "");
    auto* extr_barr_mandatory = extr_barr->add_option_group("Mandatory", "")->require_option(1);

    auto remove_trailing_zeros_from_floats = [](const std::string& s) {
      return std::string{absl::StripSuffix(s, ".0")};
    };

    // clang-format off
    io->add_option(
            "-c,--chromosome-size-file",
            this->_config.path_to_chr_sizes,
            "Path to file with chromosome size (should be in a TSV-like format with two fields per line).")
            ->check(CLI::ExistingFile)->required();

    io->add_option(
            "--chromosome-subrange-file",
            this->_config.path_to_chr_subranges,
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

    gen->add_option(
            "-b,--bin-size",
            this->_config.bin_size,
            "Bin size in base pairs.")
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
            "--number-of-randomly-generated-lefs",
            this->_config.number_of_lefs,
            "Number of loop extrusion factors (LEFs) to be randomly generated and bound.")
            ->check(CLI::NonNegativeNumber)
            ->transform((remove_trailing_zeros_from_floats))
            ->required();

    gen->add_option(
            "--avg-lef-processivity",
            this->_config.average_lef_processivity,
            "Average loop extrusion factor processivity, or in other words, average loop size in base pairs.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zeros_from_floats))
            ->capture_default_str();

    gen->add_option(
            "--lef-unloader-strength",
            this->_config.lef_unloader_strength,
            "Coefficient to control the increase in the stability of the LEF-DNA bind when a LEF is stalled in both direction by two extrusion barriers in convergent orientation.")
            ->check(CLI::Range(0.0, 1.0))
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
            "--probability-of-lef-rebind",
            this->_config.probability_of_lef_rebind,
            "Probability that an unbound loop extruding factor (LEF) will randomly rebind to DNA.")
            ->check(CLI::Range(0.0, 1.0))
            ->capture_default_str();

    prob->add_option(
            "--probability-of-lef-bypass",
            this->_config.probability_of_extrusion_unit_bypass,
            "Probability that a loop extruding factor (LEF) will not block when meeting another LEF.")
            ->check(CLI::Range(0.0, 1.0))
            ->capture_default_str();

    rand->add_option(
            "--contact-sampling-interval",
            this->_config.contact_sampling_interval,
            "Number of simulation rounds between contact sampling. When specifying this is specified, contact sampling will be performed at constant intervals.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zeros_from_floats))
            ->capture_default_str();

    rand->add_flag(
            "--randomize-contact-sampling-interval,!--uniform-contact-sampling-interval",
            this->_config.randomize_contact_sampling_interval,
            "Sample contacts using a bernoulli process with a probability of success of 1 / n, where n is the interval set through --contact-sampling-interval.")
            ->capture_default_str();

    burn->add_option(
            "--min-burnin-rounds",
            this->_config.min_n_of_burnin_rounds,
            "Minimum number of extrusion rounds to simulate during the burn-in phase (set to 0 to disable).")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zeros_from_floats))
            ->capture_default_str();

    burn->add_option(
            "--min-number-of-loops-per-lef",
            this->_config.min_n_of_loops_per_lef,
            "Minimum number of loops that each LEF must have produced before stopping the burn-in phase.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zeros_from_floats))
            ->capture_default_str();

    extr_barr_mandatory->add_option(
            "--extrusion-barrier-file",
            this->_config.path_to_extr_barriers_bed,
            "Path to BED file containing the extrusion barriers to be used in the simulation. Barriers corresponding to chromosomes that are not being simulated will be silently discarded.")
            ->check(CLI::ExistingFile);

    extr_barr_mandatory->add_option(
            "--number-of-randomly-generated-barriers",
            this->_config.number_of_randomly_gen_extr_barriers,
            "Number of extrusion barriers to be randomly generated and bound.")
            ->check(CLI::PositiveNumber)
            ->transform((remove_trailing_zeros_from_floats))
            ->capture_default_str();

    extr_barr->add_option(
            "--probability-of-barrier-block",
            this->_config.probability_of_extrusion_barrier_block,
            "Probability of extrusion block by an extrusion barrier. When --extrusion-barrier-file is passed, this setting will overwrite the probability of block specified in the BED file.")
            ->check(CLI::Range(0.0, 1.0));

    extr_barr->add_flag(
            "--exclude-chr-wo-barriers,!--keep-chr-without-barriers",
            this->_config.exclude_chr_wo_extr_barriers,
            fmt::format(FMT_STRING("Do not simulate loop extrusion on chromosomes without any extrusion barrier. Default: {}"), this->_config.exclude_chr_wo_extr_barriers ? "--exclude-chr-wo-barriers" : "--keep-chr-without-barriers"))
            ->capture_default_str();

    hidden->add_flag("--skip-output", this->_config.skip_output, "Don't write output files and plots. Useful for profiling");

    gen->get_option("--skip-burn-in")->excludes(burn->get_option("--min-burnin-rounds"));
    gen->get_option("--skip-burn-in")->excludes(burn->get_option("--min-number-of-loops-per-lef"));
    this->_cli.get_option("--extrusion-barrier-file")->excludes(this->_cli.get_option("--number-of-randomly-generated-barriers"));
    // clang-format on
  }

 public:
  inline Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) {
    this->MakeCli();
  }
  [[nodiscard]] const config& parse_arguments() {
    this->_cli.name(this->_exec_name);
    this->_cli.parse(this->_argc, this->_argv);
    this->_config.argc = _argc;
    this->_config.argv = _argv;
    return this->_config;
  }
  [[nodiscard]] inline int exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }
};

}  // namespace modle
