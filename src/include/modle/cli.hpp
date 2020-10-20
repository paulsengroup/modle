#pragma once

#include "CLI/CLI.hpp"
#include "modle/config.hpp"

namespace modle {

class Cli {
 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  config _config;
  CLI::App _cli{};

  void MakeCli() {
    this->_cli.description("in-silico modeling of DNA loop extrusion.");
    auto* io = this->_cli.add_option_group("Input/Output", "");
    auto* gen = this->_cli.add_option_group("Generic", "");
    auto* prob = this->_cli.add_option_group("Probabilities", "");
    auto* rand = this->_cli.add_option_group("Random", "");
    auto* burn = this->_cli.add_option_group("Burn-in", "");
    auto* hidden = this->_cli.add_option_group("", "");

    // clang-format off
    io->add_option("-c,--chromosome-size-file", this->_config.path_to_chr_size, "Path to file with chromosome size (should be in a TSV-like format with two fields per line).")->check(CLI::ExistingFile)->required();
    io->add_option("-o,--output-dir", this->_config.output_dir, "Path to output directory. It will be created if it doesn't already exist.")->required();
    io->add_flag("--make-heatmaps,!--no-make-heatmaps", this->_config.make_heatmaps, "Generate heatmaps from the complete and diagonal contact matrices.")->capture_default_str();
    io->add_flag("--force", this->_config.force, "Overwrite existing output files.")->capture_default_str();
    gen->add_option("--number-of-iterations", this->_config.simulation_iterations, "Number of simulation iterations.")->check(CLI::PositiveNumber)->capture_default_str();
    gen->add_option("--avg-lef-processivity", this->_config.average_lef_processivity, "Average loop extrusion factor processivity, or in other words, average loop size in base pairs.")->check(CLI::PositiveNumber)->capture_default_str();
    gen->add_option("--lef-unloader-strength", this->_config.lef_unloader_strength, "Coefficient to control the increase in the stability of the LEF-DNA bind when a LEF is stalled in both direction by two extrusion barriers in convergent orientation.")->check(CLI::Range(0.0, 1.0))->capture_default_str();
    gen->add_flag("--skip-burn-in", this->_config.skip_burnin, "Skip burn-in phase and start counting contacts from the first extrusion round.")->capture_default_str();
    gen->add_option("--seed", this->_config.seed, "Seed to use for random number generation.")->check(CLI::NonNegativeNumber)->capture_default_str();
    prob->add_option("--probability-of-barrier-block", this->_config.probability_of_barrier_block, "Probability of extrusion block by an extrusion barrier.")->check(CLI::Range(0.0, 1.0))->capture_default_str();
    prob->add_option("--probability-of-lef-rebind", this->_config.probability_of_lef_rebind, "Probability that an unbound loop extruding factor (LEF) will randomly rebind to DNA.")->check(CLI::Range(0.0, 1.0))->capture_default_str();
    prob->add_option("--probability-of-lef-bypass", this->_config.probability_of_extrusion_unit_bypass, "Probability that a loop extruding factor (LEF) will not block when meeting another LEF.")->check(CLI::Range(0.0, 1.0))->capture_default_str();
    // TODO: remove requited after writing the code to read barriers and LEFs coordinates from a BED file.
    rand->add_option("-b,--bin-size", this->_config.bin_size, "Bin size in base pairs.")->check(CLI::PositiveNumber)->capture_default_str();
    rand->add_option("--number-of-randomly-generated-barriers", this->_config.number_of_barriers, "Number of extrusion barriers to be randomly generated and bound.")->check(CLI::PositiveNumber)->required();
    rand->add_option("--number-of-randomly-generated-lefs", this->_config.number_of_lefs, "Number of loop extrusion factors (LEFs) to be randomly generated and bound.")->check(CLI::PositiveNumber)->required();
    rand->add_option("--contact-sampling-interval", this->_config.contact_sampling_interval, "Number of simulation rounds between contact sampling. When specifying this is specified, contact sampling will be performed at constant intervals.")->check(CLI::Range(1U, UINT32_MAX))->capture_default_str();
    rand->add_flag("--randomize-contact-sampling-interval,!--uniform-contact-sampling-interval", this->_config.randomize_contact_sampling, "Sample contacts using a bernoulli process with a probability of success of 1 / n, where n is the interval set through --contact-sampling-interval.")->capture_default_str();
    burn->add_option("--min-burnin-rounds", this->_config.min_n_of_burnin_rounds, "Minimum number of extrusion rounds to simulate during the burn-in phase (set to 0 to disable).")->check(CLI::PositiveNumber)->capture_default_str();
    burn->add_option("--min-number-of-loops-per-lef", this->_config.min_n_of_loops_per_lef, "Minimum number of loops that each LEF must have produced before stopping the burn-in phase.")->check(CLI::PositiveNumber)->capture_default_str();

    gen->get_option("--skip-burn-in")->excludes(burn->get_option("--min-burnin-rounds"));
    gen->get_option("--skip-burn-in")->excludes(burn->get_option("--min-number-of-loops-per-lef"));

    hidden->add_flag("--skip-output", this->_config.skip_output, "Don't write output files and plots. Useful for profiling");
    // clang-format on
  }

 public:
  Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(argv[0]) { this->MakeCli(); }
  [[nodiscard]] config parse_arguments() {
    this->_cli.name(this->_exec_name);
    this->_cli.parse(this->_argc, this->_argv);
    this->_config.argc = _argc;
    this->_config.argv = _argv;
    return this->_config;
  }
  [[nodiscard]] int exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }
};

}  // namespace modle
