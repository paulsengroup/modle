#pragma once

#include "CLI/CLI.hpp"

namespace modle {
struct config {
  std::string_view path_to_bed;
  std::string_view output_dir;

  uint32_t bin_size{1'000};
  uint32_t simulation_iterations{5'000'000};
  uint32_t number_of_barriers;
  uint32_t number_of_lefs;
  uint32_t average_lef_processivity{100'000};
  double probability_of_barrier_block{0.95};
  double probability_of_lef_rebind{1.0};
  double probability_of_extrusion_unit_bypass{0.25};
  bool make_heatmaps{true};
  uint64_t seed{0};

  void print() const {
    // clang-format off
    absl::FPrintF(
        stderr,
        "###### CONFIG SUMMARY ######\n"
        "### Input/Output##\n"
        "##    Input BED:                        '%s'\n"
        "##    Output directory:                 '%s'\n"
        "### General settings\n"
        "##    Bin size (bp):                    %lu\n"
        "##    # of iterations:                  %lu\n"
        "##    Avg. LEF processivity (bp):       %lu\n"
        "##    # of randomly generated barriers: %lu\n"
        "##    # of randomly generated LEFs:     %lu\n"
        "### Probabilities:\n"
        "##    Prob. of barrier block:           %.4g\n"
        "##    Prob. of LEF rebind:              %.4g\n"
        "##    Prob. of LEF bypass:              %.4g\n"
        "### Various:\n"
        "##    Generate heatmaps:                %s\n"
        "##    Seed:                             %lu"
        "\n\n",
        this->path_to_bed, this->output_dir, this->bin_size, this->simulation_iterations,
        this->average_lef_processivity, this->number_of_barriers, this->number_of_lefs,
        this->probability_of_barrier_block, this->probability_of_lef_rebind,
        this->probability_of_extrusion_unit_bypass, this->make_heatmaps ? "Yes" : "False",
        this->seed);
    // clang-format on
  }
};

class Cli {
 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  config _config;
  CLI::App _cli{};

  void MakeCli() {
    // clang-format off
      this->_cli.description("in-silico modeling of DNA loop extrusion.");

      this->_cli.add_option("-i,--input-bed", this->_config.path_to_bed, "Path to input BED file.")->check(CLI::ExistingFile)->required();
      this->_cli.add_option("-o,--output-dir", this->_config.output_dir, "Path to output directory. It will be created if it doesn't already exist.")->required();
      this->_cli.add_option("-b,--bin-size", this->_config.bin_size, "Bin size in base pairs.")->check(CLI::PositiveNumber)->capture_default_str();
      this->_cli.add_option("--number-of-iterations", this->_config.simulation_iterations, "Number of simulation iterations.")->check(CLI::PositiveNumber)->capture_default_str();
      // TODO: remove requited after writing the code to read barriers and LEFs coordinates from a BED file.
      this->_cli.add_option("--number-of-randomly-generated-barriers", this->_config.number_of_barriers, "Number of extrusion barriers to be randomly generated and bound.")->check(CLI::PositiveNumber)->required();
      this->_cli.add_option("--number-of-randomly-generated-lefs", this->_config.number_of_lefs, "Number of loop extrusion factors (LEFs) to be randomly generated and bound.")->check(CLI::PositiveNumber)->required();
      this->_cli.add_option("--avg-lef-processivity", this->_config.average_lef_processivity, "Average loop extrusion factor processivity, or in other words, average loop size in base pairs.")->check(CLI::PositiveNumber)->capture_default_str();
      this->_cli.add_option("--probability-of-barrier-block", this->_config.probability_of_barrier_block, "Probability of extrusion block by an extrusion barrier.")->check(CLI::Range(0.0, 1.0))->capture_default_str();
      this->_cli.add_option("--probability-of-lef-rebind", this->_config.probability_of_lef_rebind, "Probability that an unbound loop extruding factor (LEF) will randomly rebind to DNA.")->check(CLI::Range(0.0, 1.0))->capture_default_str();
      this->_cli.add_option("--probability-of-lef-bypass", this->_config.probability_of_extrusion_unit_bypass, "Probability that a loop extruding factor (LEF) will not block when meeting another LEF.")->check(CLI::Range(0.0, 1.0))->capture_default_str();
      this->_cli.add_flag("--make-heatmaps,!--no-make-heatmaps", this->_config.make_heatmaps, "Generate heatmaps from the complete and diagonal contact matrices.")->capture_default_str();
      this->_cli.add_option("--seed", this->_config.seed, "Seed to use for random number generation.")->check(CLI::NonNegativeNumber)->capture_default_str();
    // clang-format on
  }

 public:
  Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(argv[0]) { this->MakeCli(); }
  [[nodiscard]] config parse_arguments() {
    this->_cli.name(this->_exec_name);
    this->_cli.parse(this->_argc, this->_argv);
    return this->_config;
  }
  [[nodiscard]] int exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }
};

}  // namespace modle
