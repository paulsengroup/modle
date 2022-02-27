// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <CLI/CLI.hpp>  // for App
#include <string>       // for string

#include "./cli.hpp"
#include "modle/common/cli_utils.hpp"
#include "modle/common/common.hpp"             // for u8
#include "modle/common/simulation_config.hpp"  // for Config
namespace CLI {
class ParseError;
}  // namespace CLI

namespace CLI {
class ParseError;
}  // namespace CLI

namespace modle {

class Cli {
 public:
  enum subcommand : u8f { help, simulate, pertubate, replay };

  Cli(int argc, char** argv);
  [[nodiscard]] const Config& parse_arguments();
  [[nodiscard]] std::string detect_path_collisions(modle::Config& c) const;

  [[nodiscard]] int exit(const CLI::ParseError& e) const;
  [[nodiscard]] subcommand get_subcommand() const;
  void print_config(bool print_default_args = false) const;
  void write_config_file(bool write_default_args = false) const;
  [[nodiscard]] std::string to_json() const;

  [[nodiscard]] bool config_file_parsed() const;

 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  Config _config;
  CLI::App _cli{};
  subcommand _subcommand{subcommand::help};

  void make_cli();
  void make_simulation_subcommand();
  void make_perturbate_subcommand();
  void make_replay_subcommand();
  void validate_args() const;
  void transform_args();
  [[nodiscard]] CLI::App* get_subcommand_ptr();
  [[nodiscard]] const CLI::App* get_subcommand_ptr() const;

 public:
  using StoppingCriterionMappings = utils::CliEnumMappings<Config::StoppingCriterion>;
  inline static const StoppingCriterionMappings stopping_criterion_map{
      std::make_pair("contact-density", Config::StoppingCriterion::contact_density),
      std::make_pair("simulation-epochs", Config::StoppingCriterion::simulation_epochs)};

  using CS_ = Config::ContactSamplingStrategy;
  using CS_ut_ = CS_::underlying_type;
  // It is important that we use the underlying type in this map
  using ContactSamplingMethodMappings = utils::CliEnumMappings<CS_ut_>;
  inline static const ContactSamplingMethodMappings contact_sampling_strategy_map{
      // clang-format off
      std::make_pair("tad-only",                 CS_ut_(CS_::tad)),
      std::make_pair("loop-only",                CS_ut_(CS_::loop)),
      std::make_pair("tad-plus-loop",            CS_ut_(CS_::tad  | CS_::loop)),
      std::make_pair("tad-only-with-noise",      CS_ut_(CS_::tad  | CS_::noisify)),
      std::make_pair("loop-only-with-noise",     CS_ut_(CS_::loop | CS_::noisify)),
      std::make_pair("tad-plus-loop-with-noise", CS_ut_(CS_::tad  | CS_::loop    | CS_::noisify))
      // clang-format on
  };
};

}  // namespace modle
