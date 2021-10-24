// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <CLI/CLI.hpp>  // for App
#include <string>       // for string

#include "./cli.hpp"
#include "modle/common/common.hpp"  // for u8
#include "modle/common/config.hpp"  // for Config
namespace CLI {
class ParseError;
}  // namespace CLI

namespace CLI {
class ParseError;
}  // namespace CLI

namespace modle {

class Cli {
 public:
  enum subcommand : u8 { help, simulate, pertubate, replay };

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
};

}  // namespace modle
