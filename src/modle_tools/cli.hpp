// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <CLI/CLI.hpp>  // for App
#include <string>       // for string
#include <string_view>  // for string_view

#include "modle_tools/modle_tools_config.hpp"  // for config

namespace modle::tools {

class Cli {
 public:
  enum subcommand : u8f {
    help,
    eval,
    fbcl,
    noisify,
    transform,
  };
  Cli(int argc, char** argv);
  [[nodiscard]] bool is_ok() const noexcept;
  [[nodiscard]] subcommand get_subcommand() const noexcept;
  [[nodiscard]] std::string_view get_printable_subcommand() const noexcept;
  [[nodiscard]] modle_tools_config parse_arguments();
  [[nodiscard]] int exit(const CLI::ParseError& e) const;
  [[nodiscard]] std::string to_json() const;
  [[nodiscard]] static std::string_view subcommand_to_str(subcommand s) noexcept;

 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  int _exit_code{1};
  modle_tools_config _config{};
  CLI::App _cli{};
  subcommand _subcommand{subcommand::help};

  void make_eval_subcommand();
  void make_find_barrier_clusters_subcommand();
  void make_noisify_subcommand();
  void make_transform_subcommand();
  void make_cli();

  void validate_eval_subcommand() const;
  void validate_find_barrier_clusters_subcommand() const;
  void validate_noisify_subcommand() const;
  void validate_transform_subcommand() const;
  void validate() const;
};

}  // namespace modle::tools
