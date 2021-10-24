// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <CLI/CLI.hpp>  // for App
#include <cstdint>      // for uint_fast8_t
#include <string>       // for string
#include <string_view>  // for string_view

#include "modle_tools/config.hpp"  // for config

namespace modle::tools {

class Cli {
 public:
  enum subcommand : std::uint_fast8_t {
    help,
    eval,
    fbcl,
    noisify,
    stats,
  };
  Cli(int argc, char** argv);
  [[nodiscard]] bool is_ok() const noexcept;
  [[nodiscard]] subcommand get_subcommand() const noexcept;
  [[nodiscard]] std::string_view get_printable_subcommand() const noexcept;
  [[nodiscard]] config parse_arguments();
  [[nodiscard]] int exit(const CLI::ParseError& e) const;
  [[nodiscard]] std::string to_json() const;
  [[nodiscard]] static std::string_view subcommand_to_str(subcommand s) noexcept;

 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  int _exit_code{1};
  config _config{};
  CLI::App _cli{};
  subcommand _subcommand{subcommand::help};

  void make_eval_subcommand();
  void make_find_barrier_clusters_subcommand();
  void make_noisify_subcommand();
  void make_stats_subcommand();
  void make_cli();

  void validate_eval_subcommand() const;
  void validate_find_barrier_clusters_subcommand() const;
  void validate_noisify_subcommand() const;
  void validate_stats_subcommand() const;
  void validate() const;
};

}  // namespace modle::tools
