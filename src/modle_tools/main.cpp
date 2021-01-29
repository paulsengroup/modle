#include <fmt/format.h>  // for print

#include <cstdio>        // for stderr
#include <filesystem>    // for create_directories, is_empty, remove_all, exists, filesystem_error
#include <stdexcept>     // for runtime_error
#include <system_error>  // for error_code

#include "./cli.hpp"              // for config, Cli, Cli::subcommand, Cli::convert, Cli::eval
#include "modle_tools/tools.hpp"  // for convert, eval

int main(int argc, char** argv) {
  modle::tools::Cli cli(argc, argv);
  auto c = cli.parse_arguments();
  if (!cli.is_ok()) {
    return c.exit_code;
  }

  if (!c.output_base_name.empty()) {
    std::filesystem::create_directories(c.output_base_name);
  }
  try {
    switch (cli.get_subcommand()) {
      case modle::tools::Cli::subcommand::eval:
        modle::tools::eval_subcmd(c);
        break;
      case modle::tools::Cli::subcommand::stats:
        modle::tools::stats_subcmd(c);
        break;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in main.cpp should be unreachable! If you see "
            "this message, please open an issue on GitHub");
    }
  } catch (const std::runtime_error& err) {
    fmt::print(stderr, "FAILURE: {}.\n", err.what());
    std::error_code ec;
    if (!c.keep_tmp_files && std::filesystem::is_empty(c.tmp_dir)) {
      std::filesystem::remove_all(c.tmp_dir, ec);
    }
    return 1;
  }

  return 0;
}
