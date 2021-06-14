#include <fmt/format.h>  // for print

#include <cstdio>        // for stderr
#include <exception>     // for exception
#include <filesystem>    // for create_directories, exists, is_empty, remove_all, path
#include <stdexcept>     // for runtime_error
#include <system_error>  // for error_code

#include "./cli.hpp"               // for Cli, Cli::subcommand, Cli::eval, Cli::stats
#include "modle_tools/config.hpp"  // for config
#include "modle_tools/tools.hpp"   // for eval_subcmd, stats_subcmd

int main(int argc, char** argv) {
  modle::tools::config c;
  try {
    modle::tools::Cli cli(argc, argv);
    c = cli.parse_arguments();
    if (!cli.is_ok()) {
      return cli.get_exit_code();
    }

    if (!c.output_base_name.empty()) {
      std::filesystem::create_directories(c.output_base_name.parent_path());
    }

    switch (cli.get_subcommand()) {
      case modle::tools::Cli::subcommand::eval:
        modle::tools::eval_subcmd(c);
        break;
      case modle::tools::Cli::subcommand::filter_barriers:
        modle::tools::filter_barriers_subcmd(c);
        break;
      case modle::tools::Cli::subcommand::stats:
        modle::tools::stats_subcmd(c);
        break;
      case modle::tools::Cli::subcommand::noisify:
        modle::tools::noisify_subcmd(c);
        break;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in modle_tools::main() should be unreachable! If "
            "you see this message, please file an issue on GitHub");
    }
  } catch (const std::exception& err) {
    fmt::print(stderr, "FAILURE: {}.\n", err.what());
    std::error_code ec;
    if (!c.keep_tmp_files && std::filesystem::exists(c.tmp_dir) &&
        std::filesystem::is_empty(c.tmp_dir)) {
      std::filesystem::remove_all(c.tmp_dir, ec);
    }
    return 1;
  } catch (...) {
    fmt::print(stderr,
               "FAILURE: modle_tools encountered an unknown error! If you see this message, please "
               "file an issue on GitHub.");
    return 1;
  }

  return 0;
}
