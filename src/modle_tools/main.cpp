#include <fmt/format.h>  // for print

#include <boost/filesystem/operations.hpp>  // for create_directories, exists, is_empty, remove_all
#include <boost/filesystem/path.hpp>        // for path
#include <boost/system/error_code.hpp>      // for system::error_code
#include <cstdio>                           // for stderr
#include <exception>                        // for exception
#include <stdexcept>                        // for runtime_error
#include <system_error>                     // for error_code

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

    // TODO: move this inside the subcommands
    /*
    if (!c.output_base_name.empty()) {
      boost::filesystem::create_directories(c.output_base_name.parent_path());
    }
     */

    switch (cli.get_subcommand()) {
      case modle::tools::Cli::subcommand::eval:
        modle::tools::eval_subcmd(absl::get<modle::tools::eval_config>(c));
        break;
      case modle::tools::Cli::subcommand::filter_barriers:
        modle::tools::filter_barriers_subcmd(absl::get<modle::tools::filter_barrier_config>(c));
        break;
      case modle::tools::Cli::subcommand::stats:
        modle::tools::stats_subcmd(absl::get<modle::tools::stats_config>(c));
        break;
      case modle::tools::Cli::subcommand::noisify:
        modle::tools::noisify_subcmd(absl::get<modle::tools::noisify_config>(c));
        break;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in modle_tools::main() should be unreachable! If "
            "you see this message, please file an issue on GitHub");
    }
  } catch (const std::exception& err) {
    fmt::print(stderr, "FAILURE: {}.\n", err.what());
    // TODO Fixme
    /*
    boost::system::error_code ec;
    if (!c.keep_tmp_files && boost::filesystem::exists(c.tmp_dir) &&
        boost::filesystem::is_empty(c.tmp_dir)) {
      boost::filesystem::remove_all(c.tmp_dir, ec);
    }*/
    return 1;
  } catch (...) {
    fmt::print(stderr,
               "FAILURE: modle_tools encountered an unknown error! If you see this message, please "
               "file an issue on GitHub.");
    return 1;
  }

  return 0;
}
