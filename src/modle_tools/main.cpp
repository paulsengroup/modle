#include <absl/debugging/failure_signal_handler.h>
#include <absl/debugging/symbolize.h>
#include <fmt/format.h>  // for print
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <boost/filesystem/operations.hpp>  // for create_directories, exists, is_empty, remove_all
#include <boost/filesystem/path.hpp>        // for path
#include <cstdio>                           // for stderr
#include <exception>                        // for exception
#include <stdexcept>                        // for runtime_error
#include <system_error>                     // for error_code

#include "./cli.hpp"               // for Cli, Cli::subcommand, Cli::eval, Cli::stats
#include "modle_tools/config.hpp"  // for config
#include "modle_tools/tools.hpp"   // for eval_subcmd, stats_subcmd

int main(int argc, char** argv) {
  absl::InitializeSymbolizer(argv[0]);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  absl::FailureSignalHandlerOptions options;
  // TODO: figure out a way to make this callback async-signal-safe
  options.writerfn = [](const char* buff) {
    if (buff) {
      const std::string_view buff_{buff, strlen(buff)};
      spdlog::error(FMT_STRING("{}"), absl::StripSuffix(buff_, "\n"));
    } else {
      spdlog::shutdown();
    }
  };
  absl::InstallFailureSignalHandler(options);
  std::unique_ptr<modle::tools::Cli> cli{nullptr};
  spdlog::set_default_logger(std::make_shared<spdlog::logger>("main_logger"));
  {
    auto stderr_sink = spdlog::default_logger()->sinks().emplace_back(
        std::make_shared<spdlog::sinks::stderr_color_sink_mt>());
    //                        [2021-08-12 17:49:34.581] [info]: my log msg
    stderr_sink->set_pattern("[%Y-%m-%d %T.%e] %^[%l]%$: %v");
  }
  try {
    cli = std::make_unique<modle::tools::Cli>(argc, argv);
    const auto config = cli->parse_arguments();

    // TODO: move this inside the subcommands
    /*
    if (!config.output_base_name.empty()) {
      boost::filesystem::create_directories(config.output_base_name.parent_path());
    }
     */

    switch (cli->get_subcommand()) {
      case modle::tools::Cli::subcommand::eval:
        modle::tools::eval_subcmd(absl::get<modle::tools::eval_config>(config));
        return 0;
      case modle::tools::Cli::subcommand::filter_barriers:
        modle::tools::filter_barriers_subcmd(
            absl::get<modle::tools::filter_barrier_config>(config));
        return 0;
      case modle::tools::Cli::subcommand::find_barrier_clusters:
        modle::tools::find_barrier_clusters_subcmd(
            absl::get<modle::tools::find_barrier_clusters_config>(config));
        return 0;
      case modle::tools::Cli::subcommand::stats:
        modle::tools::stats_subcmd(absl::get<modle::tools::stats_config>(config));
        return 0;
      case modle::tools::Cli::subcommand::noisify:
        modle::tools::noisify_subcmd(absl::get<modle::tools::noisify_config>(config));
        return 0;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in modle_tools::main() should be unreachable! If "
            "you see this message, please file an issue on GitHub");
    }
  } catch (const CLI::ParseError& e) {
    return cli->exit(e);  //  This takes care of formatting and printing error messages (if any)
  } catch (const std::bad_alloc& err) {
    spdlog::error(FMT_STRING("FAILURE! Unable to allocate enough memory: {}"), err.what());
    return 1;
  } catch (const spdlog::spdlog_ex& e) {
    fmt::print(
        stderr,
        FMT_STRING(
            "FAILURE! An error occurred while setting up the main application logger: {}.\n"),
        e.what());
    return 1;
  } catch (const std::exception& e) {
    spdlog::error(FMT_STRING("FAILURE! An error occurred during simulation: {}."), e.what());
    return 1;
  } catch (...) {
    spdlog::error(
        FMT_STRING("FAILURE! An error occurred during simulation: Caught an unhandled exception! "
                   "If you see this message, please file an issue on GitHub."));
    return 1;
  }
  return 0;
}
