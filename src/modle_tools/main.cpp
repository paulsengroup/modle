// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/debugging/failure_signal_handler.h>  // for InstallFailureSignalHandler, FailureS...
#include <absl/debugging/symbolize.h>               // for InitializeSymbolizer
#include <absl/strings/strip.h>                     // for StripSuffix
#include <absl/types/variant.h>                     // for get
#include <fmt/format.h>                             // for make_format_args, vformat_to, FMT_STRING
#include <spdlog/common.h>                          // for sink_ptr, spdlog_ex
#include <spdlog/logger.h>                          // for logger
#include <spdlog/sinks/sink.h>                      // for sink
#include <spdlog/sinks/stdout_color_sinks.h>        // for stderr_color_sink_mt
#include <spdlog/spdlog.h>                          // for error

#include <CLI/CLI.hpp>  // for ParseError
#include <algorithm>    // for max
#include <cassert>      // for assert
#include <cstdio>       // for stderr
#include <cstring>      // for strlen
#include <exception>    // for exception
#include <memory>       // for unique_ptr, make_shared, make_unique
#include <new>          // for bad_alloc
#include <stdexcept>    // for runtime_error
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "./cli.hpp"              // for Cli, Cli::subcommand, Cli::eval, Cli:...
#include "modle_tools/tools.hpp"  // for eval_subcmd, find_barrier_clusters_su...

int main(int argc, char** argv) {
  using namespace modle::tools;
  absl::InitializeSymbolizer(argv[0]);
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
  std::unique_ptr<Cli> cli{nullptr};
  spdlog::set_default_logger(std::make_shared<spdlog::logger>("main_logger"));
  {
    auto stderr_sink = spdlog::default_logger()->sinks().emplace_back(
        std::make_shared<spdlog::sinks::stderr_color_sink_mt>());
    //                        [2021-08-12 17:49:34.581] [info]: my log msg
    stderr_sink->set_pattern("[%Y-%m-%d %T.%e] %^[%l]%$: %v");
  }
  try {
    cli = std::make_unique<Cli>(argc, argv);
    const auto config = cli->parse_arguments();

    // TODO: move this inside the subcommands
    /*
    if (!config.output_prefix.empty()) {
      boost::filesystem::create_directories(config.output_prefix.parent_path());
    }
     */

    {
      using subcmd = Cli::subcommand;
      switch (cli->get_subcommand()) {
        case subcmd::eval:
          eval_subcmd(absl::get<eval_config>(config));
          return 0;
        case subcmd::fbcl:
          find_barrier_clusters_subcmd(absl::get<find_barrier_clusters_config>(config));
          return 0;
        case subcmd::noisify:
          noisify_subcmd(absl::get<noisify_config>(config));
          return 0;
        case subcmd::transform:
          transform_subcmd(absl::get<transform_config>(config));
          return 0;
        default:
          throw std::runtime_error(
              "Default branch in switch statement in modle_tools::main() should be unreachable! If "
              "you see this message, please file an issue on GitHub");
      }
    }
  } catch (const CLI::ParseError& e) {
    assert(cli);
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
    assert(cli);
    spdlog::error(FMT_STRING("FAILURE! modle_tools {} encountered the following error: {}."),
                  cli->get_printable_subcommand(), e.what());
    return 1;
  } catch (...) {
    spdlog::error(FMT_STRING("FAILURE! modle_tools {} encountered the following error: Caught an "
                             "unhandled exception! "
                             "If you see this message, please file an issue on GitHub."),
                  cli->get_printable_subcommand());
    throw;
    return 1;
  }
  return 0;
}
