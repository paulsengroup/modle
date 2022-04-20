// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/debugging/failure_signal_handler.h>  // for InstallFailureSignalHandler, FailureS...
#include <absl/debugging/symbolize.h>               // for InitializeSymbolizer
#include <absl/strings/strip.h>                     // for StripSuffix
#include <absl/time/clock.h>                        // for Now
#include <absl/time/time.h>                         // for FormatDuration, operator-, Time
#include <fmt/format.h>                             // for make_format_args, vformat_to, FMT_STRING
#include <fmt/ostream.h>                            // for formatbuf<>::int_type
#include <spdlog/common.h>                          // for sink_ptr, spdlog_ex, err
#include <spdlog/logger.h>                          // for logger
#include <spdlog/sinks/basic_file_sink.h>           // for basic_file_sink_mt
#include <spdlog/sinks/sink.h>                      // for sink
#include <spdlog/sinks/stdout_color_sinks.h>        // for stderr_color_sink_mt
#include <spdlog/spdlog.h>                          // for error, info

#include <CLI/CLI.hpp>                      // for ParseError
#include <algorithm>                        // for max
#include <boost/filesystem/operations.hpp>  // for create_directories
#include <boost/filesystem/path.hpp>        // for path, operator<<
#include <cassert>                          // for assert
#include <cstdio>                           // for stderr
#include <cstring>                          // for strlen
#include <exception>                        // for exception
#include <iosfwd>                           // for streamsize
#include <memory>                           // for make_shared, __shared_ptr_access, uni...
#include <new>                              // for bad_alloc
#include <stdexcept>                        // for runtime_error
#include <string>                           // for basic_string
#include <string_view>                      // for string_view
#include <vector>                           // for vector

#include "./cli.hpp"                           // for Cli, Cli::subcommand
#include "modle/common/common.hpp"             // for usize
#include "modle/common/simulation_config.hpp"  // for Config
#include "modle/config/version.hpp"            // for str_long
#include "modle/simulation.hpp"                // for Simulation

void setup_logger_console(const bool quiet) {
  spdlog::set_default_logger(std::make_shared<spdlog::logger>("main_logger"));

  auto stderr_sink = spdlog::default_logger()->sinks().emplace_back(
      std::make_shared<spdlog::sinks::stderr_color_sink_mt>());
  //                        [2021-08-12 17:49:34.581] [info]: my log msg
  stderr_sink->set_pattern("[%Y-%m-%d %T.%e] %^[%l]%$: %v");
  if (quiet) {
    stderr_sink->set_level(spdlog::level::err);
  }
  spdlog::info(FMT_STRING("Running MoDLE v{}"), modle::config::version::str());
}

void setup_logger_file(const boost::filesystem::path& path_to_log_file) {
  spdlog::logger("tmp_logger", spdlog::default_logger()->sinks().front())
      .info(FMT_STRING("Complete log will be written to file {}"), path_to_log_file);

  auto file_sink = spdlog::default_logger()->sinks().emplace_back(
      std::make_shared<spdlog::sinks::basic_file_sink_mt>(path_to_log_file.string(), true));
  //                      [2021-08-12 17:49:34.581] [139797797574208] [info]: my log msg
  file_sink->set_pattern("[%Y-%m-%d %T.%e] [%t] %^[%l]%$: %v");

  spdlog::logger("tmp_logger", file_sink)
      .info(FMT_STRING("Running MoDLE v{}"), modle::config::version::str());
}

void setup_failure_signal_handler(const char* argv_0) {
  absl::InitializeSymbolizer(argv_0);
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
}

void write_param_summary_to_log(const modle::Config& c) {
  spdlog::info(FMT_STRING("Command: {}"), fmt::join(c.args, " "));
  spdlog::info(FMT_STRING("Simulation will use up to {} out of {} available CPU cores."),
               c.nthreads, std::thread::hardware_concurrency());
  if (c.stopping_criterion == modle::Config::StoppingCriterion::contact_density) {
    spdlog::info(FMT_STRING("Using --target-contact-density={:.2f} as stopping criterion."),
                 c.target_contact_density);
  } else {
    spdlog::info(FMT_STRING("Using --target-number-of-epochs={} as stropping criterion."),
                 c.target_simulation_epochs);
  }
  spdlog::info(FMT_STRING("Contact sampling strategy: {}."),
               modle::Cli::contact_sampling_strategy_map.at(c.contact_sampling_strategy));
}

int main(int argc, char** argv) noexcept {
  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  setup_failure_signal_handler(argv[0]);

  std::unique_ptr<modle::Cli> cli{nullptr};
  try {
    cli = std::make_unique<modle::Cli>(argc, argv);
    auto config = cli->parse_arguments();
    setup_logger_console(config.quiet);

    if (const auto collisions = cli->detect_path_collisions(config); !collisions.empty()) {
      spdlog::error(FMT_STRING("FAILURE! The following path collision(s) have been detected:\n{}"),
                    collisions);
      return 1;
    }

    if (!config.skip_output) {
      assert(!config.path_to_log_file.empty());
      if (const auto& output_dir = config.path_to_output_prefix.parent_path();
          !output_dir.empty()) {
        boost::filesystem::create_directories(output_dir.string());
      }
      setup_logger_file(config.path_to_log_file);

      if (!cli->config_file_parsed()) {
        spdlog::info(FMT_STRING("Writing simulation parameters to config file {}"),
                     config.path_to_config_file);
        cli->write_config_file();
      }
    }

    write_param_summary_to_log(config);
    const auto t0 = absl::Now();
    modle::Simulation sim(config);
    switch (cli->get_subcommand()) {
      using subcommand = modle::Cli::subcommand;
      case subcommand::simulate:
        sim.run_simulate();
        break;
      case subcommand::perturbate:
        if (config.compute_reference_matrix) {
          sim.run_simulate();
        }
        sim.run_perturbate();
        break;
      case subcommand::replay:
        sim.run_replay();
        break;
      default:
        MODLE_UNREACHABLE_CODE;
    }

    spdlog::info(FMT_STRING("Simulation terminated without errors in {}!\n\nBye."),
                 absl::FormatDuration(absl::Now() - t0));
  } catch (const CLI::ParseError& e) {
    return cli->exit(e);  //  This takes care of formatting and printing error messages (if any)
  } catch (const std::bad_alloc& e) {
    spdlog::error(FMT_STRING("FAILURE! Unable to allocate enough memory: {}"), e.what());
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
