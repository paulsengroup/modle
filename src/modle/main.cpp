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
#include <boost/filesystem/exception.hpp>   // for filesystem::filesystem_exception
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

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::atomic<bool> logger_ready{false};

void setup_logger_console(const bool quiet) {
  auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  //                        [2021-08-12 17:49:34.581] [info]: my log msg
  stderr_sink->set_pattern("[%Y-%m-%d %T.%e] %^[%l]%$: %v");
  if (quiet) {
    stderr_sink->set_level(spdlog::level::err);
  }

  auto main_logger = std::make_shared<spdlog::logger>("main_logger", stderr_sink);

  spdlog::set_default_logger(main_logger);
  spdlog::info(FMT_STRING("Running MoDLE v{}"), modle::config::version::str());

  logger_ready = true;
}

void setup_logger_file(const boost::filesystem::path& path_to_log_file) {
  spdlog::info(FMT_STRING("Complete log will be written to file {}"), path_to_log_file);

  auto file_sink =
      std::make_shared<spdlog::sinks::basic_file_sink_mt>(path_to_log_file.string(), true);
  //                      [2021-08-12 17:49:34.581] [139797797574208] [info]: my log msg
  file_sink->set_pattern("[%Y-%m-%d %T.%e] [%t] %^[%l]%$: %v");

  spdlog::logger("tmp_logger", file_sink)
      .info(FMT_STRING("Running MoDLE v{}"), modle::config::version::str());

  spdlog::default_logger()->sinks().emplace_back(std::move(file_sink));

  logger_ready = true;
}

void setup_failure_signal_handler(const char* argv_0) {
  absl::InitializeSymbolizer(argv_0);
  absl::FailureSignalHandlerOptions options;
  // TODO: figure out a way to make this callback async-signal-safe
  options.writerfn = [](const char* buff) {
    if (buff) {
      const std::string_view buff_sv{buff, strlen(buff)};
      if (logger_ready) {
        spdlog::error(FMT_STRING("{}"), absl::StripSuffix(buff_sv, "\n"));
      } else {
        fmt::print(stderr, FMT_STRING("{}"), buff_sv);
      }
    }
    spdlog::shutdown();
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
  spdlog::info(FMT_STRING("Contact matrix resolution: {}bp"), c.bin_size);
}

std::tuple<int, modle::Cli::subcommand, modle::Config> parse_cli_and_setup_logger(int argc,
                                                                                  char** argv) {
  std::unique_ptr<modle::Cli> cli{nullptr};
  try {
    cli = std::make_unique<modle::Cli>(argc, argv);
    auto config = cli->parse_arguments();
    const auto subcmd = cli->get_subcommand();
    setup_logger_console(config.quiet);

    if (const auto collisions = cli->detect_path_collisions(config); !collisions.empty()) {
      throw boost::filesystem::filesystem_error(
          fmt::format(FMT_STRING("Detected the following path collision(s):\n{}"), collisions.size(),
                      collisions),
          std::make_error_code(std::errc::file_exists));
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

    cli->log_warnings();

    return std::make_tuple(0, subcmd, config);
    // NOTE: GCC7 crashes if using modle::Config{} instead of modle::Config() in the catch blocks
    // below
  } catch (const CLI::ParseError& e) {
    assert(cli);
    //  This takes care of formatting and printing error messages (if any)
    return std::make_tuple(cli->exit(e), modle::Cli::subcommand::help, modle::Config());

  } catch (const boost::filesystem::filesystem_error& e) {
    spdlog::error(FMT_STRING("FAILURE! One or more filesystem error(s) occurred: {}."), e.what());
    return std::make_tuple(1, modle::Cli::subcommand::help, modle::Config());
  } catch (const spdlog::spdlog_ex& e) {
    fmt::print(
        stderr,
        FMT_STRING(
            "FAILURE! An error occurred while setting up the main application logger: {}.\n"),
        e.what());
    return std::make_tuple(1, modle::Cli::subcommand::help, modle::Config());
  }
}

template <typename... Args>
void try_log_fatal_error(fmt::format_string<Args...> fmt, Args&&... args) {
  if (logger_ready) {
    assert(spdlog::default_logger());
    spdlog::error(fmt, std::forward<Args>(args)...);
  } else {
    fmt::print(stderr, fmt, std::forward<Args>(args)...);
  }
}

int main(int argc, char** argv) noexcept {
  // No need to set up the signal handler w/ symbol support when proj. is built in Debug mode
  if constexpr (modle::utils::ndebug_defined()) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    setup_failure_signal_handler(argv[0]);
  }

  try {
    const auto [ec, subcmd, config] = parse_cli_and_setup_logger(argc, argv);
    if (ec != 0 || subcmd == modle::Cli::subcommand::help) {
      return ec;
    }

    assert(spdlog::default_logger());
    write_param_summary_to_log(config);
    const auto t0 = absl::Now();
    modle::Simulation sim(config);
    switch (subcmd) {
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
  } catch (const std::bad_alloc& e) {
    try_log_fatal_error(FMT_STRING("FAILURE! Unable to allocate enough memory: {}."), e.what());
    return 1;
  } catch (const std::exception& e) {
    try_log_fatal_error(FMT_STRING("FAILURE! An error occurred during simulation: {}."), e.what());
    return 1;
  } catch (...) {
    try_log_fatal_error(
        FMT_STRING("FAILURE! An error occurred during simulation: Caught an unhandled exception! "
                   "If you see this message, please file an issue on GitHub."));
    return 1;
  }
  return 0;
}
