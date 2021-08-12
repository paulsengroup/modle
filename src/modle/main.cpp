#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration, operator-, Time
#include <fmt/format.h>       // for print, system_error, FMT_STRING

#include <CLI/Error.hpp>  // for ParseError
#include <cstdio>         // for stderr
#include <exception>      // for current_exception, exception_ptr, rethrow_...
#include <iostream>       // for operator<<, basic_ostream, cerr, ostream
#include <memory>         // for unique_ptr, make_unique
#include <new>            // for bad_alloc
#include <stdexcept>      // for runtime_error
#include <string>         // for basic_string

#ifndef BOOST_STACKTRACE_USE_NOOP
#include <boost/exception/get_error_info.hpp>  // for get_error_info
#include <boost/stacktrace/stacktrace.hpp>     // for operator<<
#include <iostream>                            // for operator<<, basic_ostream, cerr, ostream
#endif

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "./cli.hpp"                // for Cli
#include "modle/common/config.hpp"  // for Config
#include "modle/common/utils.hpp"   // for traced
#include "modle/simulation.hpp"     // for Simulation

/*
[[nodiscard]] std::pair<spdlog::logger, spdlog::logger> setup_loggers(
    const boost::filesystem::path& path_to_log_file, const bool quiet) {
  auto file_logger = spdlog::basic_logger_mt("file_logger", path_to_log_file.string(), true);
  auto stderr_logger = spdlog::stderr_color_mt("stderr_logger");

  return std::make_pair(stderr_logger, file_logger);
}
 */

int main(int argc, char** argv) noexcept {
  std::unique_ptr<modle::Cli> cli{nullptr};

  try {
    cli = std::make_unique<modle::Cli>(argc, argv);
    auto config = cli->parse_arguments();
    if (const auto collisions = modle::Cli::process_paths_and_check_for_collisions(config);
        !collisions.empty()) {
      fmt::print(stderr, FMT_STRING("The following path collision(s) have been detected:\n{}"),
                 collisions);
      return 1;
    }
    config.print();

    const auto t0 = absl::Now();
    if (!config.skip_output) {
      boost::filesystem::create_directories(config.path_to_output_prefix.parent_path());
      std::ofstream log_file(config.path_to_log_file.string());
      if (log_file) {
        fmt::print(log_file, FMT_STRING("{}\n{}\n"),
                   absl::StrJoin(config.argv, config.argv + config.argc, " "), config.to_string());
      } else {
        fmt::print(
            stderr,
            FMT_STRING("WARNING: Unable to open log file {} for writing. Continuing anyway..."),
            config.path_to_log_file);
      }
    }
    modle::Simulation sim(config);
    switch (cli->get_subcommand()) {
      case modle::Cli::subcommand::simulate:
        sim.run_simulation();
        break;
      case modle::Cli::subcommand::pertubate:
        if (config.compute_reference_matrix) {
          sim.run_simulation();
        }
        sim.run_perturbate();
        break;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in modle::main() should be unreachable! If you see "
            "this message, please file an issue on GitHub");
    }
    fmt::print(stderr, FMT_STRING("Simulation terminated without errors in {}!\n"),
               absl::FormatDuration(absl::Now() - t0));
    fmt::print(stderr, FMT_STRING("\nBye.\n"));
  } catch (const CLI::ParseError& e) {
    return cli->exit(e);  //  This takes care of formatting and printing error messages (if any)
  } catch (const std::bad_alloc& err) {
    fmt::print(stderr, "FAILURE! Unable to allocate enough memory.\n");
    return 1;
  } catch (const std::exception& e) {
    fmt::print(stderr, "FAILURE! An error occurred during simulation: {}.\n", e.what());
#ifndef BOOST_STACKTRACE_USE_NOOP
    const auto* st = boost::get_error_info<modle::utils::traced>(e);
    if (st) {
      std::cerr << *st << '\n';
    } else {
      fmt::print(stderr, "Stack trace not available!\n");
    }
#endif
    return 1;
  } catch (...) {
    const auto err = std::current_exception();
    auto handle_except = [&]() {
      try {
        if (err) {
          std::rethrow_exception(err);
        }
      } catch (const std::exception& e) {
        fmt::print(stderr,
                   "FAILURE! An error occurred during simulation: Caught an unhandled exception! "
                   "If you see this message, please file an issue on GitHub.\nerr.what(): {}\n",
                   e.what());
        return 1;
      }
      return 0;
    };
    return handle_except();
  }
  return 0;
}
