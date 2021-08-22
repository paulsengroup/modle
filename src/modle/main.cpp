#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration, operator-, Time
#include <fmt/format.h>       // for print
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

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

#include "./cli.hpp"                // for Cli
#include "modle/common/config.hpp"  // for Config
#include "modle/common/utils.hpp"   // for traced
#include "modle/simulation.hpp"     // for Simulation

int main(int argc, char** argv) noexcept {
  std::unique_ptr<modle::Cli> cli{nullptr};
  spdlog::set_default_logger(std::make_shared<spdlog::logger>("main_logger"));
  {
    auto stderr_sink = spdlog::default_logger()->sinks().emplace_back(
        std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
    //                        [2021-08-12 17:49:34.581] [info]: my log msg
    stderr_sink->set_pattern("[%Y-%m-%d %T.%e] %^[%l]%$: %v");
  }

  try {
    cli = std::make_unique<modle::Cli>(argc, argv);
    auto config = cli->parse_arguments();
    if (const auto collisions = cli->detect_path_collisions(config); !collisions.empty()) {
      spdlog::error(FMT_STRING("FAILURE! The following path collision(s) have been detected:\n{}"),
                    collisions);
      return 1;
    }
    if (config.quiet) {
      assert(spdlog::default_logger()->sinks().size() == 1);  // NOLINT
      spdlog::default_logger()->sinks().front()->set_level(spdlog::level::err);
    }

    const auto t0 = absl::Now();
    if (!config.skip_output) {
      assert(!config.path_to_log_file.empty());  // NOLINT
      boost::filesystem::create_directories(config.path_to_output_prefix.parent_path());
      spdlog::logger("tmp_logger", spdlog::default_logger()->sinks().front())
          .info(FMT_STRING("Complete log will be written to file {}"), config.path_to_log_file);
      auto file_sink = spdlog::default_logger()->sinks().emplace_back(
          std::make_shared<spdlog::sinks::basic_file_sink_mt>(config.path_to_log_file.string(),
                                                              true));
      //                      [2021-08-12 17:49:34.581] [139797797574208] [info]: my log msg
      file_sink->set_pattern("[%Y-%m-%d %T.%e] [%t] %^[%l]%$: %v");
      if (!cli->config_file_parsed()) {
        spdlog::info(FMT_STRING("Writing simulation parameters to config file {}"),
                     config.path_to_config_file);
        cli->write_config_file();
      }
    }
    spdlog::info(FMT_STRING("Command: {}"),
                 absl::StrJoin(config.argv, config.argv + config.argc, " "));
    modle::Simulation sim(config);
    switch (cli->get_subcommand()) {
      using subcommand = modle::Cli::subcommand;
      case subcommand::simulate:
        sim.run_simulate();
        break;
      case subcommand::pertubate:
        if (config.compute_reference_matrix) {
          sim.run_simulate();
        }
        sim.run_perturbate();
        break;
      case subcommand::replay:
        sim.run_replay();
        break;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in modle::main() should be unreachable! If you see "
            "this message, please file an issue on GitHub");
    }
    spdlog::info(FMT_STRING("Simulation terminated without errors in {}!\n\nBye."),
                 absl::FormatDuration(absl::Now() - t0));
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
#ifndef BOOST_STACKTRACE_USE_NOOP
    const auto* st = boost::get_error_info<modle::utils::traced>(e);
    if (st) {
      spdlog::error(*st);
    } else {
      spdlog::error("Stack trace not available!");
    }
#endif
    return 1;
  } catch (...) {
    spdlog::error(
        FMT_STRING("FAILURE! An error occurred during simulation: Caught an unhandled exception! "
                   "If you see this message, please file an issue on GitHub."));
    return 1;
  }
  return 0;
}
