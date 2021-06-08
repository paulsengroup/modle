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

#include "./cli.hpp"                // for Cli
#include "modle/common/config.hpp"  // for Config
#include "modle/common/utils.hpp"   // for traced
#include "modle/simulation.hpp"     // for Simulation

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

    modle::Simulation{config}.run();

    fmt::print(stderr, FMT_STRING("Simulation terminated without errors in {}!\n\nBye.\n"),
               absl::FormatDuration(absl::Now() - t0));
  } catch (const CLI::ParseError& e) {
    return cli->exit(e);  //  This takes care of formatting and printing error messages (if any)
  } catch (const fmt::system_error& err) {
    fmt::print(stderr, "FAILURE! An error occurred during simulation: {}.\n", err.what());
    return 1;
  } catch (const std::runtime_error& err) {
    fmt::print(stderr, "FAILURE! An error occurred during simulation: {}.\n", err.what());
    return 1;
  } catch (const std::bad_alloc& err) {
    fmt::print(stderr, "FAILURE! Unable to allocate enough memory.\n");
    return 1;
  } catch (const std::exception& err) {
    fmt::print(stderr, FMT_STRING("{}\n"), err.what());
#ifndef BOOST_STACKTRACE_USE_NOOP
    const auto* st = boost::get_error_info<modle::utils::traced>(err);
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
                   "FAILURE! An error occurred during simulation: Caught an exception that was not "
                   "handled properly! If you see this message, please open an issue on GitHub. "
                   "err.what(): {}.\n",
                   e.what());
        return 1;
      }
      return 0;
    };
    return handle_except();
  }
  return 0;
}
