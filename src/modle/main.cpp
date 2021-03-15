#include <absl/strings/str_join.h>
#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration, Time
#include <fmt/format.h>       // for FMT_STRING, format, system_error
#include <fmt/ostream.h>      // for print

#include <CLI/Error.hpp>  // for ParseError
#include <cstdint>        // for uint64_t
#include <cstdio>         // for stderr
#include <exception>
#include <filesystem>  // for exists, weakly_canonical, path
#include <iosfwd>      // for ofstream
#include <stdexcept>   // for runtime_error
#include <vector>

#include "./cli.hpp"            // for Cli
#include "modle/chr_sizes.hpp"  // for ChrSize, Parser
#include "modle/config.hpp"     // for config
#include "modle/genome.hpp"     // for Genome

namespace modle {

std::vector<std::string> check_if_output_file_exists(const modle::config& c) {
  std::vector<std::string> fn_collisions;
  if (c.force) {
    return fn_collisions;
  }
  if (std::filesystem::exists(c.path_to_output_file)) {
    fn_collisions.push_back(std::filesystem::weakly_canonical(c.path_to_output_file));
  }
  return fn_collisions;
}

void run_simulation(const modle::config& c) {
  const auto t0 = absl::Now();
  modle::Genome genome(c);

  if (!c.skip_output) {  // Write simulation params to file
    if (c.force) {
      std::filesystem::remove_all(c.path_to_output_file);
    }
    std::filesystem::create_directories(std::filesystem::path(c.path_to_output_file).parent_path());
    std::ofstream log_file(c.path_to_log_file);
    if (log_file) {
      fmt::print(log_file, FMT_STRING("{}\n{}\n"), c.to_string(),
                 absl::StrJoin(c.argv, c.argv + c.argc, " "));
    } else {
      fmt::print(
          stderr,
          FMT_STRING("WARNING: Unable to open log file {} for writing. Continuing anyway..."),
          c.path_to_log_file);
    }
  }
  genome.simulate_extrusion(c.skip_output ? "" : c.path_to_output_file, c.ncells,
                            c.simulation_iterations, c.nthreads);
  fmt::print(stderr, FMT_STRING("Simulation took {}.\n"), absl::FormatDuration(absl::Now() - t0));
  fmt::print(stderr, "Simulation terminated without errors!\n\nBye.\n");
}

}  // namespace modle

int main(int argc, char** argv) noexcept {
  auto cli = modle::Cli(argc, argv);

  try {
    auto config = cli.parse_arguments();
    if (const auto collisions = modle::Cli::process_paths_and_check_for_collisions(config);
        !collisions.empty()) {
      fmt::print(stderr, FMT_STRING("The following path collision(s) have been detected:\n{}"),
                 collisions);
      return 1;
    }
    config.print();

    modle::run_simulation(config);
  } catch (const CLI::ParseError& e) {
    return cli.exit(e);  //  This takes care of formatting and printing error messages (if any)
  } catch (const fmt::system_error& err) {
    fmt::print(stderr, "FAILURE! An error occurred during simulation: {}.\n", err.what());
    return 1;
  } catch (const std::runtime_error& err) {
    fmt::print(stderr, "FAILURE! An error occurred during simulation: {}.\n", err.what());
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
