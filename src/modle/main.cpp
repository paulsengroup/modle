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

#include "./cli.hpp"               // for Cli
#include "modle/chr_sizes.hpp"     // for ChrSize, Parser
#include "modle/config.hpp"        // for config
#include "modle/extr_barrier.hpp"  // for ExtrusionBarrier
#include "modle/genome.hpp"        // for Genome

namespace modle {

std::vector<std::string> check_if_output_file_exists(const modle::config& c) {
  std::vector<std::string> fn_collisions;
  if (c.force) {
    return fn_collisions;
  }
  modle::chr_sizes::Parser parser(c.path_to_chr_sizes);
  for (const auto& record : parser.parse_all()) {
    const auto f1 = fmt::format("{}/{}.tsv.bz2", c.path_to_output_file, record.name);
    const auto f2 = fmt::format("{}/{}_raw.tsv.bz2", c.path_to_output_file, record.name);
    if (std::filesystem::exists(f1)) {
      fn_collisions.push_back(std::filesystem::weakly_canonical(f1));
    }
    if (std::filesystem::exists(f2)) {
      fn_collisions.push_back(std::filesystem::weakly_canonical(f2));
    }
  }

  return fn_collisions;
}

void run_simulation(const modle::config& c) {
  auto t0 = absl::Now();
  modle::Genome genome(c);

  uint64_t tot_barriers = c.number_of_randomly_gen_extr_barriers;
  uint64_t barriers_ignored = 0;
  if (c.number_of_randomly_gen_extr_barriers != 0) {
    genome.randomly_generate_extrusion_barriers(c.number_of_randomly_gen_extr_barriers, c.seed);
  }
  if (!c.path_to_extr_barriers_bed.empty()) {
    const auto tmp = genome.import_extrusion_barriers_from_bed(
        c.path_to_extr_barriers_bed, c.probability_of_extrusion_barrier_block);
    tot_barriers += tmp.first;
    barriers_ignored = tmp.second;
  }
  genome.sort_extr_barriers_by_pos();
  if (c.exclude_chr_wo_extr_barriers) {
    genome.exclude_chr_wo_extr_barriers();
  }
  const auto n_of_chr_removed =
      std::accumulate(genome.get_chromosomes().begin(), genome.get_chromosomes().end(), 0UL,
                      [](std::size_t accumulator, const auto& chr) {
                        return accumulator + chr.get_nbarriers() == 0;
                      });
  if (genome.get_nchromosomes() == 0) {
    throw std::runtime_error(  // TODO: Improve this error message
        "All the input sequences were discarded because there were no extrusion barriers mapping "
        "on them. Check your input files");
  }
  fmt::print(stderr,
             FMT_STRING("Initialization took {}\n"
                        " - # of sequences:       {} ({} ignored)\n"
                        " - Avg. sequence length: {:.3f} Mbp\n"
                        " - Genome N50:           {:.3f} Mbp\n"
                        " - # of LEFs:            {}\n"
                        " - # of extr. barriers   {} ({} ignored)\n\n"),
             absl::FormatDuration(absl::Now() - t0), genome.get_nchromosomes(), n_of_chr_removed,
             (static_cast<double>(genome.size()) / genome.get_nchromosomes()) / 1.0e6,
             static_cast<double>(genome.n50()) / 1.0e6, genome.get_nlefs(),
             tot_barriers - barriers_ignored, barriers_ignored);

  t0 = absl::Now();
  // Assign LEFs and bind them to a random pos if skip_burn in is true
  genome.assign_lefs(c.skip_burnin);

  if (!c.skip_output) {  // Write simulation params to file
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

  if (!c.skip_burnin) {
    fmt::print(stderr, "Running burnin phase...\n");
    const auto& [avg_burnin_rounds, burnin_rounds_stdev] = genome.run_burnin(
        c.probability_of_lef_rebind, c.min_n_of_loops_per_lef, c.min_n_of_burnin_rounds);
    fmt::print(stderr, FMT_STRING("Burnin completed in {:.0f} rounds{}.\n"), avg_burnin_rounds,
               burnin_rounds_stdev == 0
                   ? ""
                   : fmt::format(FMT_STRING(" on average (stdev = {:.2f})"), burnin_rounds_stdev));
  }
  fmt::print(stderr, FMT_STRING("Bound {} LEFs in {}!\n"), genome.get_n_busy_lefs(),
             absl::FormatDuration(absl::Now() - t0));

  t0 = absl::Now();
  fmt::print(stderr, "About to start simulating loop extrusion...\n");
  if (c.target_contact_density != 0) {
    genome.simulate_extrusion(c.target_contact_density);
  } else {
    genome.simulate_extrusion(c.simulation_iterations);
  }
  fmt::print(stderr, FMT_STRING("Simulation took {}.\n"), absl::FormatDuration(absl::Now() - t0));

  if (!c.skip_output) {  // Mostly useful for profiling
    if (c.force) {
      std::filesystem::remove_all(c.path_to_output_file);
    }
    if (c.write_raw_contacts) {
      genome.write_contacts_to_file(c.path_to_output_file, c.write_contacts_for_ko_chroms);
    }
    if (c.write_contacts_with_noise) {
      genome.write_contacts_w_noise_to_file(c.path_to_output_file_w_noise, c.random_noise_mean,
                                            c.random_noise_std, c.write_contacts_for_ko_chroms);
    }
  }
  fmt::print(stderr, "Simulation terminated without errors!\n\nBye.\n");
}

}  // namespace modle

int main(int argc, char** argv) noexcept {
  auto cli = modle::Cli(argc, argv);

  try {
    auto config = cli.parse_arguments();
    if (const auto collisions = cli.process_paths_and_check_for_collisions(config);
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
