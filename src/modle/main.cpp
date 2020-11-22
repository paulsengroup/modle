#include <filesystem>
#include <string>

#include "./cli.hpp"
#include "absl/time/clock.h"
#include "fmt/printf.h"
#include "modle/chr_sizes.hpp"
#include "modle/genome.hpp"

namespace modle {

std::vector<std::string> check_if_output_file_exists(const modle::config& c) {
  std::vector<std::string> fn_collisions;
  if (c.force) {
    return fn_collisions;
  }
  modle::chr_sizes::Parser parser(c.path_to_chr_sizes);
  for (const auto& record : parser.parse()) {
    const auto f1 = fmt::format("{}/{}.tsv.bz2", c.output_dir, record.name);
    const auto f2 = fmt::format("{}/{}_raw.tsv.bz2", c.output_dir, record.name);
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
    genome.randomly_generate_extrusion_barriers(c.number_of_randomly_gen_extr_barriers);
  }
  if (!c.path_to_extr_barriers_bed.empty()) {
    const auto tmp = genome.import_extrusion_barriers_from_bed(
        c.path_to_extr_barriers_bed, c.probability_of_extrusion_barrier_block);
    tot_barriers += tmp.first;
    barriers_ignored = tmp.second;
  }
  const auto n_of_chr_removed = genome.remove_chromosomes_wo_extr_barriers();
  if (genome.get_n_chromosomes() == 0) {
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
             absl::FormatDuration(absl::Now() - t0), genome.get_n_chromosomes(), n_of_chr_removed,
             (static_cast<double>(genome.size()) / genome.get_n_chromosomes()) / 1.0e6,
             static_cast<double>(genome.n50()) / 1.0e6, genome.get_n_lefs(),
             tot_barriers - barriers_ignored, barriers_ignored);

  t0 = absl::Now();
  if (c.skip_burnin) {
    genome.randomly_bind_lefs();
    fmt::print(stderr, FMT_STRING("Bound {} LEFs in {}.\n"), genome.get_n_of_busy_lefs(),
               absl::FormatDuration(absl::Now() - t0));
  } else {
    const auto burnin_rounds = genome.run_burnin(
        c.probability_of_lef_rebind, c.min_n_of_loops_per_lef, c.min_n_of_burnin_rounds);
    fmt::print(stderr, FMT_STRING("Burnin completed in {}! ({} rounds).\n"),
               absl::FormatDuration(absl::Now() - t0), burnin_rounds);
  }

  t0 = absl::Now();
  fmt::print(stderr, "About to start simulating loop extrusion...\n");
  genome.simulate_extrusion(c.simulation_iterations);
  fmt::print(stderr, FMT_STRING("Simulation took {}.\n"), absl::FormatDuration(absl::Now() - t0));

  if (!c.skip_output) {  // Mostly useful for profiling
    genome.write_contacts_to_file(c.output_dir, c.force);
    genome.write_extrusion_barriers_to_file(c.output_dir, c.force);
    std::ofstream cmd_file(fmt::format("{}/settings.log", c.output_dir));
    fmt::print(cmd_file, FMT_STRING("{}\n{}\n"), c.to_string(),
               absl::StrJoin(c.argv, c.argv + c.argc, " "));
  }
  fmt::print(stderr, "Simulation terminated without errors!\nBye.\n");
}
}  // namespace modle

int main(int argc, char** argv) noexcept {
  auto cli = modle::Cli(argc, argv);

  try {
    auto config = cli.parse_arguments();
    config.print();

    if (const auto files = check_if_output_file_exists(config); !files.empty()) {
      fmt::print(
          stderr,
          "Refusing to run the simulation because some of the output file(s) already exist. Pass "
          "--force to overwrite.\nCollision detected for the following file(s):\n - {}\n",
          absl::StrJoin(files.begin(), files.end(), "\n - "));
      return 1;
    }

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