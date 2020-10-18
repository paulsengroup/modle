#include <filesystem>
#include <string>

#include "absl/strings/str_format.h"
#include "absl/time/clock.h"
#include "modle/cli.hpp"
#include "modle/genome.hpp"
#include "modle/parsers.hpp"

namespace modle {

std::vector<std::string> check_if_output_file_exists(const modle::config& c) {
  std::vector<std::string> fn_collisions;
  if (c.force) return fn_collisions;
  modle::SimpleBEDParser parser(c.path_to_bed);
  modle::SimpleBED record;
  while (!(record = parser.parse_next()).empty()) {
    const auto f1 = absl::StrFormat("%s/%s.tsv.bz2", c.output_dir, record.chr);
    const auto f2 = absl::StrFormat("%s/%s_raw.tsv.bz2", c.output_dir, record.chr);
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

  genome.randomly_generate_barriers(c.number_of_barriers);

  absl::FPrintF(stderr, "Initialization took %s.\n", absl::FormatDuration(absl::Now() - t0));
  t0 = absl::Now();
  if (c.skip_burnin) {
    genome.randomly_bind_lefs();
    absl::FPrintF(stderr, "Bound %lu LEFs in %s.\n", genome.get_n_of_busy_lefs(),
                  absl::FormatDuration(absl::Now() - t0));
  } else {
    genome.run_burnin(c.probability_of_lef_rebind, c.min_n_of_loops_per_lef,
                      c.min_n_of_burnin_rounds);
  }

  t0 = absl::Now();
  genome.simulate_extrusion(c.simulation_iterations);
  absl::FPrintF(stderr, "Simulation took %s.\n", absl::FormatDuration(absl::Now() - t0));
  genome.write_contacts_to_file(c.output_dir.data(), c.force);

  if (c.make_heatmaps) genome.make_heatmaps(c.output_dir, c.force);
}
}  // namespace modle

int main(int argc, char** argv) {
  auto cli = modle::Cli(argc, argv);

  try {
    auto config = cli.parse_arguments();
    config.print();

    if (const auto files = check_if_output_file_exists(config); !files.empty()) {
      absl::FPrintF(
          stderr,
          "Refusing to run the simulation because some of the output files already exist. Pass "
          "--force to overwrite.\nCollision detected:\n - %s\n",
          absl::StrJoin(files.begin(), files.end(), "\n - "));
      return 1;
    }

    modle::run_simulation(config);
  } catch (const CLI::ParseError& e) {
    return cli.exit(e);  //  This takes care of formatting and printing error messages (if any)
  }
  return 0;
}