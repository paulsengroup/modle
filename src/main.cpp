#include <filesystem>
#include <string>

#include "absl/strings/str_format.h"
#include "absl/time/clock.h"
#include "modle/cli.hpp"
#include "modle/genome.hpp"

namespace modle {
void run_simulation(const modle::config& c) {
  auto t0 = absl::Now();
  modle::Genome genome(c.path_to_bed, c.bin_size, c.number_of_lefs, c.average_lef_processivity,
                       c.probability_of_barrier_block, c.probability_of_lef_rebind,
                       c.probability_of_extrusion_unit_bypass, c.seed);

  genome.randomly_generate_barriers(c.number_of_barriers);

  absl::FPrintF(stderr, "Initialization took %s.\n", absl::FormatDuration(absl::Now() - t0));
  //  genome.print_statistics();
  t0 = absl::Now();
  genome.randomly_bind_lefs();
  absl::FPrintF(stderr, "Bound %lu LEFs in %s.\n", genome.get_n_of_bound_lefs(),
                absl::FormatDuration(absl::Now() - t0));

  t0 = absl::Now();
  genome.simulate_extrusion(c.simulation_iterations);
  absl::FPrintF(stderr, "Simulation took %s.\nWriting contact matrices to directory '%s'...",
                absl::FormatDuration(absl::Now() - t0), c.output_dir);
  std::filesystem::create_directories(c.output_dir);
  for (const auto& chr : genome.get_chromosomes()) {
    t0 = absl::Now();
    auto file = absl::StrFormat("%s/%s.tsv.gz", c.output_dir, chr.name);
    absl::FPrintF(stderr, "Writing full contact matrix for '%s' to file '%s'...\n", chr.name, file);
    chr.write_contacts_to_tsv(file, true);
    file = absl::StrFormat("%s/%s_raw.tsv.gz", c.output_dir, chr.name);
    chr.write_contacts_to_tsv(file);
    absl::FPrintF(stderr, "Writing raw contact matrix for '%s' to file '%s'...\n", chr.name, file);
    absl::FPrintF(stderr, "DONE in %s!\n", absl::FormatDuration(absl::Now() - t0));
  }
  if (c.make_heatmaps) {
    t0 = absl::Now();
    absl::FPrintF(stderr, "Generating heatmaps...\n");
    absl::FPrintF(stderr, "TODO!!\n");
    absl::FPrintF(stderr, "DONE in %s!\n", absl::FormatDuration(absl::Now() - t0));
  }
}
}  // namespace modle

int main(int argc, char** argv) {
  auto cli = modle::Cli(argc, argv);

  try {
    auto config = cli.parse_arguments();
    modle::run_simulation(config);
  } catch (const CLI::ParseError& e) {
    return cli.exit(e);  //  This takes care of formatting and printing error messages (if any)
  }
  return 0;
}