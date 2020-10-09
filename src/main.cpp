#include <filesystem>
#include <string>

#include "absl/strings/str_format.h"
#include "absl/time/clock.h"
#include "modle/genome.hpp"

int main(int argc, char** argv) {
  if (argc < 4) {
    absl::FPrintF(stderr, "Usage: %s chromosomes.bed 20000 test/data/tmp_oputput/\n", argv[0]);
    return 1;
  }
  std::string_view path_to_bed = argv[1];
  uint32_t bin_size = std::stoul(argv[2]);
  std::string_view output_dir = argv[3];
  /*
  const uint32_t tot_n_barriers = 50'000;
  const uint32_t avg_lef_processivity = 25'000;
  const double probability_of_barrier_block = 0.95;
  const uint32_t tot_n_lefs = 5'000;
*/
  path_to_bed = "/home/roby/github/modle/test/data/test.bed";
  const uint32_t tot_n_barriers = 15;
  const uint32_t avg_lef_processivity = 100'000;
  const double probability_of_barrier_block = 0.95;
  const double probability_of_lef_rebind = 0.25;
  const double probability_of_extr_unit_bypass = 0.10;
      //  const double probability_of_barrier_block = 0.999975;
      const uint32_t tot_n_lefs = 45;

  auto t0 = absl::Now();
  modle::Genome genome(path_to_bed, bin_size, tot_n_lefs, avg_lef_processivity,
                       probability_of_barrier_block, probability_of_lef_rebind,
                       probability_of_extr_unit_bypass);
  genome.randomly_generate_barriers(tot_n_barriers);

  /*
  for (const auto& chr : genome.get_chromosomes()) {
    absl::FPrintF(stderr, "%s (%lu) has %lu barriers\n", chr.name, chr.length(), chr.n_barriers());
    for (const auto& barr : chr.barriers) {
      absl::FPrintF(stderr, "prob_of_block=%.5f\n", barr->get_prob_of_block());
    }
  }
   */

  /*
  uint32_t i = 0;
  for (const auto&chr: genome.get_chromosomes()) {
    for (const auto &b: chr.barriers) {
      absl::FPrintF(stderr, "Barrier %lu: %lu = %s\n", i++, b.get_abs_pos(), b. is_active() ?
  "Active" : "Inactive");
    }

  }

  sleep(5);
*/

  absl::FPrintF(stderr, "Genome size: %.2f Gbp; n_chr = %lu; n_bins = %lu\n", genome.size() / 1.0e9,
                genome.get_n_chromosomes(), genome.n_bins());
  absl::FPrintF(stderr, "Initial allocation took %s\n", absl::FormatDuration(absl::Now() - t0));

  // Associate extrusion barriers with DNA bins
  t0 = absl::Now();
  genome.randomly_bind_lefs();
  absl::FPrintF(stderr, "Initial lef binding took %s\n", absl::FormatDuration(absl::Now() - t0));

  genome.simulate_extrusion(25'000'000);

  std::filesystem::create_directories(output_dir);
  for (const auto& chr : genome.get_chromosomes()) {
    t0 = absl::Now();
    auto file = absl::StrFormat("%s/%s.tsv.gz", output_dir, chr.name);
    absl::FPrintF(stderr, "Writing full contact matrix for '%s' to file '%s'...", chr.name, file);
    chr.write_contacts_to_tsv(file, true);
    absl::FPrintF(stderr, " DONE in %s!\n", absl::FormatDuration(absl::Now() - t0));
    t0 = absl::Now();
    //    chr.write_contacts_to_tsv(file);
    file = absl::StrFormat("%s/%s_raw.tsv.gz", output_dir, chr.name);
    chr.write_contacts_to_tsv(file);
    absl::FPrintF(stderr, "Writing raw contact matrix for '%s' to file '%s'...", chr.name, file);
    absl::FPrintF(stderr, " DONE in %s!\n", absl::FormatDuration(absl::Now() - t0));
  }
  return 0;
}