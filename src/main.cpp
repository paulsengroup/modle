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
  const uint32_t avg_lef_processivity = 350'000;
  const double probability_of_barrier_block = 0.99995;
//  const double probability_of_barrier_block = 0.999975;
  const uint32_t tot_n_lefs = 250;

  auto t0 = absl::Now();
  modle::Genome genome(path_to_bed, bin_size, tot_n_lefs, avg_lef_processivity,
                       probability_of_barrier_block);
  genome.randomly_generate_barriers(tot_n_barriers);

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

  genome.simulate_extrusion(100'000'000);
  for (const auto& lef : genome.get_lefs()) {
    absl::FPrintF(stderr, "Loop size: %lu\n", lef.get_loop_size());
    if (lef.right_is_stalled() || lef.left_is_stalled()) {
      absl::FPrintF(stderr, "chr=%s; loop_size=%lu; left_stalled=%s; right_stalled=%s;\n",
                    lef.get_chr_name(), lef.get_loop_size(),
                    lef.left_is_stalled() ? "True" : "False",
                    lef.right_is_stalled() ? "True" : "False");
    }
  }

  std::filesystem::create_directories(output_dir);
  for (const auto& chr : genome.get_chromosomes()) {
    t0 = absl::Now();
    const auto file = absl::StrFormat("%s/%s.tsv.gz", output_dir, chr.name);
    absl::FPrintF(stderr, "Writing contacts for '%s' to file '%s'...", chr.name, file);
    //    chr.write_contacts_to_tsv(file);
    chr.write_contacts_to_tsv(file, true);
    absl::FPrintF(stderr, " DONE in %s!\n", absl::FormatDuration(absl::Now() - t0));
  }

  const auto& lefs = genome.get_lefs();
  absl::FPrintF(stderr, "Mean loop size: %.3f\n",
                std::accumulate(lefs.begin(), lefs.end(), 0UL,
                                [](uint32_t accumulator, const modle::Lef& lef) {
                                  return accumulator + lef.get_loop_size();
                                }) /
                    static_cast<double>(genome.n_lefs()));

  for (const auto&chr: genome.get_chromosomes()) {
    uint64_t i = 0;
    for (const auto&b: chr.barriers) {
      absl::FPrintF(stderr, "b%lu=%lu\n", i++, b.get_abs_pos());
    }
  }

  return 0;
}