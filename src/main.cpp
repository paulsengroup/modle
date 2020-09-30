#include <random>
#include <string>

#include "absl/strings/str_format.h"
#include "absl/time/clock.h"
#include "modle/genome.hpp"

int main(int argc, char** argv) {
  if (argc < 3) {
    absl::FPrintF(stderr, "Usage: %s chromosomes.bed 20000\n", argv[0]);
    return 1;
  }
  std::string_view path_to_bed = argv[1];
  uint32_t bin_size = std::stoul(argv[2]);
  const uint32_t tot_n_barriers = 50'000;
  const uint32_t avg_lef_processivity = 25'000;
  const double probability_of_barrier_block = 0.95;
  const uint32_t tot_n_lefs = 5'000;

  auto t0 = absl::Now();
  modle::Genome genome(path_to_bed, bin_size, tot_n_lefs, tot_n_barriers, avg_lef_processivity,
                       probability_of_barrier_block);
  genome.randomly_generate_barriers(tot_n_barriers, probability_of_barrier_block);
  absl::FPrintF(stderr, "Genome size: %.2f Gbp; n_chr = %lu; n_bins = %lu\n", genome.size() / 1.0e9,
                genome.get_n_chromosomes(), genome.n_bins());
  absl::FPrintF(stderr, "Initial allocation took %s\n", absl::FormatDuration(absl::Now() - t0));

  // Associate extrusion barriers with DNA bins
  t0 = absl::Now();
  genome.randomly_bind_lefs();
  absl::FPrintF(stderr, "Initial lef binding took %s\n", absl::FormatDuration(absl::Now() - t0));

  //  for (uint32_t i = 1; i <= 5'000; ++i) {
  genome.simulate_extrusion(5000);
  for (const auto& lef : genome.get_lefs()) {
    if (lef.right_is_stalled() && lef.left_is_stalled()) {
      absl::FPrintF(stderr, "chr=%s; loop_size=%lu; left_stalled=%s; right_stalled=%s;\n",
                    lef.get_chr_name(), lef.get_loop_size(),
                    lef.left_is_stalled() ? "True" : "False",
                    lef.right_is_stalled() ? "True" : "False");
    }
  }

  return 0;
}
