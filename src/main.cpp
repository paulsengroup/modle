#include <random>
#include <string>

#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_format.h"
#include "absl/time/clock.h"
#include "modle/dna.hpp"
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

  const auto t0 = absl::Now();
  modle::Genome genome(path_to_bed, bin_size, tot_n_lefs, tot_n_barriers, avg_lef_processivity,
                       probability_of_barrier_block);
  absl::FPrintF(stderr, "Genome size: %.2f Gbp; n_chr = %lu; n_bins = %lu\n", genome.size() / 1.0e9,
                genome.n_chromosomes(), genome.n_bins());
  absl::FPrintF(stderr, "Duration: %s\n", absl::FormatDuration(absl::Now() - t0));

  // Associate extrusion barriers with DNA bins
  genome.randomly_bind_lefs();
  for (auto& lef : genome.get_lefs()) {
    if (lef.is_bound()) {
      const auto& [left, right] = lef.get_bins();
      auto left_pos = (left->get_end() + left->get_start()) / 2;
      auto right_pos = (right->get_end() + right->get_start()) / 2;
      absl::FPrintF(stderr, "chr=%s; left=%lu; right=%lu; loop size=%lu\n", lef.get_chr_name(),
                    left_pos, right_pos, right_pos - left_pos);
    } else {
      absl::FPrintF(stderr, "LEF not bound to DNA\n");
    }
  }
  absl::FPrintF(stderr, "%lu\n", genome.get_barriers().size());

  return 0;
}
