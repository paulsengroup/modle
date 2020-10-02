#include <string>

#include "absl/strings/str_format.h"
#include "absl/time/clock.h"
#include "modle/genome.hpp"
#include <fstream>
#include "absl/strings/str_split.h"

void debugging_contacts() {
   uint32_t size = 5'000, diagonal_width = 50;
//  uint32_t size = 200, diagonal_width = 10;
  modle::ContactMatrix m1(diagonal_width, size);
  std::string path_to_file = "/home/roby/github/modle/test/data/symm_matrix_5000_50.tsv";
//  std::string path_to_file = "/home/roby/github/modle/test/data/symm_matrix_200_10.tsv";
  std::ifstream f(path_to_file);
  if (f.is_open()) {
    std::string line;
    std::string buff;
    uint32_t i = 0, j;
    while (std::getline(f, line)) {
      j = 0;
      for (const auto& tok : absl::StrSplit(line, '\t')) {
//        absl::FPrintF(stderr, "%lu;%lu='%s'\n", i, j, tok);
        if (tok != "0") {
          try {
            buff = tok;
            m1.set(i, j, std::stoul(buff));
          } catch (const std::invalid_argument &e) {
            absl::FPrintF(stderr, "Unable to convert '%s' to uint32_t: %s", buff, e.what());
            
          }
        }
        ++j;
      }
      ++i;
    }
  } else {
    throw std::runtime_error(absl::StrFormat("Unable to open file '%s'.", path_to_file));
  }
//  const auto m2 = m1.generate_symmetric_matrix();
  m1.print_symmetric_matrix();
}

int main(int argc, char** argv) {
  debugging_contacts();
  /*
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
  modle::Genome genome(path_to_bed, bin_size, tot_n_lefs, avg_lef_processivity,
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
  genome.simulate_extrusion(25000);
  for (const auto& lef : genome.get_lefs()) {
    absl::FPrintF(stderr, "Loop size: %lu\n", lef.get_loop_size());
    if (lef.right_is_stalled() || lef.left_is_stalled()) {
      absl::FPrintF(stderr, "chr=%s; loop_size=%lu; left_stalled=%s; right_stalled=%s;\n",
                    lef.get_chr_name(), lef.get_loop_size(),
                    lef.left_is_stalled() ? "True" : "False",
                    lef.right_is_stalled() ? "True" : "False");
    }
  }

  const auto& lefs = genome.get_lefs();
  absl::FPrintF(stderr, "Mean loop size: %.3f\n",
                std::accumulate(lefs.begin(), lefs.end(), 0UL,
                                [](uint32_t accumulator, const modle::Lef& lef) {
                                  return accumulator + lef.get_loop_size();
                                }) /
                    static_cast<double>(genome.n_lefs()));
*/




  return 0;
}