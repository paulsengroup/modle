#pragma once

#include <absl/container/btree_set.h>

#include <cstddef>
#include <filesystem>
#include <string_view>
#include <vector>

#include "modle/dna.hpp"
#include "modle/extrusion_barriers.hpp"

namespace modle::io {

class Genome {
 public:
  inline Genome() = default;
  inline Genome(const std::filesystem::path& path_to_chrom_sizes,
                const std::filesystem::path& path_to_extr_barriers,
                const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms);

  [[nodiscard]] inline auto begin();
  [[nodiscard]] inline auto end();

  [[nodiscard]] inline auto begin() const;
  [[nodiscard]] inline auto end() const;

  [[nodiscard]] inline size_t size() const;
  [[nodiscard]] inline size_t simulated_size() const;

  [[nodiscard]] inline const Chromosome& chromosome_with_longest_name() const;
  [[nodiscard]] inline const Chromosome& longest_chromosome() const;
  [[nodiscard]] inline const Chromosome& chromosome_with_max_nbarriers() const;
  [[nodiscard]] inline size_t max_target_contacts(size_t bin_size, size_t diagonal_width,
                                                  double target_contact_density,
                                                  size_t simulation_iterations,
                                                  double lef_fraction_contact_sampling,
                                                  double nlefs_per_mbp, size_t ncells) const;

  /// Allocate and return the extrusion barriers for the chromosome \p chrom.
  [[nodiscard]] inline std::vector<ExtrusionBarrier> generate_vect_of_barriers(
      std::string_view chrom_name, double ctcf_prob_occ_to_occ, double ctcf_prob_nocc_to_nocc);

  /// A simple wrapper function that imports chromosomes and extrusion barriers that comprise the
  /// genome that is being simulated.
  [[nodiscard]] static inline absl::btree_set<Chromosome, Chromosome::Comparator>
  instantiate_genome(const std::filesystem::path& path_to_chrom_sizes,
                     const std::filesystem::path& path_to_extr_barriers,
                     const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms);

 private:
  absl::btree_set<Chromosome, Chromosome::Comparator> _chromosomes{};
};

}  // namespace modle::io

#include "../../libmodle_io_impl.hpp"
