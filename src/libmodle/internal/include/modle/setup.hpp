#pragma once

#include <absl/container/btree_set.h>  // for btree_set, btree_set_container<>::const_iterator

#include <cstddef>      // for size_t
#include <cstdint>      // for uint32_t
#include <filesystem>   // for path
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "modle/dna.hpp"  // for Chromosome

namespace modle {

class ExtrusionBarrier;

class Genome {
 public:
  Genome() = default;
  Genome(const std::filesystem::path& path_to_chrom_sizes,
         const std::filesystem::path& path_to_extr_barriers,
         const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms);

  using iterator = absl::btree_set<Chromosome>::iterator;
  using const_iterator = absl::btree_set<Chromosome>::const_iterator;

  [[nodiscard]] iterator begin();
  [[nodiscard]] iterator end();

  [[nodiscard]] const_iterator begin() const;
  [[nodiscard]] const_iterator end() const;

  [[nodiscard]] size_t size() const;
  [[nodiscard]] size_t simulated_size() const;

  [[nodiscard]] const Chromosome& chromosome_with_longest_name() const;
  [[nodiscard]] const Chromosome& longest_chromosome() const;
  [[nodiscard]] const Chromosome& chromosome_with_max_nbarriers() const;
  [[nodiscard]] size_t max_target_contacts(size_t bin_size, size_t diagonal_width,
                                           double target_contact_density,
                                           size_t simulation_iterations,
                                           double lef_fraction_contact_sampling,
                                           double nlefs_per_mbp, size_t ncells) const;

  /// Allocate and return the extrusion barriers for the chromosome \p chrom.
  [[nodiscard]] std::vector<ExtrusionBarrier> generate_vect_of_barriers(
      std::string_view chrom_name, double ctcf_prob_occ_to_occ, double ctcf_prob_nocc_to_nocc);

  /// A simple wrapper function that imports chromosomes and extrusion barriers that comprise the
  /// genome that is being simulated.
  [[nodiscard]] static absl::btree_set<Chromosome> instantiate_genome(
      const std::filesystem::path& path_to_chrom_sizes,
      const std::filesystem::path& path_to_extr_barriers,
      const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms);

 private:
  absl::btree_set<Chromosome> _chromosomes{};
};

}  // namespace modle
