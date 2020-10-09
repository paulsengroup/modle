#ifndef MODLE_GENOME_HPP
#define MODLE_GENOME_HPP

#include <memory>

#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_map.h"
#include "modle/contacts.hpp"
#include "modle/dna.hpp"
#include "modle/extr_barrier.hpp"
#include "modle/lefs.hpp"

namespace modle {

class Genome {
 public:
  Genome(std::string_view path_to_bed, uint32_t bin_size, uint32_t n_lefs,
         uint32_t avg_lef_processivity, double probability_of_barrier_block,
         double probability_of_lef_rebind, double probability_of_extr_unit_bypass,
         uint64_t seed = 0);
  [[nodiscard]] const std::vector<Chromosome>& get_chromosomes() const;
  [[nodiscard]] std::vector<Chromosome>& get_chromosomes();
  [[nodiscard]] std::vector<Lef>& get_lefs();
  [[nodiscard]] const std::vector<Lef>& get_lefs() const;
  [[nodiscard]] uint32_t get_n_chromosomes() const;
  [[nodiscard]] uint64_t size() const;
  [[nodiscard]] uint32_t n_bins() const;
  [[nodiscard]] uint32_t n_lefs() const;
  [[nodiscard]] uint32_t n_barriers() const;
  [[nodiscard]] std::vector<uint32_t> get_chromosome_lengths() const;
  [[nodiscard]] uint32_t get_n_of_free_lefs() const;
  [[nodiscard]] uint32_t get_n_of_bound_lefs() const;
  [[nodiscard]] uint32_t get_n_lefs() const;

  void randomly_generate_barriers(uint32_t n_barriers);
  void randomly_bind_lefs();
  void simulate_extrusion(uint32_t iterations = 1);

 private:
  std::string _path_to_bed;
  uint32_t _bin_size;
  uint32_t _avg_lef_processivity;
  double _probability_of_barrier_block;
  double _probability_of_lef_rebind;
  double _probability_of_extr_unit_bypass;
  std::vector<Lef> _lefs;
  std::vector<Chromosome> _chromosomes;
  uint64_t _seed;
  std::mt19937 _rand_gen;

  [[nodiscard]] std::vector<Chromosome> init_chromosomes_from_bed() const;
  [[nodiscard]] std::vector<Lef> generate_lefs(uint32_t n);
};
}  // namespace modle

#endif  // MODLE_GENOME_HPP
