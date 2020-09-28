#ifndef MODLE_GENOME_HPP
#define MODLE_GENOME_HPP

#include <memory>

#include "modle/dna.hpp"
#include "modle/lefs.hpp"
#include "absl/container/flat_hash_map.h"
#include "absl/container/btree_map.h"

namespace modle {
class Genome {
 public:
  Genome(std::string_view path_to_bed, uint32_t bin_size, uint32_t n_lefs, uint32_t n_barriers,
         uint32_t avg_lef_processivity, double probability_of_barrier_block, uint64_t seed = 0);
  [[nodiscard]] const absl::flat_hash_map<const std::string, DNA>& get_chromosomes() const;
  [[nodiscard]] absl::flat_hash_map<const std::string, DNA>& get_chromosomes();
  [[nodiscard]] std::vector<Lef>& get_lefs();
  [[nodiscard]] const std::vector<Lef>& get_lefs() const;
  [[nodiscard]] absl::btree_multimap<std::string_view, ExtrusionBarrier>& get_barriers();
  [[nodiscard]] const absl::btree_multimap<std::string_view, ExtrusionBarrier>& get_barriers() const;
  [[nodiscard]] uint32_t n_chromosomes() const;
  [[nodiscard]] uint64_t size(uint32_t min_length = 0) const;
  [[nodiscard]] uint32_t n_bins() const;
  [[nodiscard]] uint32_t tot_n_lefs() const;
  [[nodiscard]] uint32_t n_barriers() const;

  void uniformly_bind_lefs(uint32_t min_length = 0);
  void randomly_bind_lefs(uint32_t min_length = 0);

 private:
  std::string _path_to_bed;
  uint32_t _bin_size;
  absl::flat_hash_map<const std::string, DNA> _chromosomes;
  std::default_random_engine _rng_dev;
  std::vector<Lef> _lefs;
  absl::flat_hash_map<std::string_view, std::shared_ptr<Lef>> _available_lefs{};
  uint32_t _avg_lef_processivity;
  absl::btree_multimap<std::string_view, ExtrusionBarrier> _barriers;
  double _probability_of_barrier_block;
  uint64_t _seed;

  absl::flat_hash_map<const std::string, DNA> init_chromosomes_from_bed();
  [[nodiscard]] static std::vector<Lef> generate_lefs(uint32_t n, uint32_t avg_processivity);
  [[nodiscard]] absl::btree_multimap<std::string_view, ExtrusionBarrier> generate_barriers(
      uint32_t n, double prob_of_block);
};
}  // namespace modle

#endif  // MODLE_GENOME_HPP
