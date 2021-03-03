#pragma once

#include <absl/container/btree_set.h>
#include <absl/types/span.h>

#include <boost/asio/thread_pool.hpp>
#include <boost/dynamic_bitset.hpp>
#include <cstdint>      // for uint*_t
#include <memory>       // for _Destroy, allocator
#include <string>       // for basic_string, string
#include <string_view>  // for string_view
#include <utility>      // for pair
#include <vector>       // for vector

#include "modle/bed.hpp"
#include "modle/dna.hpp"
#include "modle/extrusion_factors.hpp"

namespace modle {

struct config;

class Genome {
 public:
  inline void test();
  inline explicit Genome(const config& c, bool import_chroms = true);

  [[nodiscard]] std::size_t size() const;
  [[nodiscard]] std::size_t simulated_size() const;

  using Chromosomes = absl::btree_set<Chromosome, Chromosome::Comparator>;
  using lef_bar_stall_generator_t = std::geometric_distribution<Bp>;
  using lef_lef_stall_generator_t = std::geometric_distribution<Bp>;
  using lef_lifetime_generator_t = std::geometric_distribution<Bp>;
  using chrom_pos_generator_t = std::uniform_int_distribution<Bp>;
  void simulate_extrusion(std::size_t iterations);

 private:
  std::filesystem::path _path_to_chrom_sizes;
  std::filesystem::path _path_to_chrom_subranges;
  std::filesystem::path _path_to_extr_barriers;
  uint32_t _bin_size;
  uint32_t _avg_lef_lifetime;  ///< Average loop size if loop extrusion takes place unobstructed
  double _probability_of_barrier_block;
  double _probability_of_lef_rebind;
  double _probability_of_extr_unit_bypass;
  double _soft_stall_multiplier;
  double _hard_stall_multiplier;
  bool _allow_lef_lifetime_extension;
  uint32_t _sampling_interval;
  bool _randomize_contact_sampling;
  uint32_t _nthreads;
  uint64_t _seed;

  Chromosomes _chromosomes{};

  // inline void simulate_extrusion();

  // inline void simulate_extrusion(double target_contact_density);
  // inline void simulate_extrusion(uint32_t iterations, double target_contact_density);
  [[nodiscard]] inline boost::asio::thread_pool instantiate_thread_pool() const;
  template <typename I>
  [[nodiscard]] inline static boost::asio::thread_pool instantiate_thread_pool(I nthreads);

  [[nodiscard]] inline static Chromosomes import_chromosomes(
      const std::filesystem::path& path_to_chrom_sizes,
      const std::filesystem::path& path_to_extr_barriers,
      const std::filesystem::path& path_to_chrom_subranges = {});

  inline void simulate_extrusion_kernel(const Chromosome* chrom, std::size_t ncell,
                                        std::vector<Lef> lefs, std::vector<Bp> fwd_barriers,
                                        std::vector<Bp> rev_barriers,
                                        std::vector<double> fwd_barriers_pblock,
                                        std::vector<double> rev_barriers_pblock);
  template <typename I>
  inline void bind_lefs(const Chromosome* chrom, std::vector<Lef>& lefs, modle::PRNG& rand_eng,
                        const std::vector<I>& mask);

  inline void bind_all_lefs(const Chromosome* chrom, std::vector<Lef>& lefs, modle::PRNG& rand_eng);

  // Loop over lefs and identifi colliding extr. units (i.e. units that travel in opposite direction
  // and that are within <p>dist_threshold</p> bp from each other
  template <typename I>
  inline static void check_lef_lef_collisions(const std::vector<Lef>& lefs,
                                              const std::vector<Bp>& fwd_rank_buff,
                                              const std::vector<Bp>& rev_rank_buff,
                                              std::vector<I>& fwd_collision_buff,
                                              std::vector<I>& rev_collision_buff,
                                              Bp dist_threshold);
  template <typename I>
  inline void check_lef_lef_collisions(const std::vector<Lef>& lefs,
                                       const std::vector<Bp>& fwd_rank_buff,
                                       const std::vector<Bp>& rev_rank_buff,
                                       std::vector<I>& fwd_collision_buff,
                                       std::vector<I>& rev_collision_buff);

#ifdef ENABLE_TESTING
 public:
  template <typename I>
  inline void test_check_lef_lef_collisions(const std::vector<Lef>& lefs,
                                            const std::vector<Bp>& fwd_rank_buff,
                                            const std::vector<Bp>& rev_rank_buff,
                                            std::vector<I>& fwd_collision_buff,
                                            std::vector<I>& rev_collision_buff) {
    this->check_lef_lef_collisions(lefs, fwd_rank_buff, rev_rank_buff, fwd_collision_buff,
                                   rev_collision_buff);
  }

#endif
};
}  // namespace modle

#include "../../genome_impl.hpp"
