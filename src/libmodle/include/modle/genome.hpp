#pragma once

#include <absl/container/btree_set.h>
#include <absl/types/span.h>

#ifdef USE_XOSHIRO
#include <XoshiroCpp.hpp>
#else
#include <random>
#endif

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
#include "modle/extrusion_barriers.hpp"
#include "modle/extrusion_factors.hpp"

namespace modle {

#ifdef USE_XOSHIRO
using PRNG = XoshiroCpp::Xoshiro256PlusPlus;
using seeder = XoshiroCpp::SplitMix64;
#else
using PRNG = std::mt19937_64;
using seeder = std::seed_seq;
#endif

struct config;

class Genome {
 public:
  inline void test();
  inline explicit Genome(const config& c, bool import_chroms = true);

  [[nodiscard]] inline std::size_t size() const;
  [[nodiscard]] inline std::size_t simulated_size() const;

  using Chromosomes = absl::btree_set<Chromosome, Chromosome::Comparator>;
  using lef_bar_stall_generator_t = std::geometric_distribution<Bp>;
  using lef_lef_stall_generator_t = std::geometric_distribution<Bp>;
  using lef_lifetime_generator_t = std::geometric_distribution<Bp>;
  using chrom_pos_generator_t = std::uniform_int_distribution<Bp>;
  inline void simulate_extrusion(std::size_t iterations);

 private:
  std::filesystem::path _path_to_chrom_sizes;
  std::filesystem::path _path_to_chrom_subranges;
  std::filesystem::path _path_to_extr_barriers;
  uint32_t _bin_size;
  uint32_t _avg_lef_lifetime;  ///< Average loop size if loop extrusion takes place unobstructed
  double _nlefs_per_mbp;
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

  inline void simulate_extrusion_kernel(const Chromosome* chrom, std::size_t cell_id,
                                        std::size_t burnin_iters, std::size_t simulation_iters,
                                        std::vector<Lef> lefs,
                                        std::vector<ExtrusionBarrier> extr_barriers);
  template <typename I>
  inline void bind_lefs(const Chromosome* chrom, std::vector<Lef>& lefs, modle::PRNG& rand_eng,
                        const std::vector<I>& mask);

  inline void bind_all_lefs(const Chromosome* chrom, std::vector<Lef>& lefs,
                            std::vector<std::size_t>& fwd_lef_rank_buff,
                            std::vector<std::size_t>& rev_lef_rank_buff, PRNG& rand_eng);

  inline static void rank_lefs(std::vector<Lef>& lefs, std::vector<std::size_t>& fwd_lef_rank_buff,
                               std::vector<std::size_t>& rev_lef_rank_buff,
                               bool init_buffers = false);

  inline void extrude(const Chromosome* chrom, std::vector<Lef>& lefs);

  // Loop over lefs and identify colliding extr. units (i.e. units that travel in opposite
  // direction and that are within <p>dist_threshold</p> bp from each other
  template <typename I>
  inline static void check_lef_lef_collisions(const std::vector<Lef>& lefs,
                                              const std::vector<std::size_t>& rev_lef_rank_buff,
                                              const std::vector<std::size_t>& fwd_lef_rank_buff,
                                              std::vector<I>& rev_collision_buff,
                                              std::vector<I>& fwd_collision_buff,
                                              Bp dist_threshold);
  template <typename I>
  inline void check_lef_lef_collisions(const std::vector<Lef>& lefs,
                                       const std::vector<std::size_t>& rev_lef_rank_buff,
                                       const std::vector<std::size_t>& fwd_lef_rank_buff,
                                       std::vector<I>& rev_collision_buff,
                                       std::vector<I>& fwd_collision_buff);

  template <typename I>
  inline static void apply_lef_lef_stalls(std::vector<Lef>& lefs,
                                          const std::vector<I>& rev_collision_buff,
                                          const std::vector<I>& fwd_collision_buff,
                                          const std::vector<std::size_t>& rev_lef_rank_buff,
                                          const std::vector<std::size_t>& fwd_lef_rank_buff,
                                          PRNG& rand_eng, double prob_of_bypass);

  template <typename I>
  inline static void check_lef_bar_collisions(const std::vector<Lef>& lefs,
                                              const std::vector<std::size_t>& rev_lef_rank_buff,
                                              const std::vector<std::size_t>& fwd_lef_rank_buff,
                                              const std::vector<ExtrusionBarrier>& extr_barriers,
                                              std::vector<I>& rev_collision_buff,
                                              std::vector<I>& fwd_collision_buff,
                                              Bp dist_threshold);

  template <typename I>
  inline void check_lef_bar_collisions(const std::vector<Lef>& lefs,
                                       const std::vector<std::size_t>& rev_lef_rank_buff,
                                       const std::vector<std::size_t>& fwd_lef_rank_buff,
                                       const std::vector<ExtrusionBarrier>& extr_barriers,
                                       std::vector<I>& rev_collision_buff,
                                       std::vector<I>& fwd_collision_buff);
  template <typename I>
  inline void apply_lef_bar_stalls(std::vector<Lef>& lefs, const std::vector<I>& rev_collision_buff,
                                   const std::vector<I>& fwd_collision_buff,
                                   const std::vector<ExtrusionBarrier>& extr_barriers,
                                   const std::vector<std::size_t>& rev_lef_rank_buff,
                                   const std::vector<std::size_t>& fwd_lef_rank_buff,
                                   PRNG& rand_eng, double soft_stall_multiplier,
                                   double hard_stall_multiplier);

#ifdef ENABLE_TESTING
 public:
  template <typename I>
  inline void test_check_lef_lef_collisions(const std::vector<Lef>& lefs,
                                            const std::vector<std::size_t>& rev_lef_rank_buff,
                                            const std::vector<std::size_t>& fwd_lef_rank_buff,
                                            std::vector<I>& rev_collision_buff,
                                            std::vector<I>& fwd_collision_buff) {
    this->check_lef_lef_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_collision_buff,
                                   fwd_collision_buff);
  }

  template <typename I>
  inline void test_apply_lef_lef_stalls(std::vector<Lef>& lefs,
                                        const std::vector<I>& rev_collision_buff,
                                        const std::vector<I>& fwd_collision_buff,
                                        const std::vector<std::size_t>& rev_lef_rank_buff,
                                        const std::vector<std::size_t>& fwd_lef_rank_buff,
                                        PRNG& rand_eng, double prob_of_block) {
    this->apply_lef_lef_stalls(lefs, rev_collision_buff, fwd_collision_buff, rev_lef_rank_buff,
                               fwd_lef_rank_buff, rand_eng, prob_of_block);
  }

  template <typename I>
  inline void test_check_lef_bar_collisions(const std::vector<Lef>& lefs,
                                            const std::vector<std::size_t>& rev_lef_rank_buff,
                                            const std::vector<std::size_t>& fwd_lef_rank_buff,
                                            const std::vector<ExtrusionBarrier>& extr_barriers,
                                            std::vector<I>& rev_collision_buff,
                                            std::vector<I>& fwd_collision_buff) {
    this->check_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, extr_barriers,
                                   rev_collision_buff, fwd_collision_buff);
  }

  template <typename I>
  inline void test_apply_lef_bar_stalls(std::vector<Lef>& lefs,
                                        const std::vector<I>& rev_collision_buff,
                                        const std::vector<I>& fwd_collision_buff,
                                        const std::vector<ExtrusionBarrier>& extr_barriers,
                                        const std::vector<std::size_t>& rev_lef_rank_buff,
                                        const std::vector<std::size_t>& fwd_lef_rank_buff,
                                        PRNG& rand_eng, double soft_stall_multiplier,
                                        double hard_stall_multiplier) {
    this->template apply_lef_bar_stalls(lefs, rev_collision_buff, fwd_collision_buff, extr_barriers,
                                        rev_lef_rank_buff, fwd_lef_rank_buff, rand_eng,
                                        soft_stall_multiplier, hard_stall_multiplier);
  }

#endif
};
}  // namespace modle

#include "../../genome_impl.hpp"
