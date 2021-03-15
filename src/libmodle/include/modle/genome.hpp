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
#include <thread>
#include <utility>  // for pair
#include <vector>   // for vector

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
  inline explicit Genome(const config& c, bool import_chroms = true);

  [[nodiscard]] inline std::size_t size() const;
  [[nodiscard]] inline std::size_t simulated_size() const;

  using Chromosomes = absl::btree_set<Chromosome, Chromosome::Comparator>;
  using lef_bar_stall_generator_t = std::geometric_distribution<Bp>;
  using lef_lef_stall_generator_t = std::geometric_distribution<Bp>;
  using lef_lifetime_generator_t = std::geometric_distribution<Bp>;
  using chrom_pos_generator_t = std::uniform_int_distribution<Bp>;
  using collision_t = uint_fast16_t;
  inline void simulate_extrusion(const std::filesystem::path& output_path, uint32_t ncells,
                                 uint32_t simulation_rounds);
  inline void simulate_extrusion(const std::filesystem::path& output_path, uint32_t ncells,
                                 double target_contact_density);

  struct Task {
    std::size_t id;
    Chromosome* chrom;
    std::size_t cell_id;
    std::size_t nrounds;
    std::size_t nlefs;
    std::shared_ptr<std::vector<ExtrusionBarrier>> barriers;
  };

 private:
  std::filesystem::path _path_to_chrom_sizes;
  std::filesystem::path _path_to_chrom_subranges;
  std::filesystem::path _path_to_extr_barriers;
  Bp _bin_size;
  Bp _diagonal_width;
  Bp _avg_lef_lifetime;  ///< Average loop size if loop extrusion takes place unobstructed
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

  [[nodiscard]] inline boost::asio::thread_pool instantiate_thread_pool() const;
  template <typename I>
  [[nodiscard]] inline static boost::asio::thread_pool instantiate_thread_pool(
      I nthreads, bool clamp_nthreads = true);

  [[nodiscard]] inline static Chromosomes import_chromosomes(
      const std::filesystem::path& path_to_chrom_sizes,
      const std::filesystem::path& path_to_extr_barriers,
      const std::filesystem::path& path_to_chrom_subranges = {});
  [[nodiscard]] inline static std::vector<ExtrusionBarrier> allocate_barriers(
      const Chromosome* chrom, double default_prob_of_block);

  inline void simulate_extrusion(const std::filesystem::path& output_path, uint32_t ncells,
                                 uint32_t simulation_rounds, double target_contact_density);

  inline void simulate_extrusion_kernel(
      Chromosome* chrom, std::size_t cell_id, std::size_t simulation_rounds,
      std::vector<Lef> lef_buff,
      std::shared_ptr<const std::vector<ExtrusionBarrier>> extr_barrier_buff);
  template <typename MaskT>
  inline void bind_lefs(const Chromosome* chrom, absl::Span<Lef> lefs,
                        absl::Span<std::size_t> rev_lef_rank_buff,
                        absl::Span<std::size_t> fwd_lef_rank_buff, modle::PRNG& rand_eng,
                        MaskT& mask);

  inline void bind_all_lefs(const Chromosome* chrom, absl::Span<Lef> lefs,
                            absl::Span<std::size_t> rev_lef_rank_buff,
                            absl::Span<std::size_t> fwd_lef_rank_buff, PRNG& rand_eng);

  inline static void rank_lefs(absl::Span<Lef> lefs, absl::Span<std::size_t> rev_lef_rank_buff,
                               absl::Span<std::size_t> fwd_lef_rank_buff,
                               bool init_buffers = false);

  inline void extrude(const Chromosome* chrom, absl::Span<Lef> lefs);

  // Loop over lefs and identify colliding extr. units (i.e. units that travel in opposite
  // direction and that are within <p>dist_threshold</p> bp from each other
  inline static void check_lef_lef_collisions(absl::Span<const Lef> lefs,
                                              absl::Span<const std::size_t> rev_lef_rank_buff,
                                              absl::Span<const std::size_t> fwd_lef_rank_buff,
                                              absl::Span<collision_t> rev_collision_buff,
                                              absl::Span<collision_t> fwd_collision_buff,
                                              Bp dist_threshold);
  inline void check_lef_lef_collisions(absl::Span<const Lef> lefs,
                                       absl::Span<const std::size_t> rev_lef_rank_buff,
                                       absl::Span<const std::size_t> fwd_lef_rank_buff,
                                       absl::Span<collision_t> rev_collision_buff,
                                       absl::Span<collision_t> fwd_collision_buff);

  inline static void apply_lef_lef_stalls(absl::Span<Lef> lefs,
                                          absl::Span<const collision_t> rev_collision_buff,
                                          absl::Span<const collision_t> fwd_collision_buff,
                                          absl::Span<const std::size_t> rev_lef_rank_buff,
                                          absl::Span<const std::size_t> fwd_lef_rank_buff,
                                          PRNG& rand_eng, double prob_of_bypass);

  inline static void check_lef_bar_collisions(absl::Span<const Lef> lefs,
                                              absl::Span<const std::size_t> rev_lef_rank_buff,
                                              absl::Span<const std::size_t> fwd_lef_rank_buff,
                                              absl::Span<const ExtrusionBarrier> extr_barriers,
                                              absl::Span<collision_t> rev_collision_buff,
                                              absl::Span<collision_t> fwd_collision_buff,
                                              Bp dist_threshold);

  inline void check_lef_bar_collisions(absl::Span<const Lef> lefs,
                                       absl::Span<const std::size_t> rev_lef_rank_buff,
                                       absl::Span<const std::size_t> fwd_lef_rank_buff,
                                       absl::Span<const ExtrusionBarrier> extr_barriers,
                                       absl::Span<collision_t> rev_collision_buff,
                                       absl::Span<collision_t> fwd_collision_buff);

  inline void apply_lef_bar_stalls(absl::Span<Lef> lefs,
                                   absl::Span<const collision_t> rev_collision_buff,
                                   absl::Span<const collision_t> fwd_collision_buff,
                                   absl::Span<const ExtrusionBarrier> extr_barriers, PRNG& rand_eng,
                                   double soft_stall_multiplier, double hard_stall_multiplier);

  inline void register_contacts(Chromosome* chrom, absl::Span<const Lef> lefs);

  template <typename MaskT>
  inline void select_lefs_to_bind(absl::Span<const Lef> lefs, MaskT& mask);

#ifdef ENABLE_TESTING
 public:
  inline void test_check_lef_lef_collisions(absl::Span<const Lef> lefs,
                                            absl::Span<const std::size_t> rev_lef_rank_buff,
                                            absl::Span<const std::size_t> fwd_lef_rank_buff,
                                            absl::Span<collision_t> rev_collision_buff,
                                            absl::Span<collision_t> fwd_collision_buff) {
    this->check_lef_lef_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_collision_buff,
                                   fwd_collision_buff);
  }

  inline void test_apply_lef_lef_stalls(absl::Span<Lef> lefs,
                                        const absl::Span<collision_t> rev_collision_buff,
                                        const absl::Span<collision_t> fwd_collision_buff,
                                        absl::Span<const std::size_t> rev_lef_rank_buff,
                                        absl::Span<const std::size_t> fwd_lef_rank_buff,
                                        PRNG& rand_eng, double prob_of_block) {
    this->apply_lef_lef_stalls(lefs, rev_collision_buff, fwd_collision_buff, rev_lef_rank_buff,
                               fwd_lef_rank_buff, rand_eng, prob_of_block);
  }

  inline void test_check_lef_bar_collisions(absl::Span<const Lef> lefs,
                                            absl::Span<const std::size_t> rev_lef_rank_buff,
                                            absl::Span<const std::size_t> fwd_lef_rank_buff,
                                            absl::Span<const ExtrusionBarrier> extr_barriers,
                                            absl::Span<collision_t> rev_collision_buff,
                                            absl::Span<collision_t> fwd_collision_buff) {
    this->check_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, extr_barriers,
                                   rev_collision_buff, fwd_collision_buff);
  }

  inline void test_apply_lef_bar_stalls(absl::Span<Lef> lefs,
                                        const absl::Span<collision_t> rev_collision_buff,
                                        const absl::Span<collision_t> fwd_collision_buff,
                                        absl::Span<const ExtrusionBarrier> extr_barriers,
                                        PRNG& rand_eng, double soft_stall_multiplier,
                                        double hard_stall_multiplier) {
    this->apply_lef_bar_stalls(lefs, rev_collision_buff, fwd_collision_buff, extr_barriers,
                               rand_eng, soft_stall_multiplier, hard_stall_multiplier);
  }

#endif
};
}  // namespace modle

#include "../../genome_impl.hpp"
