#pragma once

#include <absl/container/btree_set.h>
#include <absl/container/btree_set.h>  // IWYU pragma: keep for btree_set
#include <absl/types/span.h>
#include <absl/types/span.h>  // IWYU pragma: keep for Span

#include <boost/asio/thread_pool.hpp>               // for thread_pool
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cstddef>                                  // for size_t
#include <cstdint>                                  // for uint_fast16_t
#include <filesystem>                               // for path
#include <limits>                                   // for numeric_limits
#include <string>                                   // for string
#include <vector>                                   // for vector

#include "modle/common.hpp"  // for Bp
#include "modle/config.hpp"  // for Config
#include "modle/dna.hpp"     // for Chromosome

#ifdef USE_XOSHIRO
#include <Xoshiro-cpp/XoshiroCpp.hpp>  // for XoshiroCpp::Xoshiro256PlusPlus XoshiroCpp::SplitMix64
#else
#include <random>  // for geometric_distribution, mt19937_64
#endif

namespace modle {

#ifdef USE_XOSHIRO
using PRNG = XoshiroCpp::Xoshiro256PlusPlus;
using seeder = XoshiroCpp::SplitMix64;
#else
using PRNG = std::mt19937_64;
using seeder = std::seed_seq;
#endif

class ExtrusionBarrier;
struct Lef;

class Simulation : Config {
 public:
  inline explicit Simulation(const Config& c, bool import_chroms = true);
  inline std::string to_string() = delete;
  inline void print() = delete;

  [[nodiscard]] inline std::size_t size() const;
  [[nodiscard]] inline std::size_t simulated_size() const;

  using Genome = absl::btree_set<Chromosome, Chromosome::Comparator>;
  using lef_move_generator_t = std::normal_distribution<double>;
  using chrom_pos_generator_t = std::uniform_int_distribution<bp_t>;
  using collision_t = uint_fast16_t;
  static constexpr auto NO_COLLISION = std::numeric_limits<collision_t>::max();

  struct Task {
    std::size_t id;
    Chromosome* chrom;
    std::size_t cell_id;
    std::size_t n_target_epochs;
    std::size_t nlefs;
    absl::Span<const ExtrusionBarrier> barriers;
  };

  struct State : Task {
    std::vector<Lef> lef_buff;
    std::vector<double> lef_unloader_affinity;
    std::vector<std::size_t> rank_buff1;
    std::vector<std::size_t> rank_buff2;
    boost::dynamic_bitset<> barrier_mask;
    std::vector<bp_t> moves_buff1;
    std::vector<bp_t> moves_buff2;
    std::vector<std::size_t> idx_buff1;
    std::vector<std::size_t> idx_buff2;
    std::vector<std::size_t> epoch_buff;
    modle::PRNG rand_eng;
    uint64_t seed;

    inline void operator=(const Task& task);
    inline void resize(std::size_t size = std::numeric_limits<std::size_t>::max());
    inline void reset();
  };

  inline void run();

 private:
  Genome _chromosomes{};

  [[nodiscard]] inline boost::asio::thread_pool instantiate_thread_pool() const;
  template <typename I>
  [[nodiscard]] inline static boost::asio::thread_pool instantiate_thread_pool(
      I nthreads, bool clamp_nthreads = true);

  [[nodiscard]] inline static Genome import_chromosomes(
      const std::filesystem::path& path_to_chrom_sizes,
      const std::filesystem::path& path_to_extr_barriers,
      const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms);
  [[nodiscard]] inline std::vector<ExtrusionBarrier> allocate_barriers(const Chromosome* chrom);

  inline void simulate_extrusion_kernel(State& state);

  template <typename MaskT>
  inline void bind_lefs(const Chromosome* chrom, absl::Span<Lef> lefs,
                        absl::Span<std::size_t> rev_lef_ranks,
                        absl::Span<std::size_t> fwd_lef_ranks, modle::PRNG& rand_eng, MaskT& mask);

  inline void generate_ctcf_states(absl::Span<const ExtrusionBarrier> extr_barriers,
                                   boost::dynamic_bitset<>& mask, modle::PRNG& rand_eng);

  inline void generate_moves(const Chromosome* chrom, absl::Span<const Lef> lefs,
                             absl::Span<const std::size_t> rev_lef_ranks,
                             absl::Span<const std::size_t> fwd_lef_ranks,
                             absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
                             modle::PRNG& rand_eng);

  inline static void rank_lefs(absl::Span<const Lef> lefs,
                               absl::Span<std::size_t> rev_lef_rank_buff,
                               absl::Span<std::size_t> fwd_lef_rank_buff,
                               bool init_buffers = false);

  inline void extrude(const Chromosome* chrom, absl::Span<Lef> lefs,
                      absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
                      absl::Span<const std::size_t> rev_lef_mask,
                      absl::Span<const std::size_t> fwd_lef_mask);

  // Loop over lefs and identify colliding extr. units (i.e. units that travel in opposite
  // direction and that are within <p>dist_threshold</p> bp from each other
  inline void check_lef_lef_collisions(absl::Span<Lef> lefs,
                                       absl::Span<const std::size_t> rev_lef_rank_buff,
                                       absl::Span<const std::size_t> fwd_lef_rank_buff,
                                       absl::Span<const bp_t> rev_move_buff,
                                       absl::Span<const bp_t> fwd_move_buff,
                                       absl::Span<std::size_t> rev_collision_buff,
                                       absl::Span<std::size_t> fwd_collision_buff, PRNG& rand_eng);

  inline void check_lef_bar_collisions(
      absl::Span<Lef> lefs, absl::Span<const std::size_t> rev_lef_rank_buff,
      absl::Span<const std::size_t> fwd_lef_rank_buff, absl::Span<const bp_t> rev_move_buff,
      absl::Span<const bp_t> fwd_move_buff, absl::Span<const ExtrusionBarrier> extr_barriers,
      const boost::dynamic_bitset<>& barrier_mask, absl::Span<std::size_t> rev_collisions,
      absl::Span<std::size_t> fwd_collisions);

  inline void register_contacts(Chromosome* chrom, absl::Span<const Lef> lefs,
                                absl::Span<const std::size_t> selected_lef_idx);

  template <typename MaskT>
  inline void select_lefs_to_bind(absl::Span<const Lef> lefs, MaskT& mask);

  inline void generate_lef_unloader_affinities(absl::Span<const Lef> lefs,
                                               absl::Span<const ExtrusionBarrier> barriers,
                                               absl::Span<const std::size_t> rev_collisions,
                                               absl::Span<const std::size_t> fwd_collisions,
                                               absl::Span<double> lef_unloader_affinity);

  inline void select_lefs_to_release(absl::Span<std::size_t> lef_idx,
                                     absl::Span<const double> lef_unloader_affinity,
                                     modle::PRNG& rand_eng);

  inline void release_lefs(absl::Span<Lef> lefs, absl::Span<const std::size_t> lef_idx);

#ifdef ENABLE_TESTING
 public:
  inline void test_check_lef_lef_collisions(absl::Span<const Lef> lefs,
                                            absl::Span<const std::size_t> rev_lef_rank_buff,
                                            absl::Span<const std::size_t> fwd_lef_rank_buff,
                                            absl::Span<const bp_t> rev_move_buff,
                                            absl::Span<const bp_t> fwd_move_buff,
                                            absl::Span<collision_t> rev_collision_buff,
                                            absl::Span<collision_t> fwd_collision_buff) {
    this->check_lef_lef_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff,
                                   fwd_move_buff, rev_collision_buff, fwd_collision_buff);
  }

  inline void test_apply_lef_lef_stalls(absl::Span<Lef> lefs,
                                        absl::Span<const collision_t> rev_collision_buff,
                                        absl::Span<const collision_t> fwd_collision_buff,
                                        absl::Span<const std::size_t> rev_lef_rank_buff,
                                        absl::Span<const std::size_t> fwd_lef_rank_buff,
                                        PRNG& rand_eng) {
    this->apply_lef_lef_stalls(lefs, rev_collision_buff, fwd_collision_buff, rev_lef_rank_buff,
                               fwd_lef_rank_buff, rand_eng);
  }

  inline void test_check_lef_bar_collisions(
      absl::Span<const Lef> lefs, absl::Span<const std::size_t> rev_lef_rank_buff,
      absl::Span<const std::size_t> fwd_lef_rank_buff, absl::Span<const bp_t> rev_move_buff,
      absl::Span<const bp_t> fwd_move_buff, absl::Span<const ExtrusionBarrier> extr_barriers,
      absl::Span<collision_t> rev_collision_buff, absl::Span<collision_t> fwd_collision_buff) {
    this->check_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff,
                                   fwd_move_buff, extr_barriers, rev_collision_buff,
                                   fwd_collision_buff);
  }

  inline void test_apply_lef_bar_stalls(absl::Span<Lef> lefs,
                                        absl::Span<const collision_t> rev_collision_buff,
                                        absl::Span<const collision_t> fwd_collision_buff,
                                        absl::Span<const ExtrusionBarrier> extr_barriers,
                                        PRNG& rand_eng) {
    this->apply_lef_bar_stalls(lefs, rev_collision_buff, fwd_collision_buff, extr_barriers,
                               rand_eng);
  }
#endif
};
}  // namespace modle

#include "../../simulation_impl.hpp"  // IWYU pragma: keep
