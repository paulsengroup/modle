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

#include "modle/common.hpp"  // for bp_t, PRNG, seeder
#include "modle/config.hpp"  // for Config
#include "modle/dna.hpp"     // for Chromosome
#include "modle/utils.hpp"   // for ndebug_defined

namespace modle {

class ExtrusionBarrier;
class ExtrusionUnit;
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

  static constexpr auto NO_COLLISION = std::numeric_limits<collision_t>::max();
  static constexpr auto LEF_LEF_COLLISION = std::numeric_limits<collision_t>::max() - 1;
  static constexpr auto REACHED_CHROM_BOUNDARY = std::numeric_limits<collision_t>::max() - 2;

  [[nodiscard]] inline static constexpr bool is_lef_bar_collision(collision_t cause_of_collision);

  struct Task {
    std::size_t id;
    Chromosome* chrom;
    std::size_t cell_id;
    std::size_t n_target_epochs;
    std::size_t n_target_contacts;
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
    modle::PRNG_t rand_eng;
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
                        absl::Span<std::size_t> fwd_lef_ranks, MaskT& mask,
                        modle::PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

  inline void generate_ctcf_states(absl::Span<const ExtrusionBarrier> extr_barriers,
                                   boost::dynamic_bitset<>& mask,
                                   modle::PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

  [[nodiscard]] inline bp_t generate_rev_move(const Chromosome* const chrom,
                                              const ExtrusionUnit& unit, modle::PRNG_t& rand_eng);
  [[nodiscard]] inline bp_t generate_fwd_move(const Chromosome* const chrom,
                                              const ExtrusionUnit& unit, modle::PRNG_t& rand_eng);

  inline void generate_moves(const Chromosome* chrom, absl::Span<const Lef> lefs,
                             absl::Span<const std::size_t> rev_lef_ranks,
                             absl::Span<const std::size_t> fwd_lef_ranks,
                             absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
                             modle::PRNG_t& rand_eng,
                             bool adjust_moves_ = true) noexcept(utils::ndebug_defined());

  inline void adjust_moves(const Chromosome* chrom, absl::Span<const Lef> lefs,
                           absl::Span<const std::size_t> rev_lef_ranks,
                           absl::Span<const std::size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
                           absl::Span<bp_t> fwd_moves) noexcept(utils::ndebug_defined());

  inline static void rank_lefs(absl::Span<const Lef> lefs,
                               absl::Span<std::size_t> rev_lef_rank_buff,
                               absl::Span<std::size_t> fwd_lef_rank_buff,
                               bool init_buffers = false) noexcept(utils::ndebug_defined());

  inline void extrude(const Chromosome* chrom, absl::Span<Lef> lefs,
                      absl::Span<const bp_t> rev_moves,
                      absl::Span<const bp_t> fwd_moves) noexcept(utils::ndebug_defined());

  template <typename I>
  inline void process_lef_bar_collisions(absl::Span<const Lef> lefs,
                                         absl::Span<const std::size_t> rev_lef_ranks,
                                         absl::Span<const std::size_t> fwd_lef_ranks,
                                         absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
                                         absl::Span<const ExtrusionBarrier> extr_barriers,
                                         const boost::dynamic_bitset<>& barrier_mask,
                                         absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions,
                                         modle::PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

  template <typename I>
  inline void process_lef_lef_collisions(const Chromosome* chrom, absl::Span<const Lef> lefs,
                                         absl::Span<const ExtrusionBarrier> barriers,
                                         absl::Span<const std::size_t> rev_lef_ranks,
                                         absl::Span<const std::size_t> fwd_lef_ranks,
                                         absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
                                         absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions,
                                         PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

  template <typename I>
  inline void process_primary_lef_lef_collisions(
      const Chromosome* chrom, absl::Span<const Lef> lefs,
      absl::Span<const ExtrusionBarrier> barriers, absl::Span<const std::size_t> rev_lef_ranks,
      absl::Span<const std::size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
      absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions,
      PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

  template <typename I>
  inline void process_secondary_lef_lef_collisions(
      const Chromosome* chrom, absl::Span<const Lef> lefs,
      absl::Span<const std::size_t> rev_lef_ranks, absl::Span<const std::size_t> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions, PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

  inline std::size_t register_contacts(
      Chromosome* chrom, absl::Span<const Lef> lefs,
      absl::Span<const std::size_t> selected_lef_idx) noexcept(utils::ndebug_defined());

  template <typename MaskT>
  inline void select_lefs_to_bind(absl::Span<const Lef> lefs,
                                  MaskT& mask) noexcept(utils::ndebug_defined());

  inline void generate_lef_unloader_affinities(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<const collision_t> rev_collisions, absl::Span<const collision_t> fwd_collisions,
      absl::Span<double> lef_unloader_affinity) noexcept(utils::ndebug_defined());

  inline void select_lefs_to_release(absl::Span<std::size_t> lef_idx,
                                     absl::Span<const double> lef_unloader_affinity,
                                     modle::PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

  inline void release_lefs(absl::Span<Lef> lefs, absl::Span<const std::size_t> lef_idx) noexcept;

  [[nodiscard]] inline std::pair<bp_t /*rev*/, bp_t /*fwd*/> compute_lef_lef_collision_pos(
      const ExtrusionUnit& rev_unit, const ExtrusionUnit& fwd_unit, bp_t rev_move, bp_t fwd_move);

#ifdef ENABLE_TESTING
 public:
  template <typename MaskT>
  inline void test_bind_lefs(const Chromosome* chrom, absl::Span<Lef> lefs,
                             absl::Span<std::size_t> rev_lef_ranks,
                             absl::Span<std::size_t> fwd_lef_ranks, MaskT& mask,
                             modle::PRNG_t& rand_eng) {
    this->bind_lefs(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, mask, rand_eng);
  }

  inline void test_adjust_moves(const Chromosome* const chrom, absl::Span<const Lef> lefs,
                                absl::Span<const std::size_t> rev_lef_ranks,
                                absl::Span<const std::size_t> fwd_lef_ranks,
                                absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves) {
    this->adjust_moves(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves);
  }

  inline void test_generate_moves(const Chromosome* const chrom, absl::Span<const Lef> lefs,
                                  absl::Span<const std::size_t> rev_lef_ranks,
                                  absl::Span<const std::size_t> fwd_lef_ranks,
                                  absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
                                  modle::PRNG_t& rand_eng, bool adjust_moves_ = false) {
    this->generate_moves(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, rand_eng,
                         adjust_moves_);
  }

  inline static void test_rank_lefs(absl::Span<const Lef> lefs,
                                    absl::Span<std::size_t> rev_lef_rank_buff,
                                    absl::Span<std::size_t> fwd_lef_rank_buff, bool init_buffers) {
    Simulation::rank_lefs(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, init_buffers);
  }

  template <typename I>
  inline void test_process_lef_lef_collisions(
      const Chromosome* const chrom, absl::Span<const Lef> lefs,
      absl::Span<const ExtrusionBarrier> barriers, absl::Span<const std::size_t> rev_lef_rank_buff,
      absl::Span<const std::size_t> fwd_lef_rank_buff, absl::Span<bp_t> rev_move_buff,
      absl::Span<bp_t> fwd_move_buff, absl::Span<I> rev_collision_buff,
      absl::Span<I> fwd_collision_buff, modle::PRNG_t& rand_eng) {
    this->process_lef_lef_collisions(chrom, lefs, barriers, rev_lef_rank_buff, fwd_lef_rank_buff,
                                     rev_move_buff, fwd_move_buff, rev_collision_buff,
                                     fwd_collision_buff, rand_eng);
  }

  template <typename I>
  inline void test_process_lef_bar_collisions(
      absl::Span<const Lef> lefs, absl::Span<const std::size_t> rev_lef_rank_buff,
      absl::Span<const std::size_t> fwd_lef_rank_buff, absl::Span<bp_t> rev_move_buff,
      absl::Span<bp_t> fwd_move_buff, absl::Span<const ExtrusionBarrier> extr_barriers,
      const boost::dynamic_bitset<>& barrier_mask, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions, modle::PRNG_t& rand_eng) {
    this->process_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff,
                                     fwd_move_buff, extr_barriers, barrier_mask, rev_collisions,
                                     fwd_collisions, rand_eng);
  }
#endif
};
}  // namespace modle

#include "../../simulation_impl.hpp"  // IWYU pragma: keep
