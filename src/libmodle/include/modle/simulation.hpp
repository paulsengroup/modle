#pragma once

#include <absl/container/btree_set.h>
#include <absl/container/btree_set.h>  // IWYU pragma: keep for btree_set
#include <absl/types/span.h>
#include <absl/types/span.h>  // IWYU pragma: keep for Span

#include <boost/asio/thread_pool.hpp>               // for thread_pool
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cstddef>                                  // for size_t
#include <cstdint>                                  // for uint64_t
#include <filesystem>                               // for path
#include <limits>                                   // for numeric_limits
#include <string>                                   // for string
#include <utility>                                  // for pair
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

  struct Task {
    inline Task() = default;
    std::size_t id{};
    Chromosome* chrom{};
    std::size_t cell_id{};
    std::size_t n_target_epochs{};
    std::size_t n_target_contacts{};
    std::size_t nlefs{};
    absl::Span<const ExtrusionBarrier> barriers{};
  };

  struct State : Task {
    inline State() = default;
    std::vector<Lef> lef_buff{};
    std::vector<double> lef_unloader_affinity{};
    std::vector<std::size_t> rank_buff1{};
    std::vector<std::size_t> rank_buff2{};
    boost::dynamic_bitset<> barrier_mask{};
    std::vector<bp_t> moves_buff1{};
    std::vector<bp_t> moves_buff2{};
    std::vector<std::size_t> idx_buff1{};
    std::vector<std::size_t> idx_buff2{};
    std::vector<std::size_t> epoch_buff{};
    modle::PRNG_t rand_eng{};
    uint64_t seed{};

    inline State& operator=(const Task& task);
    inline void resize(std::size_t size = std::numeric_limits<std::size_t>::max());
    inline void reset();
  };

  inline void run();

 private:
  Genome _chromosomes{};

  [[nodiscard]] [[maybe_unused]] inline boost::asio::thread_pool instantiate_thread_pool() const;
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
                        absl::Span<std::size_t> fwd_lef_ranks, MaskT& mask, modle::PRNG_t& rand_eng,
                        bool first_epoch = false) noexcept(utils::ndebug_defined());

  //! Update CTCF states for the current iteration based on the states from the previous iteration.

  //! \param extr_barriers
  //! \param mask bitset used to store CTCF states. States will be updated inplace. Bits set to 0
  //!        represent CTCFs in NOT_OCCUPIED state, while bits set to 1 represent CTCFs in OCCUPIED
  //!        state
  //! \param rand_eng
  inline void generate_ctcf_states(absl::Span<const ExtrusionBarrier> extr_barriers,
                                   boost::dynamic_bitset<>& mask,
                                   modle::PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

  // clang-format off
  //! Generate moves in reverse direction

  //! This function ensures that generated moves will never cause an extrusion unit to move past
  //! chrom. boundaries.
  //! The callee must ensure that only extrusion units moving in rev. directions are given as the
  //! \p unit parameter
  //! \return move
  // clang-format on
  [[nodiscard]] inline bp_t generate_rev_move(const Chromosome* chrom, const ExtrusionUnit& unit,
                                              modle::PRNG_t& rand_eng);

  //! Generate moves in forward direction. See doc for Simulation::generate_rev_move for more
  //! details
  [[nodiscard]] inline bp_t generate_fwd_move(const Chromosome* chrom, const ExtrusionUnit& unit,
                                              modle::PRNG_t& rand_eng);
  // clang-format off
  //! Generate moves for all active LEFs

  //! Moves are generated by calling Simulation::generate_rev_move and Simulation::generate_fwd_move.
  //! When adjust_moves_ is true, adjust moves to make consecutive LEFs behave in a more realistic way.
  //! See Simulation::adjust_moves for more details
  // clang-format on
  inline void generate_moves(const Chromosome* chrom, absl::Span<const Lef> lefs,
                             absl::Span<const std::size_t> rev_lef_ranks,
                             absl::Span<const std::size_t> fwd_lef_ranks,
                             absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
                             modle::PRNG_t& rand_eng,
                             bool adjust_moves_ = true) noexcept(utils::ndebug_defined());

  // clang-format off
  //! Adjust moves of consecutive extr. units moving in the same direction to improve simulation realism

  //! Make sure that consecutive extr. units that are moving in the same direction do not bypass
  //! each other.
  //!
  //! Consider the following example:
  //!  - LEF1.fwd_unit.pos() = 0
  //!  - LEF2.fwd_unit.pos() = 10
  //!  - fwd_extrusion_speed = 1000
  //!  - fwd_extrusion_speed_std > 0
  //!
  //! For the current iteration, the fwd_unit of LEF1 is set to move 1050 bp, while that of LEF2
  //! is set to move only 950 bp.
  //! In this scenario, the fwd unit of LEF1 would bypass the fwd unit of LEF2.
  //! In a real system, what would likely happen, is that the fwd_unit of LEF1 would push
  //! against the fwd_unit of LEF2, temporarily increasing the fwd extr. speed of LEF2.
  // clang-format on
  inline void adjust_moves(const Chromosome* chrom, absl::Span<const Lef> lefs,
                           absl::Span<const std::size_t> rev_lef_ranks,
                           absl::Span<const std::size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
                           absl::Span<bp_t> fwd_moves) noexcept(utils::ndebug_defined());

  /// Sort extr. units by index based on their position whilst properly dealing with ties
  inline static void rank_lefs(absl::Span<const Lef> lefs,
                               absl::Span<std::size_t> rev_lef_rank_buff,
                               absl::Span<std::size_t> fwd_lef_rank_buff,
                               bool ranks_are_partially_sorted = true,
                               bool init_buffers = false) noexcept(utils::ndebug_defined());

  /// Extrude LEFs by applying the respective moves

  //! This function assumes that the move vectors that are passed as argument have already been
  //! processed by:
  //! - Simulation::adjust_moves
  //! - Simulation::detect_lef_bar_collisions
  //! - Simulation::detect_lef_lef_collisions
  inline void extrude(const Chromosome* chrom, absl::Span<Lef> lefs,
                      absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
                      std::size_t nrev_units_at_5prime = 0,
                      std::size_t nfwd_units_at_3prime = 0) noexcept(utils::ndebug_defined());

  //! This is just a wrapper function used to make sure that process_* functions are always called
  //! in the right order

  //! \return TODO
  template <typename I, typename MaskT>
  inline std::pair<std::size_t, std::size_t> process_collisions(
      const Chromosome* chrom, absl::Span<const Lef> lefs,
      absl::Span<const ExtrusionBarrier> barriers, const MaskT& barrier_mask,
      absl::Span<const std::size_t> rev_lef_ranks, absl::Span<const std::size_t> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions, PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

  template <typename I>
  static inline std::pair<std::size_t, std::size_t> detect_collisions_at_chrom_boundaries(
      const Chromosome* chrom, absl::Span<const Lef> lefs,
      absl::Span<const std::size_t> rev_lef_ranks, absl::Span<const std::size_t> fwd_lef_ranks,
      absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
      absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions);

  /// Detect and process collisions between LEFs and extrusion barriers

  //! After calling this function, the entries in \p rev_collisions and \p fwd_collisions
  //! corresponding to reverse or forward extrusion units that will collide with an extrusion
  //! barrier in the current iteration, will be set to the index corresponding to the barrier that
  //! is causing the collision. Furthermore, the moves of the extrusion units involved in collision
  //! events are updated so that after calling Simulation::extrude, these extrusion units will be
  //! located 1bp up/downstream (depending on extrusion direction) of the extrusion barrier that is
  //! blocking them.
  template <typename I>
  inline void detect_lef_bar_collisions(
      absl::Span<const Lef> lefs, absl::Span<const std::size_t> rev_lef_ranks,
      absl::Span<const std::size_t> fwd_lef_ranks, absl::Span<const bp_t> rev_moves,
      absl::Span<const bp_t> fwd_moves, absl::Span<const ExtrusionBarrier> extr_barriers,
      const boost::dynamic_bitset<>& barrier_mask, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions, modle::PRNG_t& rand_eng, std::size_t nrev_units_at_5prime = 0,
      std::size_t nfwd_units_at_3prime = 0) noexcept(utils::ndebug_defined());

  /// Detect and process collisions between LEFs

  //! This is just a wrapper that checks whether the first rev unit and the last fwd units are at
  //! chromosomal boundaries (and sets the respective collision mask to REACHED_CHROM_BOUNDARIES
  //! when appropriate), then it calls Simulation::process_major_lef_lef_collisions and
  //! Simulation::process_minor_lef_lef_collisions. See documentation of those two functions for
  //! more details on the actions that are performed.
  template <typename I>
  inline void detect_lef_lef_collisions(
      const Chromosome* chrom, absl::Span<const Lef> lefs,
      absl::Span<const ExtrusionBarrier> barriers, absl::Span<const std::size_t> rev_lef_ranks,
      absl::Span<const std::size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
      absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions,
      PRNG_t& rand_eng, std::size_t nrev_units_at_5prime = 0,
      std::size_t nfwd_units_at_3prime = 0) noexcept(utils::ndebug_defined());

  template <typename I>
  inline void detect_primary_lef_lef_collisions(
      const Chromosome* chrom, absl::Span<const Lef> lefs,
      absl::Span<const ExtrusionBarrier> barriers, absl::Span<const std::size_t> rev_lef_ranks,
      absl::Span<const std::size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
      absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions,
      PRNG_t& rand_eng, std::size_t nrev_units_at_5prime = 0,
      std::size_t nfwd_units_at_3prime = 0) noexcept(utils::ndebug_defined());

  template <typename I>
  inline void detect_secondary_lef_lef_collisions(
      const Chromosome* chrom, absl::Span<const Lef> lefs, std::size_t nbarriers,
      absl::Span<const std::size_t> rev_lef_ranks, absl::Span<const std::size_t> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions, PRNG_t& rand_eng, std::size_t nrev_units_at_5prime = 0,
      std::size_t nfwd_units_at_3prime = 0) noexcept(utils::ndebug_defined());

  // TODO: Figure out how to write a template that takes Span<const int> as template argument
  template <typename I>
  static inline void adjust_moves_for_lef_bar_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined());

  template <typename I>
  static inline void adjust_moves_for_primary_lef_lef_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<const std::size_t> rev_ranks, absl::Span<const std::size_t> fwd_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined());

  template <typename I>
  static inline void adjust_moves_for_secondary_lef_lef_collisions(
      absl::Span<const Lef> lefs, std::size_t nbarriers, absl::Span<const std::size_t> rev_ranks,
      absl::Span<const std::size_t> fwd_ranks, absl::Span<bp_t> rev_moves,
      absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined());

  template <typename I>
  static inline void adjust_moves_based_on_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<const std::size_t> rev_ranks, absl::Span<const std::size_t> fwd_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions) noexcept(utils::ndebug_defined());

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

  [[nodiscard]] static inline std::pair<bp_t /*rev*/, bp_t /*fwd*/> compute_lef_lef_collision_pos(
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
                                    absl::Span<std::size_t> fwd_lef_rank_buff,
                                    bool ranks_are_partially_sorted = true,
                                    bool init_buffers = false) {
    Simulation::rank_lefs(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, ranks_are_partially_sorted,
                          init_buffers);
  }

  template <typename I>
  inline void test_process_units_at_chrom_boundaries(
      const Chromosome* chrom, absl::Span<const Lef> lefs,
      absl::Span<const std::size_t> rev_lef_ranks, absl::Span<const std::size_t> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions) {
    this->detect_collisions_at_chrom_boundaries(chrom, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                rev_moves, fwd_moves, rev_collisions,
                                                fwd_collisions);
  }

  template <typename I>
  inline void test_process_lef_lef_collisions(
      const Chromosome* const chrom, absl::Span<const Lef> lefs,
      absl::Span<const ExtrusionBarrier> barriers, absl::Span<const std::size_t> rev_lef_rank_buff,
      absl::Span<const std::size_t> fwd_lef_rank_buff, absl::Span<bp_t> rev_move_buff,
      absl::Span<bp_t> fwd_move_buff, absl::Span<I> rev_collision_buff,
      absl::Span<I> fwd_collision_buff, modle::PRNG_t& rand_eng) {
    this->detect_lef_lef_collisions(chrom, lefs, barriers, rev_lef_rank_buff, fwd_lef_rank_buff,
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
    this->detect_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff,
                                    fwd_move_buff, extr_barriers, barrier_mask, rev_collisions,
                                    fwd_collisions, rand_eng);
  }

  template <typename I>
  static inline void test_adjust_moves_for_lef_bar_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions) {
    Simulation::adjust_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                    rev_collisions, fwd_collisions);
  }

  template <typename I>
  static inline void test_adjust_moves_for_primary_lef_lef_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<const std::size_t> rev_ranks, absl::Span<const std::size_t> fwd_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions) {
    Simulation::adjust_moves_for_primary_lef_lef_collisions(
        lefs, barriers, rev_ranks, fwd_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);
  }

  template <typename I>
  static inline void test_adjust_moves_for_secondary_lef_lef_collisions(
      absl::Span<const Lef> lefs, std::size_t nbarriers, absl::Span<const std::size_t> rev_ranks,
      absl::Span<const std::size_t> fwd_ranks, absl::Span<bp_t> rev_moves,
      absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions, absl::Span<I> fwd_collisions) {
    Simulation::adjust_moves_for_secondary_lef_lef_collisions(lefs, nbarriers, rev_ranks, fwd_ranks,
                                                              rev_moves, fwd_moves, rev_collisions,
                                                              fwd_collisions);
  }

  template <typename I>
  static inline void test_adjust_moves_based_on_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<const std::size_t> rev_ranks, absl::Span<const std::size_t> fwd_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<I> rev_collisions,
      absl::Span<I> fwd_collisions) {
    Simulation::adjust_moves_based_on_collisions(lefs, barriers, rev_ranks, fwd_ranks, rev_moves,
                                                 fwd_moves, rev_collisions, fwd_collisions);
  }
#endif
};
}  // namespace modle

#include "../../simulation_impl.hpp"  // IWYU pragma: keep
