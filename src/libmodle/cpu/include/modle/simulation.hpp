#pragma once

#include <absl/container/btree_set.h>            // IWYU pragma: keep for btree_set
#include <absl/container/flat_hash_map.h>        // for flat_hash_map
#include <absl/types/span.h>                     // IWYU pragma: keep for Span
#include <moodycamel/blockingconcurrentqueue.h>  // for BlockingConcurrntQueue

#include <atomic>                                   // for atomic
#include <boost/asio/thread_pool.hpp>               // for thread_pool
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cstddef>                                  // for size_t
#include <cstdint>                                  // for uint64_t
#include <deque>                                    // for deque
#include <filesystem>                               // for path
#include <limits>                                   // for numeric_limits
#include <memory>                                   // for unique_ptr
#include <mutex>                                    // for mutex
#include <random>                                   // for normal_distribution, uniform_int_dist...
#include <string>                                   // for string
#include <utility>                                  // for pair
#include <vector>                                   // for vector

#include "modle/common.hpp"             // for bp_t, PRNG, seeder
#include "modle/config.hpp"             // for Config
#include "modle/dna.hpp"                // for Chromosome
#include "modle/extrusion_factors.hpp"  // for Lef, ExtrusionUnit (ptr only)
#include "modle/setup.hpp"              // for Genome
#include "modle/suppress_compiler_warnings.hpp"
#include "modle/utils.hpp"  // for ndebug_defined

namespace modle {

class ExtrusionBarrier;
class ExtrusionUnit;
struct Lef;

class Simulation : Config {
 public:
  explicit Simulation(const Config& c, bool import_chroms = true);
  std::string to_string() = delete;
  void print() = delete;

  [[nodiscard]] size_t size() const;
  [[nodiscard]] size_t simulated_size() const;

  using lef_move_generator_t = std::normal_distribution<double>;
  using chrom_pos_generator_t = std::uniform_int_distribution<bp_t>;

  static constexpr auto NO_COLLISION = std::numeric_limits<collision_t>::max();
  static constexpr auto LEF_LEF_COLLISION = std::numeric_limits<collision_t>::max() - 1;
  static constexpr auto REACHED_CHROM_BOUNDARY = std::numeric_limits<collision_t>::max() - 2;

  struct Task {
    Task() = default;
    size_t id{};
    Chromosome* chrom{};
    size_t cell_id{};
    size_t n_target_epochs{};
    size_t n_target_contacts{};
    size_t nlefs{};
    absl::Span<const ExtrusionBarrier> barriers{};
    // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  } __attribute__((aligned(64)));

  struct State : Task {
    State() = default;
    std::vector<Lef> lef_buff{};
    std::vector<double> lef_unloader_affinity{};
    std::vector<size_t> rank_buff1{};
    std::vector<size_t> rank_buff2{};
    boost::dynamic_bitset<> barrier_mask{};
    std::vector<bp_t> moves_buff1{};
    std::vector<bp_t> moves_buff2{};
    std::vector<size_t> idx_buff1{};
    std::vector<size_t> idx_buff2{};
    std::vector<size_t> epoch_buff{};
    modle::PRNG_t rand_eng{};
    uint64_t seed{};

    State& operator=(const Task& task);
    void resize(size_t size = std::numeric_limits<size_t>::max());
    void reset();
    [[nodiscard]] std::string to_string() const noexcept;
    // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  } __attribute__((aligned(128)));

  void run();

 private:
  const Config* _config;
  Genome _genome{};

  [[nodiscard]] [[maybe_unused]] boost::asio::thread_pool instantiate_thread_pool() const;
  template <typename I>
  [[nodiscard]] static boost::asio::thread_pool instantiate_thread_pool(
      I nthreads, bool clamp_nthreads = true) {
    static_assert(std::is_integral_v<I>, "nthreads should have an integral type.");
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_USELESS_CAST
    if (clamp_nthreads) {
      return boost::asio::thread_pool(
          std::min(std::thread::hardware_concurrency(), static_cast<unsigned int>(nthreads)));
    }
    assert(nthreads > 0);
    return boost::asio::thread_pool(static_cast<unsigned int>(nthreads));
    DISABLE_WARNING_POP
  }

  /// Simulate loop extrusion using the parameters and buffers passed through \p state
  void simulate_extrusion_kernel(State& state) const;

  /// Monitor the progress queue and write contacts to a .cool file when all the cells belonging to
  /// a given chromosome have been simulated.

  //! IMPORTANT: this function is meant to be run in a dedicated thread.
  void write_contacts_to_disk(std::deque<std::pair<Chromosome*, size_t>>& progress_queue,
                              std::mutex& progress_queue_mutex,
                              std::atomic<bool>& end_of_simulation) const;

  /// Worker function used to run an instance of the simulation

  //! Worker function used to consume Simulation::Tasks from a task queue, setup a
  //! Simulation::State, then run an instance of the simulation of a specific chromosome (i.e.
  //! simulate loop extrusion on a single chromosome in a cell).
  //! This function is also responsible for allocating and clearing the buffers used throughout the
  //! simulation.
  //! IMPORTANT: this function is meant to be run in a dedicated thread.
  void worker(
      moodycamel::BlockingConcurrentQueue<Simulation::Task>& task_queue,
      std::deque<std::pair<Chromosome*, size_t>>& progress_queue, std::mutex& progress_queue_mutex,
      std::mutex& barrier_mutex,
      absl::flat_hash_map<Chromosome*, std::unique_ptr<std::vector<ExtrusionBarrier>>>& barriers,
      std::atomic<bool>& end_of_simulation, size_t task_batch_size = 32) const;  // NOLINT

  /// Bind inactive LEFs, then sort them by their genomic coordinates.

  //! Bind the LEFs passed through \p lefs whose corresponding entry in \p mask is set to true, then
  //! rank rev and fwd units based on their genomic coordinates.
  //! In order to properly handle ties, the sorting criterion is actually a bit more complex. See
  //! documentation and comments for Simulation::rank_lefs for more details).
  template <typename MaskT>
  static void bind_lefs(const Chromosome* chrom, absl::Span<Lef> lefs,
                        absl::Span<size_t> rev_lef_ranks, absl::Span<size_t> fwd_lef_ranks,
                        const MaskT& mask, modle::PRNG_t& rand_eng,
                        bool first_epoch = false) noexcept(utils::ndebug_defined());

  // clang-format off
  //! Generate moves in reverse direction

  //! This function ensures that generated moves will never cause an extrusion unit to move past
  //! chrom. boundaries.
  //! The callee must ensure that only extrusion units moving in rev. directions are given as the
  //! \p unit parameter
  //! \return move
  // clang-format on
  [[nodiscard]] bp_t generate_rev_move(const Chromosome* chrom, const ExtrusionUnit& unit,
                                       modle::PRNG_t& rand_eng) const;

  //! Generate moves in forward direction. See doc for Simulation::generate_rev_move for more
  //! details
  [[nodiscard]] bp_t generate_fwd_move(const Chromosome* chrom, const ExtrusionUnit& unit,
                                       modle::PRNG_t& rand_eng) const;
  // clang-format off
  //! Generate moves for all active LEFs

  //! Moves are generated by calling Simulation::generate_rev_move and Simulation::generate_fwd_move.
  //! When adjust_moves_ is true, adjust moves to make consecutive LEFs behave in a more realistic way.
  //! See Simulation::adjust_moves_of_consecutive_extr_units for more details
  // clang-format on
  void generate_moves(const Chromosome* chrom, absl::Span<const Lef> lefs,
                      absl::Span<const size_t> rev_lef_ranks,
                      absl::Span<const size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
                      absl::Span<bp_t> fwd_moves, modle::PRNG_t& rand_eng,
                      bool adjust_moves_ = true) const noexcept(utils::ndebug_defined());

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
  void adjust_moves_of_consecutive_extr_units(const Chromosome* chrom, absl::Span<const Lef> lefs,
                                              absl::Span<const size_t> rev_lef_ranks,
                                              absl::Span<const size_t> fwd_lef_ranks,
                                              absl::Span<bp_t> rev_moves,
                                              absl::Span<bp_t> fwd_moves) const
      noexcept(utils::ndebug_defined());

  /// Sort extr. units by index based on their position whilst properly dealing with ties. See
  /// comments in \p simulation_impl.hpp file for more details.
  static void rank_lefs(absl::Span<const Lef> lefs, absl::Span<size_t> rev_lef_rank_buff,
                        absl::Span<size_t> fwd_lef_rank_buff,
                        bool ranks_are_partially_sorted = true,
                        bool init_buffers = false) noexcept(utils::ndebug_defined());

  /// Extrude LEFs by applying the respective moves

  //! This function assumes that the move vectors that are passed as argument have already been
  //! processed by:
  //! - Simulation::adjust_moves_of_consecutive_extr_units
  //! - Simulation::process_collisions
  static void extrude(const Chromosome* chrom, absl::Span<Lef> lefs,
                      absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
                      size_t nrev_units_at_5prime = 0,
                      size_t nfwd_units_at_3prime = 0) noexcept(utils::ndebug_defined());

  //! This is just a wrapper function used to make sure that process_* functions are always called
  //! in the right order

  //! \return number of rev units at the 5'-end, number of fwd units at the 3'-end
  std::pair<size_t, size_t> process_collisions(
      const Chromosome* chrom, absl::Span<const Lef> lefs,
      absl::Span<const ExtrusionBarrier> barriers, const boost::dynamic_bitset<>& barrier_mask,
      absl::Span<const size_t> rev_lef_ranks, absl::Span<const size_t> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<collision_t> rev_collisions, absl::Span<collision_t> fwd_collisions,
      PRNG_t& rand_eng) const noexcept(utils::ndebug_defined());

  /// Detect and stall LEFs with one or more extrusion units located at chromosomal boundaries.

  //! \return a pair of numbers consisting in the number of rev units located at the 5'-end and
  //! number of fwd units located at the 3'-end.
  static std::pair<size_t, size_t> detect_collisions_at_chrom_boundaries(
      const Chromosome* chrom, absl::Span<const Lef> lefs, absl::Span<const size_t> rev_lef_ranks,
      absl::Span<const size_t> fwd_lef_ranks, absl::Span<const bp_t> rev_moves,
      absl::Span<const bp_t> fwd_moves, absl::Span<collision_t> rev_collisions,
      absl::Span<collision_t> fwd_collisions);

  /// Detect collisions between LEFs and extrusion barriers.

  //! After calling this function, the entries in \p rev_collisions and \p fwd_collisions
  //! corresponding to reverse or forward extrusion units that will collide with an extrusion
  //! barrier in the current iteration, will be set to the index corresponding to the barrier that
  //! is causing the collision.
  void detect_lef_bar_collisions(absl::Span<const Lef> lefs, absl::Span<const size_t> rev_lef_ranks,
                                 absl::Span<const size_t> fwd_lef_ranks,
                                 absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
                                 absl::Span<const ExtrusionBarrier> extr_barriers,
                                 const boost::dynamic_bitset<>& barrier_mask,
                                 absl::Span<collision_t> rev_collisions,
                                 absl::Span<collision_t> fwd_collisions, modle::PRNG_t& rand_eng,
                                 size_t nrev_units_at_5prime = 0,
                                 size_t nfwd_units_at_3prime = 0) const
      noexcept(utils::ndebug_defined());

  /// Detect collisions between LEFs moving in opposite directions.

  //! Primary LEF-LEF collisions can occur when two extrusion units moving in opposite direction
  //! collide with each other.
  //! After calling this function, the entries in \p rev_collisions and \p fwd_collisions
  //! corresponding to units stalled due to primary LEF-LEF collisions will be set to the index
  //! pointing to the extrusion unit that is causing the collisions.
  //! The index i is encoded as nbarriers + i.
  void detect_primary_lef_lef_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<const size_t> rev_lef_ranks, absl::Span<const size_t> fwd_lef_ranks,
      absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
      absl::Span<collision_t> rev_collisions, absl::Span<collision_t> fwd_collisions,
      PRNG_t& rand_eng, size_t nrev_units_at_5prime = 0, size_t nfwd_units_at_3prime = 0) const
      noexcept(utils::ndebug_defined());

  // clang-format off
  /// Process collisions between consecutive extrusion units that are moving in the same direction.

  //! Secondary LEF-LEF collisions can occur when consecutive extrusion units moving in the same
  //! direction collide with each other.
  //! This can happen in the following scenario:
  //! Extrusion units U1 and U2 are located one after the other in 5'-3' direction with U1 preceding U2.
  //! Both units are extruding in fwd direction. U2 is blocked by an extrusion barrier.
  //! If U1 + move1 is greater than U2 + move2 (after adjusting move2 to comply with the LEF-BAR collision constraint)
  //! then U1 will incur in a secondary LEF-LEF collision, and its move will be adjusted so that after calling
  //! Simulation::extrude, U1 will be located 1bp upstream of U2.
  //! When a secondary LEF-LEF collision is detected, the appropriate entry in \p rev_collisions or
  //! \p fwd_collisions will be set to the index corresponding to the extrusion unit that is causing the collision.
  //! The index i is encoded as nbarriers + nlefs + i.
  // clang-format on
  void process_secondary_lef_lef_collisions(
      const Chromosome* chrom, absl::Span<const Lef> lefs, size_t nbarriers,
      absl::Span<const size_t> rev_lef_ranks, absl::Span<const size_t> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<collision_t> rev_collisions, absl::Span<collision_t> fwd_collisions,
      PRNG_t& rand_eng, size_t nrev_units_at_5prime = 0, size_t nfwd_units_at_3prime = 0) const
      noexcept(utils::ndebug_defined());

  /// Correct moves to comply with the constraints imposed by LEF-BAR collisions.
  static void correct_moves_for_lef_bar_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<const collision_t> rev_collisions,
      absl::Span<const collision_t> fwd_collisions) noexcept(utils::ndebug_defined());

  /// Correct moves to comply with the constraints imposed by primary LEF-LEF collisions.
  static void correct_moves_for_primary_lef_lef_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<const size_t> rev_ranks, absl::Span<const size_t> fwd_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<const collision_t> rev_collisions,
      absl::Span<const collision_t> fwd_collisions) noexcept(utils::ndebug_defined());
/*
  /// Correct moves to comply with the constraints imposed by secondary LEF-LEF collisions.
  static void correct_moves_for_secondary_lef_lef_collisions(
      absl::Span<const Lef> lefs, size_t nbarriers, absl::Span<const size_t> rev_ranks,
      absl::Span<const size_t> fwd_ranks, absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<const collision_t> rev_collisions, absl::Span<const collision_t> fwd_collisions,
      size_t nrev_units_at_5prime = 0,
      size_t nfwd_units_at_3prime = 0) noexcept(utils::ndebug_defined());
*/
  /// Register contacts for chromosome \p chrom using the position the extrusion units of the LEFs
  /// in \p lefs whose index is present in \p selected_lef_idx.
  size_t register_contacts(Chromosome* chrom, absl::Span<const Lef> lefs,
                           absl::Span<const size_t> selected_lef_idx) const
      noexcept(utils::ndebug_defined());

  static void select_lefs_to_bind(absl::Span<const Lef> lefs,
                                  boost::dynamic_bitset<>& mask) noexcept(utils::ndebug_defined());
  static void select_lefs_to_bind(absl::Span<const Lef> lefs,
                                  absl::Span<size_t> mask) noexcept(utils::ndebug_defined());

  void generate_lef_unloader_affinities(absl::Span<const Lef> lefs,
                                        absl::Span<const ExtrusionBarrier> barriers,
                                        absl::Span<const collision_t> rev_collisions,
                                        absl::Span<const collision_t> fwd_collisions,
                                        absl::Span<double> lef_unloader_affinity) const
      noexcept(utils::ndebug_defined());

  static void select_lefs_to_release(absl::Span<size_t> lef_idx,
                                     absl::Span<const double> lef_unloader_affinity,
                                     modle::PRNG_t& rand_eng) noexcept(utils::ndebug_defined());

  static void release_lefs(absl::Span<Lef> lefs, absl::Span<const size_t> lef_idx) noexcept;

  [[nodiscard]] static std::pair<bp_t /*rev*/, bp_t /*fwd*/> compute_lef_lef_collision_pos(
      const ExtrusionUnit& rev_unit, const ExtrusionUnit& fwd_unit, bp_t rev_move, bp_t fwd_move);

#ifdef ENABLE_TESTING
 public:
  template <typename MaskT>
  inline void test_bind_lefs(const Chromosome* chrom, absl::Span<Lef> lefs,
                             absl::Span<size_t> rev_lef_ranks, absl::Span<size_t> fwd_lef_ranks,
                             const MaskT& mask, modle::PRNG_t& rand_eng) {
    modle::Simulation::bind_lefs(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, mask, rand_eng);
  }

  inline void test_adjust_moves(const Chromosome* const chrom, absl::Span<const Lef> lefs,
                                absl::Span<const size_t> rev_lef_ranks,
                                absl::Span<const size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
                                absl::Span<bp_t> fwd_moves) {
    this->adjust_moves_of_consecutive_extr_units(chrom, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                 rev_moves, fwd_moves);
  }

  inline void test_generate_moves(const Chromosome* const chrom, absl::Span<const Lef> lefs,
                                  absl::Span<const size_t> rev_lef_ranks,
                                  absl::Span<const size_t> fwd_lef_ranks,
                                  absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
                                  modle::PRNG_t& rand_eng, bool adjust_moves_ = false) {
    this->generate_moves(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, rand_eng,
                         adjust_moves_);
  }

  inline static void test_rank_lefs(absl::Span<const Lef> lefs,
                                    absl::Span<size_t> rev_lef_rank_buff,
                                    absl::Span<size_t> fwd_lef_rank_buff,
                                    bool ranks_are_partially_sorted = true,
                                    bool init_buffers = false) {
    Simulation::rank_lefs(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, ranks_are_partially_sorted,
                          init_buffers);
  }

  inline void test_process_units_at_chrom_boundaries(
      const Chromosome* chrom, absl::Span<const Lef> lefs, absl::Span<const size_t> rev_lef_ranks,
      absl::Span<const size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
      absl::Span<bp_t> fwd_moves, absl::Span<collision_t> rev_collisions,
      absl::Span<collision_t> fwd_collisions) {
    Simulation::detect_collisions_at_chrom_boundaries(chrom, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                      rev_moves, fwd_moves, rev_collisions,
                                                      fwd_collisions);
  }

  inline void test_process_lef_bar_collisions(
      absl::Span<const Lef> lefs, absl::Span<const size_t> rev_lef_rank_buff,
      absl::Span<const size_t> fwd_lef_rank_buff, absl::Span<bp_t> rev_move_buff,
      absl::Span<bp_t> fwd_move_buff, absl::Span<const ExtrusionBarrier> extr_barriers,
      const boost::dynamic_bitset<>& barrier_mask, absl::Span<collision_t> rev_collisions,
      absl::Span<collision_t> fwd_collisions, modle::PRNG_t& rand_eng) {
    Simulation::detect_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff,
                                          fwd_move_buff, extr_barriers, barrier_mask,
                                          rev_collisions, fwd_collisions, rand_eng);
  }

  static inline void test_adjust_moves_for_lef_bar_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<collision_t> rev_collisions, absl::Span<collision_t> fwd_collisions) {
    Simulation::correct_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                     rev_collisions, fwd_collisions);
  }

  static inline void test_adjust_moves_for_primary_lef_lef_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<const size_t> rev_ranks, absl::Span<const size_t> fwd_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<collision_t> rev_collisions, absl::Span<collision_t> fwd_collisions) {
    Simulation::correct_moves_for_primary_lef_lef_collisions(
        lefs, barriers, rev_ranks, fwd_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);
  }

  static inline void test_adjust_moves_for_secondary_lef_lef_collisions(
      absl::Span<const Lef> lefs, size_t nbarriers, absl::Span<const size_t> rev_ranks,
      absl::Span<const size_t> fwd_ranks, absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<collision_t> rev_collisions, absl::Span<collision_t> fwd_collisions) {
    Simulation::correct_moves_for_secondary_lef_lef_collisions(lefs, nbarriers, rev_ranks,
                                                               fwd_ranks, rev_moves, fwd_moves,
                                                               rev_collisions, fwd_collisions);
  }
#endif
};
}  // namespace modle

#include "../../simulation_impl.hpp"  // IWYU pragma: keep
