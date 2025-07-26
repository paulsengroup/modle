// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

#include <atomic>
#include <deque>
#include <exception>
#include <filesystem>
#include <limits>
#include <memory>
#include <mutex>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "modle/bed/bed.hpp"
#include "modle/collision_encoding.hpp"
#include "modle/common/common.hpp"
#include "modle/common/random.hpp"
#include "modle/common/simulation_config.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"
#include "modle/common/utils.hpp"
#include "modle/contact_matrix_dense.hpp"
#include "modle/context_manager.hpp"
#include "modle/extrusion_barriers.hpp"
#include "modle/extrusion_factors.hpp"
#include "modle/genome.hpp"

namespace modle {

namespace compressed_io {
class Writer;
}

template <typename I, typename T>
class IITree;

class Simulation {
 public:
  explicit Simulation(Config config_, bool import_intervals = true);
  void print() = delete;

  [[nodiscard]] std::size_t size() const;
  [[nodiscard]] std::size_t simulated_size() const;

  using lef_move_generator_t = random::normal_distribution<double>;
  using chrom_pos_generator_t = random::uniform_int_distribution<bp_t>;

  static constexpr auto Mbp = 1.0e6;

  using CollisionT = Collision<std::uint_fast32_t>;
  struct Task {  // NOLINT(altera-struct-pack-align)
    enum class Status : std::uint_fast8_t { PENDING, RUNNING, COMPLETED, FAILED };
    std::size_t id{};
    GenomicInterval* interval{};
    std::size_t cell_id{};
    std::size_t num_target_epochs{};
    std::size_t num_target_contacts{};
    std::size_t num_lefs{};
    random::PRNG_t rand_eng{};
    Status status{Status::PENDING};
  };

  struct State : Task {  // NOLINT(altera-struct-pack-align)
    State() = default;
    std::size_t epoch{};             // NOLINT
    bool burnin_completed{false};    // NOLINT
    std::size_t num_active_lefs{0};  // NOLINT

    std::size_t num_burnin_epochs{0};  // NOLINT
    std::size_t num_contacts{0};       // NOLINT

    std::uint64_t seed{};                                                // NOLINT
    std::unique_ptr<compressed_io::Writer> model_state_logger{nullptr};  // NOLINT

    ExtrusionBarriers barriers{};

   protected:
    std::vector<Lef> lef_buff{};                 // NOLINT
    std::vector<std::size_t> rank_buff1{};       // NOLINT
    std::vector<std::size_t> rank_buff2{};       // NOLINT
    ExtrusionBarriers _barriers{};               // NOLINT
    std::vector<bp_t> moves_buff1{};             // NOLINT
    std::vector<bp_t> moves_buff2{};             // NOLINT
    std::vector<std::size_t> idx_buff{};         // NOLINT
    std::vector<CollisionT> collision_buff1{};   // NOLINT
    std::vector<CollisionT> collision_buff2{};   // NOLINT
    std::deque<double> cfx_of_variation_buff{};  // NOLINT
    std::deque<double> avg_loop_size_buff{};     // NOLINT

    static constexpr std::size_t npos =
        std::numeric_limits<std::span<std::size_t>::size_type>::max();

   public:
    [[nodiscard]] std::span<Lef> get_lefs(std::size_t size = npos) noexcept;
    [[nodiscard]] std::span<std::size_t> get_rev_ranks(std::size_t size = npos) noexcept;
    [[nodiscard]] std::span<std::size_t> get_fwd_ranks(std::size_t size = npos) noexcept;
    [[nodiscard]] std::span<bp_t> get_rev_moves(std::size_t size = npos) noexcept;
    [[nodiscard]] std::span<bp_t> get_fwd_moves(std::size_t size = npos) noexcept;
    [[nodiscard]] std::span<std::size_t> get_idx_buff(std::size_t size = npos) noexcept;
    [[nodiscard]] std::span<CollisionT> get_rev_collisions(std::size_t size = npos) noexcept;
    [[nodiscard]] std::span<CollisionT> get_fwd_collisions(std::size_t size = npos) noexcept;
    [[nodiscard]] std::deque<double>& get_cfx_of_variation() noexcept;
    [[nodiscard]] std::deque<double>& get_avg_loop_sizes() noexcept;

    [[nodiscard]] std::span<const Lef> get_lefs(std::size_t size = npos) const noexcept;
    [[nodiscard]] std::span<const std::size_t> get_rev_ranks(
        std::size_t size = npos) const noexcept;
    [[nodiscard]] std::span<const std::size_t> get_fwd_ranks(
        std::size_t size = npos) const noexcept;
    [[nodiscard]] std::span<const bp_t> get_rev_moves(std::size_t size = npos) const noexcept;
    [[nodiscard]] std::span<const bp_t> get_fwd_moves(std::size_t size = npos) const noexcept;
    [[nodiscard]] std::span<const std::size_t> get_idx_buff(std::size_t size = npos) const noexcept;
    [[nodiscard]] std::span<const CollisionT> get_rev_collisions(
        std::size_t size = npos) const noexcept;
    [[nodiscard]] std::span<const CollisionT> get_fwd_collisions(
        std::size_t size = npos) const noexcept;
    [[nodiscard]] const std::deque<double>& get_cfx_of_variation() const noexcept;
    [[nodiscard]] const std::deque<double>& get_avg_loop_sizes() const noexcept;

    std::shared_ptr<std::mutex> contacts_mtx{nullptr};                  // NOLINT
    std::shared_ptr<ContactMatrixDense<contacts_t>> contacts{nullptr};  // NOLINT

    State& operator=(const Task& task);

    void resize_buffers(std::size_t size = (std::numeric_limits<std::size_t>::max)());
    void reset_buffers();
  };

 private:
  Config _config{};
  Genome _genome{};
  ContextManager<Task> _ctx{0};

  static constexpr auto& model_internal_state_log_header = Config::model_internal_state_log_header;

 public:
  [[nodiscard]] constexpr const Config& config() const noexcept;
  [[nodiscard]] constexpr const Config& c() const noexcept;

  void spawn_worker_threads(std::size_t batch_size);
  void spawn_io_threads();
  void run_simulate();

 private:
  /// Simulate loop extrusion using the parameters and buffers passed through \p state
  void simulate_one_cell(std::uint64_t tid, State& s) const;

  //! IMPORTANT: this function is meant to be run in a dedicated thread.
  void simulate_io(std::chrono::milliseconds wait_time = std::chrono::milliseconds(50));

  /// Worker function used to run an instance of the simulation

  //! Worker function used to consume Simulation::Tasks from a task queue, setup a
  //! Simulation::State, then run an instance of the simulation of a specific genomic interval (i.e.
  //! simulate loop extrusion on a single chromosome in a cell).
  //! This function is also responsible for allocating and clearing the buffers used throughout the
  //! simulation.
  //! IMPORTANT: this function is meant to be run in a dedicated thread.
  void simulate_worker(std::uint64_t tid, std::size_t task_batch_size);

  /// Bind inactive LEFs, then sort them by their genomic coordinates.

  //! Bind the LEFs passed through \p lefs whose corresponding entry in \p mask is set to true, then
  //! rank rev and fwd units based on their genomic coordinates.
  //! In order to properly handle ties, the sorting criterion is actually a bit more complex. See
  //! documentation and comments for Simulation::rank_lefs for more details).
  template <typename MaskT>
  inline static void bind_lefs(const GenomicInterval& interval, std::span<Lef> lefs,
                               std::span<std::size_t> rev_lef_ranks,
                               std::span<std::size_t> fwd_lef_ranks, const MaskT& mask,
                               random::PRNG_t& rand_eng,
                               std::size_t current_epoch) noexcept(utils::ndebug_defined());

  template <typename MaskT>
  inline static void bind_lefs(bp_t start_pos, bp_t end_pos, std::span<Lef> lefs,
                               std::span<std::size_t> rev_lef_ranks,
                               std::span<std::size_t> fwd_lef_ranks, const MaskT& mask,
                               random::PRNG_t& rand_eng,
                               std::size_t current_epoch) noexcept(utils::ndebug_defined());

  static void select_and_bind_lefs(State& s) noexcept(utils::ndebug_defined());

  // clang-format off
  //! Generate moves for all active LEFs

  //! Moves are generated by calling Simulation::generate_rev_move and Simulation::generate_fwd_move.
  //! When adjust_moves_ is true, adjust moves to make consecutive LEFs behave in a more realistic way.
  //! See Simulation::adjust_moves_of_consecutive_extr_units for more details
  // clang-format on
  void generate_moves(const GenomicInterval& interval, std::span<const Lef> lefs,
                      std::span<const std::size_t> rev_lef_ranks,
                      std::span<const std::size_t> fwd_lef_ranks, std::span<bp_t> rev_moves,
                      std::span<bp_t> fwd_moves, bool burnin_completed, random::PRNG_t& rand_eng,
                      bool adjust_moves_ = true) const noexcept(utils::ndebug_defined());

  // clang-format off
  //! Adjust moves of consecutive extr. units moving in the same direction to improve simulation realism

  //! Make sure that consecutive extr. units that are moving in the same direction do not bypass
  //! each other. Also ensure that moves won't cause extrusion units to move past chromosomal boundaries.
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
  static void adjust_moves_of_consecutive_extr_units(
      const GenomicInterval& interval, std::span<const Lef> lefs,
      std::span<const std::size_t> rev_lef_ranks, std::span<const std::size_t> fwd_lef_ranks,
      std::span<bp_t> rev_moves, std::span<bp_t> fwd_moves) noexcept(utils::ndebug_defined());

  /// Clamp moves to prevent LEFs from falling off chromosomal boundaries
  static void clamp_moves(const GenomicInterval& interval, std::span<const Lef> lefs,
                          std::span<bp_t> rev_moves, std::span<bp_t> fwd_moves) noexcept;

  /// Sort extr. units by index based on their position whilst properly dealing with ties. See
  /// comments in \p simulation_impl.hpp file for more details.
  static void rank_lefs(std::span<const Lef> lefs, std::span<std::size_t> rev_lef_ranks,
                        std::span<std::size_t> fwd_lef_ranks,
                        bool ranks_are_partially_sorted = true,
                        bool init_buffers = false) noexcept(utils::ndebug_defined());

  /// Extrude LEFs by applying the respective moves

  //! This function assumes that the move vectors that are passed as argument have already been
  //! processed by:
  //! - Simulation::adjust_moves_of_consecutive_extr_units
  //! - Simulation::process_collisions
  static void extrude(const GenomicInterval& interval, std::span<Lef> lefs,
                      std::span<const bp_t> rev_moves,
                      std::span<const bp_t> fwd_moves) noexcept(utils::ndebug_defined());

  //! This is just a wrapper function used to make sure that process_* functions are always called
  //! in the right order

  //! \return number of rev units at the 5'-end, number of fwd units at the 3'-end
  std::pair<std::size_t, std::size_t> process_collisions(
      const GenomicInterval& interval, std::span<Lef> lefs, ExtrusionBarriers& barriers,
      std::span<std::size_t> rev_lef_ranks, std::span<std::size_t> fwd_lef_ranks,
      std::span<bp_t> rev_moves, std::span<bp_t> fwd_moves, std::span<CollisionT> rev_collisions,
      std::span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng) const
      noexcept(utils::ndebug_defined());

  /// Detect and stall LEFs with one or more extrusion units located at chromosomal boundaries.

  //! \return a pair of numbers consisting in the number of rev units located at the 5'-end and
  //! number of fwd units located at the 3'-end.
  static std::pair<std::size_t, std::size_t> detect_units_at_interval_boundaries(
      const GenomicInterval& interval, std::span<const Lef> lefs,
      std::span<const std::size_t> rev_lef_ranks, std::span<const std::size_t> fwd_lef_ranks,
      std::span<const bp_t> rev_moves, std::span<const bp_t> fwd_moves,
      std::span<CollisionT> rev_collisions, std::span<CollisionT> fwd_collisions);

  /// Detect collisions between LEFs and extrusion barriers.

  //! After calling this function, the entries in \p rev_collisions and \p fwd_collisions
  //! corresponding to reverse or forward extrusion units that will collide with an extrusion
  //! barrier in the current iteration, will be set to the index corresponding to the barrier that
  //! is causing the collision.
  void detect_lef_bar_collisions(
      std::span<const Lef> lefs, std::span<const std::size_t> rev_lef_ranks,
      std::span<const std::size_t> fwd_lef_ranks, std::span<const bp_t> rev_moves,
      std::span<const bp_t> fwd_moves, const ExtrusionBarriers& barriers,
      std::span<CollisionT> rev_collisions, std::span<CollisionT> fwd_collisions,
      random::PRNG_t& rand_eng, std::size_t num_rev_units_at_5prime = 0,
      std::size_t num_fwd_units_at_3prime = 0) const noexcept(utils::ndebug_defined());

  /// Detect collisions between LEFs moving in opposite directions.

  //! Primary LEF-LEF collisions can occur when two extrusion units moving in opposite direction
  //! collide with each other.
  //! After calling this function, the entries in \p rev_collisions and \p fwd_collisions
  //! corresponding to units stalled due to primary LEF-LEF collisions will be set to the index
  //! pointing to the extrusion unit that is causing the collisions.
  //! The index i is encoded as num_barriers + i.
  void detect_primary_lef_lef_collisions(
      std::span<const Lef> lefs, const ExtrusionBarriers& barriers,
      std::span<const std::size_t> rev_lef_ranks, std::span<const std::size_t> fwd_lef_ranks,
      std::span<const bp_t> rev_moves, std::span<const bp_t> fwd_moves,
      std::span<CollisionT> rev_collisions, std::span<CollisionT> fwd_collisions,
      random::PRNG_t& rand_eng, std::size_t num_rev_units_at_5prime = 0,
      std::size_t num_fwd_units_at_3prime = 0) const noexcept(utils::ndebug_defined());

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
  //! The index i is encoded as num_barriers + num_lefs + i.
  // clang-format on
  void process_secondary_lef_lef_collisions(
      const GenomicInterval& interval, std::span<const Lef> lefs,
      std::span<const std::size_t> rev_lef_ranks, std::span<const std::size_t> fwd_lef_ranks,
      std::span<bp_t> rev_moves, std::span<bp_t> fwd_moves, std::span<CollisionT> rev_collisions,
      std::span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng,
      std::size_t num_rev_units_at_5prime = 0, std::size_t num_fwd_units_at_3prime = 0) const
      noexcept(utils::ndebug_defined());

  static void fix_secondary_lef_lef_collisions(
      const GenomicInterval& interval, std::span<Lef> lefs, std::span<std::size_t> rev_lef_ranks,
      std::span<std::size_t> fwd_lef_ranks, std::span<bp_t> rev_moves, std::span<bp_t> fwd_moves,
      std::span<CollisionT> rev_collisions, std::span<CollisionT> fwd_collisions,
      std::size_t num_rev_units_at_5prime,
      std::size_t num_fwd_units_at_3prime) noexcept(utils::ndebug_defined());

  /// Correct moves to comply with the constraints imposed by LEF-BAR collisions.
  static void correct_moves_for_lef_bar_collisions(
      std::span<const Lef> lefs, const ExtrusionBarriers& barriers, std::span<bp_t> rev_moves,
      std::span<bp_t> fwd_moves, std::span<const CollisionT> rev_collisions,
      std::span<const CollisionT> fwd_collisions) noexcept(utils::ndebug_defined());

  /// Correct moves to comply with the constraints imposed by primary LEF-LEF collisions.
  static void correct_moves_for_primary_lef_lef_collisions(
      std::span<const Lef> lefs, std::span<const std::size_t> rev_ranks,
      std::span<const std::size_t> fwd_ranks, std::span<bp_t> rev_moves, std::span<bp_t> fwd_moves,
      std::span<const CollisionT> rev_collisions,
      std::span<const CollisionT> fwd_collisions) noexcept(utils::ndebug_defined());

  /// Register contacts for chromosome \p interval using the position the extrusion units of the
  /// LEFs in \p lefs whose index is present in \p selected_lef_idx.
  /// Return the actual number of contacts registered
  std::size_t register_contacts_loop(GenomicInterval& interval, std::span<const Lef> lefs,
                                     std::size_t num_sampling_events,
                                     random::PRNG_t& rand_eng) const;
  std::size_t register_contacts_loop(bp_t start_pos, bp_t end_pos,
                                     ContactMatrixDense<contacts_t>& contacts,
                                     std::span<const Lef> lefs, std::size_t num_sampling_events,
                                     random::PRNG_t& rand_eng) const;
  std::size_t register_contacts_tad(GenomicInterval& interval, std::span<const Lef> lefs,
                                    std::size_t num_sampling_events,
                                    random::PRNG_t& rand_eng) const;
  std::size_t register_contacts_tad(bp_t start_pos, bp_t end_pos,
                                    ContactMatrixDense<contacts_t>& contacts,
                                    std::span<const Lef> lefs, std::size_t num_sampling_events,
                                    random::PRNG_t& rand_eng) const;

  std::size_t register_1d_lef_occupancy(GenomicInterval& interval, std::span<const Lef> lefs,
                                        std::size_t num_sampling_events,
                                        random::PRNG_t& rand_eng) const;

  std::size_t register_1d_lef_occupancy(bp_t start_pos, bp_t end_pos,
                                        std::vector<std::atomic<std::uint64_t>>& occupancy_buff,
                                        std::span<const Lef> lefs, std::size_t num_sampling_events,
                                        random::PRNG_t& rand_eng) const;

  template <typename MaskT>
  inline static void select_lefs_to_bind(std::span<const Lef> lefs,
                                         MaskT& mask) noexcept(utils::ndebug_defined());

  std::size_t release_lefs(std::span<Lef> lefs, const ExtrusionBarriers& barriers,
                           std::span<const CollisionT> rev_collisions,
                           std::span<const CollisionT> fwd_collisions, random::PRNG_t& rand_eng,
                           bool burnin_completed) const noexcept;

  [[nodiscard]] static std::pair<bp_t /*rev*/, bp_t /*fwd*/> compute_lef_lef_collision_pos(
      const ExtrusionUnit& rev_unit, const ExtrusionUnit& fwd_unit, bp_t rev_move, bp_t fwd_move);

  static void compute_loop_size_stats(std::span<const Lef> lefs,
                                      std::deque<double>& cfx_of_variations,
                                      std::deque<double>& avg_loop_sizes,
                                      std::size_t buff_capacity) noexcept;

  static bool evaluate_burnin(const std::deque<double>& cfx_of_variations,
                              const std::deque<double>& avg_loop_sizes, std::size_t buff_capacity,
                              std::size_t window_size) noexcept;

  void run_burnin(State& s, double lef_binding_rate_burnin) const;

  void sample_and_register_contacts(State& s, std::size_t num_sampling_events) const;
  void dump_stats(std::size_t task_id, std::size_t epoch, std::size_t cell_id, bool burnin,
                  const GenomicInterval& interval, std::span<const Lef> lefs,
                  const ExtrusionBarriers& barriers, std::span<const CollisionT> rev_lef_collisions,
                  std::span<const CollisionT> fwd_lef_collisions,
                  compressed_io::Writer& log_writer) const noexcept(utils::ndebug_defined());

  [[nodiscard]] std::size_t compute_tot_target_epochs(std::size_t nlefs,
                                                      std::size_t npixels) const noexcept;
  [[nodiscard]] std::size_t compute_contacts_per_epoch(std::size_t nlefs) const noexcept;
  [[nodiscard]] std::size_t compute_num_lefs(std::size_t size_bp) const noexcept;

  void print_status_update(const Task& t) const noexcept;
  [[nodiscard]] constexpr bool run_lef_lef_collision_trial(random::PRNG_t& rand_eng) const noexcept;
  [[nodiscard]] constexpr bool run_lef_bar_collision_trial(double pblock,
                                                           random::PRNG_t& rand_eng) const noexcept;

  [[nodiscard]] static std::size_t compute_num_worker_threads(std::size_t num_threads,
                                                              std::size_t num_intervals,
                                                              std::size_t num_cells) noexcept;

#ifdef MODLE_ENABLE_TESTING
 public:
  template <class LefsT, class UsizeT, class MaskT>
  inline void test_bind_lefs(const GenomicInterval& interval, LefsT& lefs, UsizeT& rev_lef_ranks,
                             UsizeT& fwd_lef_ranks, const MaskT& mask, random::PRNG_t& rand_eng,
                             std::size_t current_epoch) {
    modle::Simulation::bind_lefs(interval, lefs, rev_lef_ranks, fwd_lef_ranks, mask, rand_eng,
                                 current_epoch);
  }

  static inline void test_adjust_moves(const GenomicInterval& interval,
                                       const std::span<const Lef> lefs,
                                       const std::span<const std::size_t> rev_lef_ranks,
                                       const std::span<const std::size_t> fwd_lef_ranks,
                                       const std::span<bp_t> rev_moves,
                                       const std::span<bp_t> fwd_moves) {
    Simulation::adjust_moves_of_consecutive_extr_units(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                       rev_moves, fwd_moves);
  }
  static inline void test_adjust_and_clamp_moves(const GenomicInterval& interval,
                                                 const std::span<const Lef> lefs,
                                                 const std::span<const std::size_t> rev_lef_ranks,
                                                 const std::span<const std::size_t> fwd_lef_ranks,
                                                 const std::span<bp_t> rev_moves,
                                                 const std::span<bp_t> fwd_moves) {
    Simulation::adjust_moves_of_consecutive_extr_units(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                       rev_moves, fwd_moves);
    Simulation::clamp_moves(interval, lefs, rev_moves, fwd_moves);
  }

  inline void test_generate_moves(const GenomicInterval& interval, const std::span<const Lef> lefs,
                                  const std::span<const std::size_t> rev_lef_ranks,
                                  const std::span<const std::size_t> fwd_lef_ranks,
                                  const std::span<bp_t> rev_moves, const std::span<bp_t> fwd_moves,
                                  random::PRNG_t& rand_eng, bool adjust_moves_ = false) {
    generate_moves(interval, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, false,
                   rand_eng, adjust_moves_);
  }

  inline static void test_rank_lefs(const std::span<const Lef> lefs,
                                    const std::span<std::size_t> rev_lef_ranks,
                                    const std::span<std::size_t> fwd_lef_ranks,
                                    bool ranks_are_partially_sorted = true,
                                    bool init_buffers = false) {
    Simulation::rank_lefs(lefs, rev_lef_ranks, fwd_lef_ranks, ranks_are_partially_sorted,
                          init_buffers);
  }

  inline static void test_detect_units_at_interval_boundaries(
      const GenomicInterval& interval, const std::span<const Lef> lefs,
      const std::span<const std::size_t> rev_lef_ranks,
      const std::span<const std::size_t> fwd_lef_ranks, const std::span<const bp_t> rev_moves,
      const std::span<const bp_t> fwd_moves, const std::span<CollisionT> rev_collisions,
      const std::span<CollisionT> fwd_collisions) {
    Simulation::detect_units_at_interval_boundaries(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                    rev_moves, fwd_moves, rev_collisions,
                                                    fwd_collisions);
  }

  inline void test_detect_lef_bar_collisions(
      const std::span<const Lef> lefs, const std::span<const std::size_t> rev_lef_ranks,
      const std::span<const std::size_t> fwd_lef_ranks, const std::span<const bp_t> rev_moves,
      const std::span<const bp_t> fwd_moves, const ExtrusionBarriers& barriers,
      const std::span<CollisionT> rev_collisions, const std::span<CollisionT> fwd_collisions,
      random::PRNG_t& rand_eng) {
    Simulation::detect_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                                          barriers, rev_collisions, fwd_collisions, rand_eng);
  }

  inline static void test_correct_moves_for_lef_bar_collisions(
      const std::span<const Lef> lefs, const ExtrusionBarriers& barriers,
      const std::span<bp_t> rev_moves, const std::span<bp_t> fwd_moves,
      const std::span<const CollisionT> rev_collisions,
      const std::span<const CollisionT> fwd_collisions) {
    Simulation::correct_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                     rev_collisions, fwd_collisions);
  }

  inline static void test_adjust_moves_for_primary_lef_lef_collisions(
      std::span<const Lef> lefs, std::span<const std::size_t> rev_ranks,
      std::span<const std::size_t> fwd_ranks, std::span<bp_t> rev_moves, std::span<bp_t> fwd_moves,
      std::span<const CollisionT> rev_collisions, std::span<const CollisionT> fwd_collisions) {
    Simulation::correct_moves_for_primary_lef_lef_collisions(
        lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);
  }

  inline void test_process_collisions(
      const GenomicInterval& interval, const std::span<const Lef> lefs,
      const std::span<const std::size_t> rev_lef_ranks,
      const std::span<const std::size_t> fwd_lef_ranks, const std::span<bp_t> rev_moves,
      const std::span<bp_t> fwd_moves, const ExtrusionBarriers& barriers,
      const std::span<CollisionT> rev_collisions, const std::span<CollisionT> fwd_collisions,
      random::PRNG_t& rand_eng) {
    const auto [num_rev_units_at_5prime, num_fwd_units_at_3prime] =
        Simulation::detect_units_at_interval_boundaries(interval, lefs, rev_lef_ranks,
                                                        fwd_lef_ranks, rev_moves, fwd_moves,
                                                        rev_collisions, fwd_collisions);

    detect_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, barriers,
                              rev_collisions, fwd_collisions, rand_eng, num_rev_units_at_5prime,
                              num_fwd_units_at_3prime);

    detect_primary_lef_lef_collisions(lefs, barriers, rev_lef_ranks, fwd_lef_ranks, rev_moves,
                                      fwd_moves, rev_collisions, fwd_collisions, rand_eng,
                                      num_rev_units_at_5prime, num_fwd_units_at_3prime);

    Simulation::correct_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                     rev_collisions, fwd_collisions);

    Simulation::correct_moves_for_primary_lef_lef_collisions(
        lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);

    process_secondary_lef_lef_collisions(interval, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves,
                                         fwd_moves, rev_collisions, fwd_collisions, rand_eng,
                                         num_rev_units_at_5prime, num_fwd_units_at_3prime);
  }

  static inline void test_fix_secondary_lef_lef_collisions(
      const GenomicInterval& interval, const std::span<Lef> lefs,
      const std::span<std::size_t> rev_lef_ranks, const std::span<std::size_t> fwd_lef_ranks,
      const std::span<bp_t> rev_moves, const std::span<bp_t> fwd_moves,
      const std::span<CollisionT> rev_collisions, const std::span<CollisionT> fwd_collisions) {
    Simulation::fix_secondary_lef_lef_collisions(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                 rev_moves, fwd_moves, rev_collisions,
                                                 fwd_collisions, 0, 0);
  }

  inline void test_process_lef_lef_collisions(
      const GenomicInterval& interval, std::span<const Lef> lefs, const ExtrusionBarriers& barriers,
      std::span<const std::size_t> rev_lef_ranks, std::span<const std::size_t> fwd_lef_ranks,
      std::span<bp_t> rev_moves, std::span<bp_t> fwd_moves, std::span<CollisionT> rev_collisions,
      std::span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng) {
    Simulation::detect_primary_lef_lef_collisions(lefs, barriers, rev_lef_ranks, fwd_lef_ranks,
                                                  rev_moves, fwd_moves, rev_collisions,
                                                  fwd_collisions, rand_eng);

    Simulation::correct_moves_for_primary_lef_lef_collisions(
        lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);

    Simulation::process_secondary_lef_lef_collisions(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                     rev_moves, fwd_moves, rev_collisions,
                                                     fwd_collisions, rand_eng);
  }

  inline void test_detect_primary_lef_lef_collisions(
      std::span<const Lef> lefs, const ExtrusionBarriers& barriers,
      std::span<const std::size_t> rev_lef_ranks, std::span<const std::size_t> fwd_lef_ranks,
      std::span<bp_t> rev_moves, std::span<bp_t> fwd_moves, std::span<CollisionT> rev_collisions,
      std::span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng) {
    Simulation::detect_primary_lef_lef_collisions(lefs, barriers, rev_lef_ranks, fwd_lef_ranks,
                                                  rev_moves, fwd_moves, rev_collisions,
                                                  fwd_collisions, rand_eng);
  }

#endif
};
}  // namespace modle

template <>
struct fmt::formatter<modle::Simulation::Task> {
  inline constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  inline auto format(const modle::Simulation::Task& t, FormatContext& ctx) const
      -> decltype(ctx.out());
};

template <>
struct fmt::formatter<modle::Simulation::State> {
  inline constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  inline auto format(const modle::Simulation::State& s, FormatContext& ctx) const
      -> decltype(ctx.out());
};

#include "../../simulation_impl.hpp"  // IWYU pragma: export
