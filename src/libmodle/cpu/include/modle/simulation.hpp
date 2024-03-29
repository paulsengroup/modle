// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/container/flat_hash_set.h>  // for flat_hash_set
#include <absl/types/span.h>               // for Span
#include <fmt/format.h>                    // for format_parse_context, formatter

#include <atomic>       // for atomic
#include <deque>        // for deque
#include <exception>    // for exception_ptr
#include <filesystem>   // for path
#include <limits>       // for numeric_limits
#include <memory>       // for shared_ptr, allocator, unique_ptr
#include <mutex>        // for mutex
#include <string>       // for string
#include <string_view>  // for string_view
#include <utility>      // for pair
#include <vector>       // for vector

#include "modle/bed/bed.hpp"                            // for BED (ptr only), BED_tree
#include "modle/collision_encoding.hpp"                 // for Collision<>
#include "modle/common/common.hpp"                      // for usize, bp_t, contacts_t
#include "modle/common/random.hpp"                      // for PRNG_t, normal_distribution, unif...
#include "modle/common/simulation_config.hpp"           // for Config
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...
#include "modle/common/utils.hpp"                       // for ndebug_defined
#include "modle/contact_matrix_dense.hpp"               // for ContactMatrixDense
#include "modle/context_manager.hpp"
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier
#include "modle/extrusion_factors.hpp"   // for Lef, ExtrusionUnit (ptr only)
#include "modle/genome.hpp"              // for GenomicInterval (ptr only), Genome

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

  [[nodiscard]] usize size() const;
  [[nodiscard]] usize simulated_size() const;

  using lef_move_generator_t = random::normal_distribution<double>;
  using chrom_pos_generator_t = random::uniform_int_distribution<bp_t>;

  static constexpr auto Mbp = 1.0e6;

  using CollisionT = Collision<u32f>;
  struct Task {  // NOLINT(altera-struct-pack-align)
    enum class Status : u8f { PENDING, RUNNING, COMPLETED, FAILED };
    usize id{};
    GenomicInterval* interval{};
    usize cell_id{};
    usize num_target_epochs{};
    usize num_target_contacts{};
    usize num_lefs{};
    random::PRNG_t rand_eng{};
    Status status{Status::PENDING};
  };

  struct State : Task {  // NOLINT(altera-struct-pack-align)
    State() = default;
    usize epoch{};                 // NOLINT
    bool burnin_completed{false};  // NOLINT
    usize num_active_lefs{0};      // NOLINT

    usize num_burnin_epochs{0};  // NOLINT
    usize num_contacts{0};       // NOLINT

    u64 seed{};                                                          // NOLINT
    std::unique_ptr<compressed_io::Writer> model_state_logger{nullptr};  // NOLINT

    ExtrusionBarriers barriers{};

   protected:
    std::vector<Lef> lef_buff{};                 // NOLINT
    std::vector<usize> rank_buff1{};             // NOLINT
    std::vector<usize> rank_buff2{};             // NOLINT
    ExtrusionBarriers _barriers{};               // NOLINT
    std::vector<bp_t> moves_buff1{};             // NOLINT
    std::vector<bp_t> moves_buff2{};             // NOLINT
    std::vector<usize> idx_buff{};               // NOLINT
    std::vector<CollisionT> collision_buff1{};   // NOLINT
    std::vector<CollisionT> collision_buff2{};   // NOLINT
    std::deque<double> cfx_of_variation_buff{};  // NOLINT
    std::deque<double> avg_loop_size_buff{};     // NOLINT

    static constexpr usize npos = absl::Span<usize>::npos;

   public:
    [[nodiscard]] absl::Span<Lef> get_lefs(usize size = npos) noexcept;
    [[nodiscard]] absl::Span<usize> get_rev_ranks(usize size = npos) noexcept;
    [[nodiscard]] absl::Span<usize> get_fwd_ranks(usize size = npos) noexcept;
    [[nodiscard]] absl::Span<bp_t> get_rev_moves(usize size = npos) noexcept;
    [[nodiscard]] absl::Span<bp_t> get_fwd_moves(usize size = npos) noexcept;
    [[nodiscard]] absl::Span<usize> get_idx_buff(usize size = npos) noexcept;
    [[nodiscard]] absl::Span<CollisionT> get_rev_collisions(usize size = npos) noexcept;
    [[nodiscard]] absl::Span<CollisionT> get_fwd_collisions(usize size = npos) noexcept;
    [[nodiscard]] std::deque<double>& get_cfx_of_variation() noexcept;
    [[nodiscard]] std::deque<double>& get_avg_loop_sizes() noexcept;

    [[nodiscard]] absl::Span<const Lef> get_lefs(usize size = npos) const noexcept;
    [[nodiscard]] absl::Span<const usize> get_rev_ranks(usize size = npos) const noexcept;
    [[nodiscard]] absl::Span<const usize> get_fwd_ranks(usize size = npos) const noexcept;
    [[nodiscard]] absl::Span<const bp_t> get_rev_moves(usize size = npos) const noexcept;
    [[nodiscard]] absl::Span<const bp_t> get_fwd_moves(usize size = npos) const noexcept;
    [[nodiscard]] absl::Span<const usize> get_idx_buff(usize size = npos) const noexcept;
    [[nodiscard]] absl::Span<const CollisionT> get_rev_collisions(usize size = npos) const noexcept;
    [[nodiscard]] absl::Span<const CollisionT> get_fwd_collisions(usize size = npos) const noexcept;
    [[nodiscard]] const std::deque<double>& get_cfx_of_variation() const noexcept;
    [[nodiscard]] const std::deque<double>& get_avg_loop_sizes() const noexcept;

    std::shared_ptr<std::mutex> contacts_mtx{nullptr};                  // NOLINT
    std::shared_ptr<ContactMatrixDense<contacts_t>> contacts{nullptr};  // NOLINT

    State& operator=(const Task& task);

    void resize_buffers(usize size = (std::numeric_limits<usize>::max)());
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

  void spawn_worker_threads(usize batch_size);
  void spawn_io_threads();
  void run_simulate();

 private:
  /// Simulate loop extrusion using the parameters and buffers passed through \p state
  void simulate_one_cell(u64 tid, State& s) const;

  //! IMPORTANT: this function is meant to be run in a dedicated thread.
  void simulate_io(std::chrono::milliseconds wait_time = std::chrono::milliseconds(50));

  /// Worker function used to run an instance of the simulation

  //! Worker function used to consume Simulation::Tasks from a task queue, setup a
  //! Simulation::State, then run an instance of the simulation of a specific genomic interval (i.e.
  //! simulate loop extrusion on a single chromosome in a cell).
  //! This function is also responsible for allocating and clearing the buffers used throughout the
  //! simulation.
  //! IMPORTANT: this function is meant to be run in a dedicated thread.
  void simulate_worker(u64 tid, usize task_batch_size);

  /// Bind inactive LEFs, then sort them by their genomic coordinates.

  //! Bind the LEFs passed through \p lefs whose corresponding entry in \p mask is set to true, then
  //! rank rev and fwd units based on their genomic coordinates.
  //! In order to properly handle ties, the sorting criterion is actually a bit more complex. See
  //! documentation and comments for Simulation::rank_lefs for more details).
  template <typename MaskT>
  inline static void bind_lefs(const GenomicInterval& interval, absl::Span<Lef> lefs,
                               absl::Span<usize> rev_lef_ranks, absl::Span<usize> fwd_lef_ranks,
                               const MaskT& mask, random::PRNG_t& rand_eng,
                               usize current_epoch) noexcept(utils::ndebug_defined());

  template <typename MaskT>
  inline static void bind_lefs(bp_t start_pos, bp_t end_pos, absl::Span<Lef> lefs,
                               absl::Span<usize> rev_lef_ranks, absl::Span<usize> fwd_lef_ranks,
                               const MaskT& mask, random::PRNG_t& rand_eng,
                               usize current_epoch) noexcept(utils::ndebug_defined());

  static void select_and_bind_lefs(State& s) noexcept(utils::ndebug_defined());

  // clang-format off
  //! Generate moves for all active LEFs

  //! Moves are generated by calling Simulation::generate_rev_move and Simulation::generate_fwd_move.
  //! When adjust_moves_ is true, adjust moves to make consecutive LEFs behave in a more realistic way.
  //! See Simulation::adjust_moves_of_consecutive_extr_units for more details
  // clang-format on
  void generate_moves(const GenomicInterval& interval, absl::Span<const Lef> lefs,
                      absl::Span<const usize> rev_lef_ranks, absl::Span<const usize> fwd_lef_ranks,
                      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, bool burnin_completed,
                      random::PRNG_t& rand_eng, bool adjust_moves_ = true) const
      noexcept(utils::ndebug_defined());

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
      const GenomicInterval& interval, absl::Span<const Lef> lefs,
      absl::Span<const usize> rev_lef_ranks, absl::Span<const usize> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves) noexcept(utils::ndebug_defined());

  /// Clamp moves to prevent LEFs from falling off chromosomal boundaries
  static void clamp_moves(const GenomicInterval& interval, absl::Span<const Lef> lefs,
                          absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves) noexcept;

  /// Sort extr. units by index based on their position whilst properly dealing with ties. See
  /// comments in \p simulation_impl.hpp file for more details.
  static void rank_lefs(absl::Span<const Lef> lefs, absl::Span<usize> rev_lef_ranks,
                        absl::Span<usize> fwd_lef_ranks, bool ranks_are_partially_sorted = true,
                        bool init_buffers = false) noexcept(utils::ndebug_defined());

  /// Extrude LEFs by applying the respective moves

  //! This function assumes that the move vectors that are passed as argument have already been
  //! processed by:
  //! - Simulation::adjust_moves_of_consecutive_extr_units
  //! - Simulation::process_collisions
  static void extrude(const GenomicInterval& interval, absl::Span<Lef> lefs,
                      absl::Span<const bp_t> rev_moves,
                      absl::Span<const bp_t> fwd_moves) noexcept(utils::ndebug_defined());

  //! This is just a wrapper function used to make sure that process_* functions are always called
  //! in the right order

  //! \return number of rev units at the 5'-end, number of fwd units at the 3'-end
  std::pair<usize, usize> process_collisions(
      const GenomicInterval& interval, absl::Span<Lef> lefs, ExtrusionBarriers& barriers,
      absl::Span<usize> rev_lef_ranks, absl::Span<usize> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
      absl::Span<bp_t> fwd_moves, absl::Span<CollisionT> rev_collisions,
      absl::Span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng) const
      noexcept(utils::ndebug_defined());

  /// Detect and stall LEFs with one or more extrusion units located at chromosomal boundaries.

  //! \return a pair of numbers consisting in the number of rev units located at the 5'-end and
  //! number of fwd units located at the 3'-end.
  static std::pair<usize, usize> detect_units_at_interval_boundaries(
      const GenomicInterval& interval, absl::Span<const Lef> lefs,
      absl::Span<const usize> rev_lef_ranks, absl::Span<const usize> fwd_lef_ranks,
      absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
      absl::Span<CollisionT> rev_collisions, absl::Span<CollisionT> fwd_collisions);

  /// Detect collisions between LEFs and extrusion barriers.

  //! After calling this function, the entries in \p rev_collisions and \p fwd_collisions
  //! corresponding to reverse or forward extrusion units that will collide with an extrusion
  //! barrier in the current iteration, will be set to the index corresponding to the barrier that
  //! is causing the collision.
  void detect_lef_bar_collisions(absl::Span<const Lef> lefs, absl::Span<const usize> rev_lef_ranks,
                                 absl::Span<const usize> fwd_lef_ranks,
                                 absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
                                 const ExtrusionBarriers& barriers,
                                 absl::Span<CollisionT> rev_collisions,
                                 absl::Span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng,
                                 usize num_rev_units_at_5prime = 0,
                                 usize num_fwd_units_at_3prime = 0) const
      noexcept(utils::ndebug_defined());

  /// Detect collisions between LEFs moving in opposite directions.

  //! Primary LEF-LEF collisions can occur when two extrusion units moving in opposite direction
  //! collide with each other.
  //! After calling this function, the entries in \p rev_collisions and \p fwd_collisions
  //! corresponding to units stalled due to primary LEF-LEF collisions will be set to the index
  //! pointing to the extrusion unit that is causing the collisions.
  //! The index i is encoded as num_barriers + i.
  void detect_primary_lef_lef_collisions(
      absl::Span<const Lef> lefs, const ExtrusionBarriers& barriers,
      absl::Span<const usize> rev_lef_ranks, absl::Span<const usize> fwd_lef_ranks,
      absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
      absl::Span<CollisionT> rev_collisions, absl::Span<CollisionT> fwd_collisions,
      random::PRNG_t& rand_eng, usize num_rev_units_at_5prime = 0,
      usize num_fwd_units_at_3prime = 0) const noexcept(utils::ndebug_defined());

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
      const GenomicInterval& interval, absl::Span<const Lef> lefs,
      absl::Span<const usize> rev_lef_ranks, absl::Span<const usize> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<CollisionT> rev_collisions,
      absl::Span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng,
      usize num_rev_units_at_5prime = 0, usize num_fwd_units_at_3prime = 0) const
      noexcept(utils::ndebug_defined());

  static void fix_secondary_lef_lef_collisions(
      const GenomicInterval& interval, absl::Span<Lef> lefs, absl::Span<usize> rev_lef_ranks,
      absl::Span<usize> fwd_lef_ranks, absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<CollisionT> rev_collisions, absl::Span<CollisionT> fwd_collisions,
      usize num_rev_units_at_5prime,
      usize num_fwd_units_at_3prime) noexcept(utils::ndebug_defined());

  /// Correct moves to comply with the constraints imposed by LEF-BAR collisions.
  static void correct_moves_for_lef_bar_collisions(
      absl::Span<const Lef> lefs, const ExtrusionBarriers& barriers, absl::Span<bp_t> rev_moves,
      absl::Span<bp_t> fwd_moves, absl::Span<const CollisionT> rev_collisions,
      absl::Span<const CollisionT> fwd_collisions) noexcept(utils::ndebug_defined());

  /// Correct moves to comply with the constraints imposed by primary LEF-LEF collisions.
  static void correct_moves_for_primary_lef_lef_collisions(
      absl::Span<const Lef> lefs, absl::Span<const usize> rev_ranks,
      absl::Span<const usize> fwd_ranks, absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<const CollisionT> rev_collisions,
      absl::Span<const CollisionT> fwd_collisions) noexcept(utils::ndebug_defined());

  /// Register contacts for chromosome \p interval using the position the extrusion units of the
  /// LEFs in \p lefs whose index is present in \p selected_lef_idx.
  /// Return the actual number of contacts registered
  usize register_contacts_loop(GenomicInterval& interval, absl::Span<const Lef> lefs,
                               usize num_sampling_events, random::PRNG_t& rand_eng) const;
  usize register_contacts_loop(bp_t start_pos, bp_t end_pos,
                               ContactMatrixDense<contacts_t>& contacts, absl::Span<const Lef> lefs,
                               usize num_sampling_events, random::PRNG_t& rand_eng) const;
  usize register_contacts_tad(GenomicInterval& interval, absl::Span<const Lef> lefs,
                              usize num_sampling_events, random::PRNG_t& rand_eng) const;
  usize register_contacts_tad(bp_t start_pos, bp_t end_pos,
                              ContactMatrixDense<contacts_t>& contacts, absl::Span<const Lef> lefs,
                              usize num_sampling_events, random::PRNG_t& rand_eng) const;

  usize register_1d_lef_occupancy(GenomicInterval& interval, absl::Span<const Lef> lefs,
                                  usize num_sampling_events, random::PRNG_t& rand_eng) const;

  usize register_1d_lef_occupancy(bp_t start_pos, bp_t end_pos,
                                  std::vector<std::atomic<u64>>& occupancy_buff,
                                  absl::Span<const Lef> lefs, usize num_sampling_events,
                                  random::PRNG_t& rand_eng) const;

  template <typename MaskT>
  inline static void select_lefs_to_bind(absl::Span<const Lef> lefs,
                                         MaskT& mask) noexcept(utils::ndebug_defined());

  usize release_lefs(absl::Span<Lef> lefs, const ExtrusionBarriers& barriers,
                     absl::Span<const CollisionT> rev_collisions,
                     absl::Span<const CollisionT> fwd_collisions, random::PRNG_t& rand_eng,
                     bool burnin_completed) const noexcept;

  [[nodiscard]] static std::pair<bp_t /*rev*/, bp_t /*fwd*/> compute_lef_lef_collision_pos(
      const ExtrusionUnit& rev_unit, const ExtrusionUnit& fwd_unit, bp_t rev_move, bp_t fwd_move);

  static void compute_loop_size_stats(absl::Span<const Lef> lefs,
                                      std::deque<double>& cfx_of_variations,
                                      std::deque<double>& avg_loop_sizes,
                                      usize buff_capacity) noexcept;

  static bool evaluate_burnin(const std::deque<double>& cfx_of_variations,
                              const std::deque<double>& avg_loop_sizes, usize buff_capacity,
                              usize window_size) noexcept;

  void run_burnin(State& s, double lef_binding_rate_burnin) const;

  void sample_and_register_contacts(State& s, usize num_sampling_events) const;
  void dump_stats(usize task_id, usize epoch, usize cell_id, bool burnin,
                  const GenomicInterval& interval, absl::Span<const Lef> lefs,
                  const ExtrusionBarriers& barriers,
                  absl::Span<const CollisionT> rev_lef_collisions,
                  absl::Span<const CollisionT> fwd_lef_collisions,
                  compressed_io::Writer& log_writer) const noexcept(utils::ndebug_defined());

  [[nodiscard]] usize compute_tot_target_epochs(usize nlefs, usize npixels) const noexcept;
  [[nodiscard]] usize compute_contacts_per_epoch(usize nlefs) const noexcept;
  [[nodiscard]] usize compute_num_lefs(usize size_bp) const noexcept;

  void print_status_update(const Task& t) const noexcept;
  [[nodiscard]] constexpr bool run_lef_lef_collision_trial(random::PRNG_t& rand_eng) const noexcept;
  [[nodiscard]] constexpr bool run_lef_bar_collision_trial(double pblock,
                                                           random::PRNG_t& rand_eng) const noexcept;

  [[nodiscard]] static usize compute_num_worker_threads(usize num_threads, usize num_intervals,
                                                        usize num_cells) noexcept;

#ifdef MODLE_ENABLE_TESTING
 public:
  template <class LefsT, class UsizeT, class MaskT>
  inline void test_bind_lefs(const GenomicInterval& interval, LefsT& lefs, UsizeT& rev_lef_ranks,
                             UsizeT& fwd_lef_ranks, const MaskT& mask, random::PRNG_t& rand_eng,
                             usize current_epoch) {
    modle::Simulation::bind_lefs(interval, absl::MakeSpan(lefs), absl::MakeSpan(rev_lef_ranks),
                                 absl::MakeSpan(fwd_lef_ranks), mask, rand_eng, current_epoch);
  }

  static inline void test_adjust_moves(const GenomicInterval& interval,
                                       const absl::Span<const Lef> lefs,
                                       const absl::Span<const usize> rev_lef_ranks,
                                       const absl::Span<const usize> fwd_lef_ranks,
                                       const absl::Span<bp_t> rev_moves,
                                       const absl::Span<bp_t> fwd_moves) {
    Simulation::adjust_moves_of_consecutive_extr_units(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                       rev_moves, fwd_moves);
  }
  static inline void test_adjust_and_clamp_moves(const GenomicInterval& interval,
                                                 const absl::Span<const Lef> lefs,
                                                 const absl::Span<const usize> rev_lef_ranks,
                                                 const absl::Span<const usize> fwd_lef_ranks,
                                                 const absl::Span<bp_t> rev_moves,
                                                 const absl::Span<bp_t> fwd_moves) {
    Simulation::adjust_moves_of_consecutive_extr_units(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                       rev_moves, fwd_moves);
    Simulation::clamp_moves(interval, lefs, absl::MakeSpan(rev_moves), absl::MakeSpan(fwd_moves));
  }

  inline void test_generate_moves(const GenomicInterval& interval, const absl::Span<const Lef> lefs,
                                  const absl::Span<const usize> rev_lef_ranks,
                                  const absl::Span<const usize> fwd_lef_ranks,
                                  const absl::Span<bp_t> rev_moves,
                                  const absl::Span<bp_t> fwd_moves, random::PRNG_t& rand_eng,
                                  bool adjust_moves_ = false) {
    this->generate_moves(interval, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, false,
                         rand_eng, adjust_moves_);
  }

  inline static void test_rank_lefs(const absl::Span<const Lef> lefs,
                                    const absl::Span<usize> rev_lef_ranks,
                                    const absl::Span<usize> fwd_lef_ranks,
                                    bool ranks_are_partially_sorted = true,
                                    bool init_buffers = false) {
    Simulation::rank_lefs(lefs, rev_lef_ranks, fwd_lef_ranks, ranks_are_partially_sorted,
                          init_buffers);
  }

  inline static void test_detect_units_at_interval_boundaries(
      const GenomicInterval& interval, const absl::Span<const Lef> lefs,
      const absl::Span<const usize> rev_lef_ranks, const absl::Span<const usize> fwd_lef_ranks,
      const absl::Span<const bp_t> rev_moves, const absl::Span<const bp_t> fwd_moves,
      const absl::Span<CollisionT> rev_collisions, const absl::Span<CollisionT> fwd_collisions) {
    Simulation::detect_units_at_interval_boundaries(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                    rev_moves, fwd_moves, rev_collisions,
                                                    fwd_collisions);
  }

  inline void test_detect_lef_bar_collisions(
      const absl::Span<const Lef> lefs, const absl::Span<const usize> rev_lef_ranks,
      const absl::Span<const usize> fwd_lef_ranks, const absl::Span<const bp_t> rev_moves,
      const absl::Span<const bp_t> fwd_moves, const ExtrusionBarriers& barriers,
      const absl::Span<CollisionT> rev_collisions, const absl::Span<CollisionT> fwd_collisions,
      random::PRNG_t& rand_eng) {
    Simulation::detect_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                                          barriers, rev_collisions, fwd_collisions, rand_eng);
  }

  inline static void test_correct_moves_for_lef_bar_collisions(
      const absl::Span<const Lef> lefs, const ExtrusionBarriers& barriers,
      const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
      const absl::Span<const CollisionT> rev_collisions,
      const absl::Span<const CollisionT> fwd_collisions) {
    Simulation::correct_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                     rev_collisions, fwd_collisions);
  }

  inline static void test_adjust_moves_for_primary_lef_lef_collisions(
      absl::Span<const Lef> lefs, absl::Span<const usize> rev_ranks,
      absl::Span<const usize> fwd_ranks, absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<const CollisionT> rev_collisions, absl::Span<const CollisionT> fwd_collisions) {
    Simulation::correct_moves_for_primary_lef_lef_collisions(
        lefs, rev_ranks, fwd_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);
  }

  inline void test_process_collisions(
      const GenomicInterval& interval, const absl::Span<const Lef> lefs,
      const absl::Span<const usize> rev_lef_ranks, const absl::Span<const usize> fwd_lef_ranks,
      const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
      const ExtrusionBarriers& barriers, const absl::Span<CollisionT> rev_collisions,
      const absl::Span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng) {
    const auto [num_rev_units_at_5prime, num_fwd_units_at_3prime] =
        Simulation::detect_units_at_interval_boundaries(interval, lefs, rev_lef_ranks,
                                                        fwd_lef_ranks, rev_moves, fwd_moves,
                                                        rev_collisions, fwd_collisions);

    this->detect_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                                    barriers, rev_collisions, fwd_collisions, rand_eng,
                                    num_rev_units_at_5prime, num_fwd_units_at_3prime);

    this->detect_primary_lef_lef_collisions(lefs, barriers, rev_lef_ranks, fwd_lef_ranks, rev_moves,
                                            fwd_moves, rev_collisions, fwd_collisions, rand_eng,
                                            num_rev_units_at_5prime, num_fwd_units_at_3prime);

    Simulation::correct_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                     rev_collisions, fwd_collisions);

    Simulation::correct_moves_for_primary_lef_lef_collisions(
        lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);

    this->process_secondary_lef_lef_collisions(
        interval, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, rev_collisions,
        fwd_collisions, rand_eng, num_rev_units_at_5prime, num_fwd_units_at_3prime);
  }

  static inline void test_fix_secondary_lef_lef_collisions(
      const GenomicInterval& interval, const absl::Span<Lef> lefs,
      const absl::Span<usize> rev_lef_ranks, const absl::Span<usize> fwd_lef_ranks,
      const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
      const absl::Span<CollisionT> rev_collisions, const absl::Span<CollisionT> fwd_collisions) {
    Simulation::fix_secondary_lef_lef_collisions(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                 rev_moves, fwd_moves, rev_collisions,
                                                 fwd_collisions, 0, 0);
  }

  inline void test_process_lef_lef_collisions(
      const GenomicInterval& interval, absl::Span<const Lef> lefs,
      const ExtrusionBarriers& barriers, absl::Span<const usize> rev_lef_ranks,
      absl::Span<const usize> fwd_lef_ranks, absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<CollisionT> rev_collisions, absl::Span<CollisionT> fwd_collisions,
      random::PRNG_t& rand_eng) {
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
      absl::Span<const Lef> lefs, const ExtrusionBarriers& barriers,
      absl::Span<const usize> rev_lef_ranks, absl::Span<const usize> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves, absl::Span<CollisionT> rev_collisions,
      absl::Span<CollisionT> fwd_collisions, random::PRNG_t& rand_eng) {
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
