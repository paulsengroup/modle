// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/container/btree_set.h>            // IWYU pragma: keep for btree_set
#include <absl/container/flat_hash_map.h>        // for flat_hash_map
#include <absl/container/flat_hash_set.h>        // for flat_hash_set
#include <absl/types/span.h>                     // for Span
#include <fmt/format.h>                          // for format_parse_context, formatter
#include <moodycamel/blockingconcurrentqueue.h>  // for BlockingConcurrentQueue
#include <xxh3.h>                                // for XXH_INLINE_XXH3_createState, XXH3...

#include <atomic>                                   // for atomic
#include <boost/circular_buffer.hpp>                // for circular_buffer
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/filesystem/path.hpp>                // for path
#include <cstddef>                                  // for size_t
#include <cstdint>                                  // for uint64_t
#include <deque>                                    // for deque
#include <exception>                                // for exception_ptr
#include <limits>                                   // for numeric_limits
#include <memory>                                   // for allocator, shared_ptr, unique_ptr
#include <mutex>                                    // for mutex
#include <string>                                   // for string
#include <string_view>                              // for string_view
#include <thread_pool/thread_pool.hpp>              // for thread_pool
#include <utility>                                  // for pair
#include <vector>                                   // for vector

#include "modle/bed.hpp"                                // for BED (ptr only), BED_tree
#include "modle/common/common.hpp"                      // for bp_t, collision_t, contacts_t
#include "modle/common/config.hpp"                      // for Config
#include "modle/common/random.hpp"                      // for PRNG_t, normal_distribution, unif...
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...
#include "modle/common/utils.hpp"                       // for ndebug_defined, XXH3_Deleter, XXH...
#include "modle/compressed_io.hpp"                      // for Writer
#include "modle/contacts.hpp"                           // for ContactMatrix
#include "modle/extrusion_barriers.hpp"                 // for ExtrusionBarrier
#include "modle/extrusion_factors.hpp"                  // for Lef, ExtrusionUnit (ptr only)
#include "modle/genome.hpp"                             // for Chromosome (ptr only), Genome

template <typename I, typename T>
class IITree;

namespace modle {

class Simulation : Config {
 public:
  explicit Simulation(const Config& c, bool import_chroms = true);
  std::string to_string() = delete;
  void print() = delete;

  [[nodiscard]] size_t size() const;
  [[nodiscard]] size_t simulated_size() const;

  using lef_move_generator_t = random::normal_distribution<double>;
  using chrom_pos_generator_t = random::uniform_int_distribution<bp_t>;

  static constexpr auto NO_COLLISION = (std::numeric_limits<collision_t>::max)();
  static constexpr auto REACHED_CHROM_BOUNDARY = (std::numeric_limits<collision_t>::max)() - 1;
  static constexpr auto Mbp = 1.0e6;

 private:
  struct BaseTask {  // NOLINT(altera-struct-pack-align)
    size_t id{};
    Chromosome* chrom{};
    size_t cell_id{};
    size_t num_target_epochs{};
    size_t num_target_contacts{};
    size_t num_lefs{};
    absl::Span<const ExtrusionBarrier> barriers{};
  };

 public:
  struct Task : BaseTask {  // NOLINT(altera-struct-pack-align)
    static Task from_string(std::string_view serialized_task, Genome& genome);
  };

  struct TaskPW : BaseTask {  // NOLINT(altera-struct-pack-align)
    TaskPW() = default;
    static TaskPW from_string(std::string_view serialized_task, Genome& genome);

    bp_t deletion_begin{};
    bp_t deletion_size{};
    bp_t window_start{};
    bp_t window_end{};
    bp_t active_window_start{};
    bp_t active_window_end{};

    std::shared_ptr<const ContactMatrix<contacts_t>> reference_contacts{};

    absl::Span<const bed::BED> feats1{};
    absl::Span<const bed::BED> feats2{};
  };

  struct State : BaseTask {  // NOLINT(altera-struct-pack-align)
    State() = default;
    size_t epoch{};
    bool burnin_completed{false};
    size_t num_active_lefs{0};
    size_t num_burnin_epochs{0};
    size_t num_contacts{0};

    random::PRNG_t rand_eng{};
    uint64_t seed{};
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_USED_BUT_MARKED_UNUSED
    std::unique_ptr<XXH3_state_t, utils::XXH3_Deleter> xxh_state{XXH3_createState()};
    DISABLE_WARNING_POP

   protected:
    std::vector<Lef> lef_buff{};
    std::vector<size_t> rank_buff1{};
    std::vector<size_t> rank_buff2{};
    boost::dynamic_bitset<> barrier_mask{};
    std::vector<bp_t> moves_buff1{};
    std::vector<bp_t> moves_buff2{};
    std::vector<size_t> idx_buff{};
    std::vector<collision_t> collision_buff1{};
    std::vector<collision_t> collision_buff2{};
    std::deque<double> cfx_of_variation_buff{};
    std::deque<double> avg_loop_size_buff{};

   public:
    [[nodiscard]] absl::Span<Lef> get_lefs(size_t size = 0) noexcept;
    [[nodiscard]] absl::Span<size_t> get_rev_ranks(size_t size = 0) noexcept;
    [[nodiscard]] absl::Span<size_t> get_fwd_ranks(size_t size = 0) noexcept;
    [[nodiscard]] boost::dynamic_bitset<>& get_barrier_mask() noexcept;
    [[nodiscard]] absl::Span<bp_t> get_rev_moves(size_t size = 0) noexcept;
    [[nodiscard]] absl::Span<bp_t> get_fwd_moves(size_t size = 0) noexcept;
    [[nodiscard]] absl::Span<size_t> get_idx_buff(size_t size = 0) noexcept;
    [[nodiscard]] absl::Span<collision_t> get_rev_collisions(size_t size = 0) noexcept;
    [[nodiscard]] absl::Span<collision_t> get_fwd_collisions(size_t size = 0) noexcept;
    [[nodiscard]] std::deque<double>& get_cfx_of_variation() noexcept;
    [[nodiscard]] std::deque<double>& get_avg_loop_sizes() noexcept;

    [[nodiscard]] absl::Span<const Lef> get_lefs(size_t size = 0) const noexcept;
    [[nodiscard]] absl::Span<const size_t> get_rev_ranks(size_t size = 0) const noexcept;
    [[nodiscard]] absl::Span<const size_t> get_fwd_ranks(size_t size = 0) const noexcept;
    [[nodiscard]] const boost::dynamic_bitset<>& get_barrier_mask() const noexcept;
    [[nodiscard]] absl::Span<const bp_t> get_rev_moves(size_t size = 0) const noexcept;
    [[nodiscard]] absl::Span<const bp_t> get_fwd_moves(size_t size = 0) const noexcept;
    [[nodiscard]] absl::Span<const size_t> get_idx_buff(size_t size = 0) const noexcept;
    [[nodiscard]] absl::Span<const collision_t> get_rev_collisions(size_t size = 0) const noexcept;
    [[nodiscard]] absl::Span<const collision_t> get_fwd_collisions(size_t size = 0) const noexcept;
    [[nodiscard]] const std::deque<double>& get_cfx_of_variation() const noexcept;
    [[nodiscard]] const std::deque<double>& get_avg_loop_sizes() const noexcept;

    [[nodiscard]] bool is_modle_pert_state() const noexcept;
    [[nodiscard]] bool is_modle_sim_state() const noexcept;

    // These fields are specific to modle pert
    bp_t deletion_begin{};
    bp_t deletion_size{};
    bp_t window_start{};
    bp_t window_end{};
    bp_t active_window_start{};
    bp_t active_window_end{};

    absl::Span<const bed::BED> feats1{};
    absl::Span<const bed::BED> feats2{};

    std::shared_ptr<std::mutex> contacts_mtx{nullptr};
    std::shared_ptr<const ContactMatrix<contacts_t>> reference_contacts{nullptr};
    std::shared_ptr<ContactMatrix<contacts_t>> contacts{nullptr};

    std::vector<ExtrusionBarrier> barrier_tmp_buff{};

    State& operator=(const Task& task);
    State& operator=(const TaskPW& task);
    [[nodiscard]] std::string to_string() const noexcept;

    void resize_buffers(size_t size = (std::numeric_limits<size_t>::max)());
    void reset_buffers();
  };

  void run_simulate();
  void run_perturbate();
  void run_replay();

 private:
  Genome _genome{};
  std::atomic<bool> _end_of_simulation{false};
  std::atomic<bool> _exception_thrown{false};
  std::vector<std::exception_ptr> _exceptions{};
  std::mutex _exceptions_mutex{};
  thread_pool _tpool;

  [[nodiscard]] bool ok() const noexcept;

  [[nodiscard]] [[maybe_unused]] thread_pool instantiate_thread_pool() const;
  template <typename I>
  [[nodiscard]] inline static thread_pool instantiate_thread_pool(I nthreads, bool clamp_nthreads);

  /// Simulate loop extrusion using the parameters and buffers passed through \p state
  void simulate_one_cell(State& s) const;

  /// Simulate loop extrusion on a Chromosome window using the parameters and buffers passed through
  /// \p state

  //! This function will call simulate_one_cell multiple times to compute the number of
  //! contacts between a pair of features given a specific barrier configuration.
  //! The different configurations are generated by this function based on the parameters from
  //! \p state
  void simulate_window(State& state, compressed_io::Writer& out_stream, std::mutex& cooler_mutex,
                       bool write_contacts_to_cooler = false) const;

  /// Advance the simulation window by one diagonal width.

  //! Return false if the new window extends past the end of \p chrom
  bool advance_window(TaskPW& base_task, const Chromosome& chrom) const;

  /// Map barrier or features to the window specified by \p base_task.

  //! Return false if the window doesn't have any barrier or is missing one or more type of
  //! features
  static bool map_barriers_to_window(TaskPW& base_task, const Chromosome& chrom);
  static bool map_features_to_window(TaskPW& base_task, const Chromosome& chrom);

  /// Monitor the progress queue and write contacts to a .cool file when all the cells belonging to
  /// a given chromosome have been simulated.

  //! IMPORTANT: this function is meant to be run in a dedicated thread.
  void write_contacts_to_disk(std::deque<std::pair<Chromosome*, size_t>>& progress_queue,
                              std::mutex& progress_queue_mutex);

  /// Worker function used to run an instance of the simulation

  //! Worker function used to consume Simulation::Tasks from a task queue, setup a
  //! Simulation::State, then run an instance of the simulation of a specific chromosome (i.e.
  //! simulate loop extrusion on a single chromosome in a cell).
  //! This function is also responsible for allocating and clearing the buffers used throughout the
  //! simulation.
  //! IMPORTANT: this function is meant to be run in a dedicated thread.
  void simulate_worker(uint64_t tid,
                       moodycamel::BlockingConcurrentQueue<Simulation::Task>& task_queue,
                       std::deque<std::pair<Chromosome*, size_t>>& progress_queue,
                       std::mutex& progress_queue_mutex,
                       size_t task_batch_size = 32);  // NOLINT

  void perturbate_worker(uint64_t tid,
                         moodycamel::BlockingConcurrentQueue<Simulation::TaskPW>& task_queue,
                         const boost::filesystem::path& tmp_output_path, std::mutex& cooler_mutex,
                         size_t task_batch_size = 1);  // NOLINT

  void replay_worker(uint64_t tid,
                     moodycamel::BlockingConcurrentQueue<Simulation::TaskPW>& task_queue,
                     std::mutex& cooler_mutex,
                     size_t task_batch_size = 32);  // NOLINT

  /// Bind inactive LEFs, then sort them by their genomic coordinates.

  //! Bind the LEFs passed through \p lefs whose corresponding entry in \p mask is set to true, then
  //! rank rev and fwd units based on their genomic coordinates.
  //! In order to properly handle ties, the sorting criterion is actually a bit more complex. See
  //! documentation and comments for Simulation::rank_lefs for more details).
  template <typename MaskT>
  inline static void bind_lefs(const Chromosome& chrom, absl::Span<Lef> lefs,
                               absl::Span<size_t> rev_lef_ranks, absl::Span<size_t> fwd_lef_ranks,
                               const MaskT& mask, random::PRNG_t& rand_eng,
                               size_t current_epoch) noexcept(utils::ndebug_defined());

  template <typename MaskT>
  inline static void bind_lefs(bp_t start_pos, bp_t end_pos, absl::Span<Lef> lefs,
                               absl::Span<size_t> rev_lef_ranks, absl::Span<size_t> fwd_lef_ranks,
                               const MaskT& mask, random::PRNG_t& rand_eng,
                               size_t current_epoch) noexcept(utils::ndebug_defined());

  static void select_and_bind_lefs(State& s) noexcept(utils::ndebug_defined());

  // clang-format off
  //! Generate moves in reverse direction

  //! This function ensures that generated moves will never cause an extrusion unit to move past
  //! chrom. boundaries.
  //! The callee must ensure that only extrusion units moving in rev. directions are given as the
  //! \p unit parameter
  //! \return move
  // clang-format on
  [[nodiscard]] static bp_t generate_rev_move(const Chromosome& chrom, const ExtrusionUnit& unit,
                                              double avg_extr_speed, double extr_speed_std,
                                              random::PRNG_t& rand_eng);

  //! Generate moves in forward direction. See doc for Simulation::generate_rev_move for more
  //! details
  [[nodiscard]] static bp_t generate_fwd_move(const Chromosome& chrom, const ExtrusionUnit& unit,
                                              double avg_extr_speed, double extr_speed_std,
                                              random::PRNG_t& rand_eng);
  // clang-format off
  //! Generate moves for all active LEFs

  //! Moves are generated by calling Simulation::generate_rev_move and Simulation::generate_fwd_move.
  //! When adjust_moves_ is true, adjust moves to make consecutive LEFs behave in a more realistic way.
  //! See Simulation::adjust_moves_of_consecutive_extr_units for more details
  // clang-format on
  void generate_moves(const Chromosome& chrom, absl::Span<const Lef> lefs,
                      absl::Span<const size_t> rev_lef_ranks,
                      absl::Span<const size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
                      absl::Span<bp_t> fwd_moves, bool burnin_completed, random::PRNG_t& rand_eng,
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
  static void adjust_moves_of_consecutive_extr_units(
      const Chromosome& chrom, absl::Span<const Lef> lefs, absl::Span<const size_t> rev_lef_ranks,
      absl::Span<const size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
      absl::Span<bp_t> fwd_moves) noexcept(utils::ndebug_defined());

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
  static void extrude(const Chromosome& chrom, absl::Span<Lef> lefs,
                      absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
                      size_t num_rev_units_at_5prime = 0,
                      size_t num_fwd_units_at_3prime = 0) noexcept(utils::ndebug_defined());

  //! This is just a wrapper function used to make sure that process_* functions are always called
  //! in the right order

  //! \return number of rev units at the 5'-end, number of fwd units at the 3'-end
  std::pair<size_t, size_t> process_collisions(
      const Chromosome& chrom, absl::Span<const Lef> lefs,
      absl::Span<const ExtrusionBarrier> barriers, const boost::dynamic_bitset<>& barrier_mask,
      absl::Span<const size_t> rev_lef_ranks, absl::Span<const size_t> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<collision_t> rev_collisions, absl::Span<collision_t> fwd_collisions,
      random::PRNG_t& rand_eng) const noexcept(utils::ndebug_defined());

  /// Detect and stall LEFs with one or more extrusion units located at chromosomal boundaries.

  //! \return a pair of numbers consisting in the number of rev units located at the 5'-end and
  //! number of fwd units located at the 3'-end.
  static std::pair<size_t, size_t> detect_units_at_chrom_boundaries(
      const Chromosome& chrom, absl::Span<const Lef> lefs, absl::Span<const size_t> rev_lef_ranks,
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
                                 absl::Span<collision_t> fwd_collisions, random::PRNG_t& rand_eng,
                                 size_t num_rev_units_at_5prime = 0,
                                 size_t num_fwd_units_at_3prime = 0) const
      noexcept(utils::ndebug_defined());

  /// Detect collisions between LEFs moving in opposite directions.

  //! Primary LEF-LEF collisions can occur when two extrusion units moving in opposite direction
  //! collide with each other.
  //! After calling this function, the entries in \p rev_collisions and \p fwd_collisions
  //! corresponding to units stalled due to primary LEF-LEF collisions will be set to the index
  //! pointing to the extrusion unit that is causing the collisions.
  //! The index i is encoded as num_barriers + i.
  void detect_primary_lef_lef_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<const size_t> rev_lef_ranks, absl::Span<const size_t> fwd_lef_ranks,
      absl::Span<const bp_t> rev_moves, absl::Span<const bp_t> fwd_moves,
      absl::Span<collision_t> rev_collisions, absl::Span<collision_t> fwd_collisions,
      random::PRNG_t& rand_eng, size_t num_rev_units_at_5prime = 0,
      size_t num_fwd_units_at_3prime = 0) const noexcept(utils::ndebug_defined());

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
      const Chromosome& chrom, absl::Span<const Lef> lefs, size_t nbarriers,
      absl::Span<const size_t> rev_lef_ranks, absl::Span<const size_t> fwd_lef_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<collision_t> rev_collisions, absl::Span<collision_t> fwd_collisions,
      random::PRNG_t& rand_eng, size_t num_rev_units_at_5prime = 0,
      size_t num_fwd_units_at_3prime = 0) const noexcept(utils::ndebug_defined());

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
        size_t num_rev_units_at_5prime = 0,
        size_t num_fwd_units_at_3prime = 0) noexcept(utils::ndebug_defined());
  */
  /// Register contacts for chromosome \p chrom using the position the extrusion units of the LEFs
  /// in \p lefs whose index is present in \p selected_lef_idx.
  size_t register_contacts(Chromosome& chrom, absl::Span<const Lef> lefs,
                           absl::Span<const size_t> selected_lef_idx) const;

  size_t register_contacts(bp_t start_pos, bp_t end_pos, ContactMatrix<contacts_t>& contacts,
                           absl::Span<const Lef> lefs,
                           absl::Span<const size_t> selected_lef_idx) const;

  size_t register_contacts_w_randomization(Chromosome& chrom, absl::Span<const Lef> lefs,
                                           absl::Span<const size_t> selected_lef_idx,
                                           random::PRNG_t& rand_eng) const;

  size_t register_contacts_w_randomization(bp_t start_pos, bp_t end_pos,
                                           ContactMatrix<contacts_t>& contacts,
                                           absl::Span<const Lef> lefs,
                                           absl::Span<const size_t> selected_lef_idx,
                                           random::PRNG_t& rand_eng) const;

  template <typename MaskT>
  inline static void select_lefs_to_bind(absl::Span<const Lef> lefs,
                                         MaskT& mask) noexcept(utils::ndebug_defined());

  void generate_lef_unloader_affinities(absl::Span<const Lef> lefs,
                                        absl::Span<const ExtrusionBarrier> barriers,
                                        absl::Span<const collision_t> rev_collisions,
                                        absl::Span<const collision_t> fwd_collisions,
                                        absl::Span<double> lef_unloader_affinity) const
      noexcept(utils::ndebug_defined());

  size_t release_lefs(const absl::Span<Lef> lefs, const absl::Span<const ExtrusionBarrier> barriers,
                      const absl::Span<const collision_t> rev_collisions,
                      const absl::Span<const collision_t> fwd_collisions,
                      random::PRNG_t& rand_eng) const noexcept;

  [[nodiscard]] static std::pair<bp_t /*rev*/, bp_t /*fwd*/> compute_lef_lef_collision_pos(
      const ExtrusionUnit& rev_unit, const ExtrusionUnit& fwd_unit, bp_t rev_move, bp_t fwd_move);

  [[nodiscard]] bed::BED_tree<> import_deletions() const;
  [[nodiscard]] bed::BED_tree<> generate_deletions() const;
  [[nodiscard]] static absl::flat_hash_set<size_t> import_task_filter(
      const boost::filesystem::path& path_to_task_filter);

  [[noreturn]] void rethrow_exceptions() const;
  [[noreturn]] void handle_exceptions();

  static void compute_loop_size_stats(absl::Span<const Lef> lefs,
                                      std::deque<double>& cfx_of_variation_buff,
                                      std::deque<double>& avg_loop_size_buff,
                                      size_t buff_capacity) noexcept;

  static bool evaluate_burnin(const std::deque<double>& cfx_of_variation_buff,
                              const std::deque<double>& avg_loop_size_buff, size_t buff_capacity,
                              size_t window_size) noexcept;

  void run_burnin(State& s, double lef_binding_rate_burnin) const;

  void sample_and_register_contacts(State& s, const double avg_nlefs_to_sample) const;

#ifdef ENABLE_TESTING
 public:
  template <typename MaskT>
  inline void test_bind_lefs(const Chromosome& chrom, const absl::Span<Lef> lefs,
                             const absl::Span<size_t> rev_lef_ranks,
                             const absl::Span<size_t> fwd_lef_ranks, const MaskT& mask,
                             random::PRNG_t& rand_eng, size_t current_epoch) {
    modle::Simulation::bind_lefs(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, mask, rand_eng,
                                 current_epoch);
  }

  static inline void test_adjust_moves(const Chromosome& chrom, const absl::Span<const Lef> lefs,
                                       const absl::Span<const size_t> rev_lef_ranks,
                                       const absl::Span<const size_t> fwd_lef_ranks,
                                       const absl::Span<bp_t> rev_moves,
                                       const absl::Span<bp_t> fwd_moves) {
    Simulation::adjust_moves_of_consecutive_extr_units(chrom, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                       rev_moves, fwd_moves);
  }

  inline void test_generate_moves(const Chromosome& chrom, const absl::Span<const Lef> lefs,
                                  const absl::Span<const size_t> rev_lef_ranks,
                                  const absl::Span<const size_t> fwd_lef_ranks,
                                  const absl::Span<bp_t> rev_moves,
                                  const absl::Span<bp_t> fwd_moves, random::PRNG_t& rand_eng,
                                  bool adjust_moves_ = false) {
    this->generate_moves(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, false,
                         rand_eng, adjust_moves_);
  }

  inline static void test_rank_lefs(const absl::Span<const Lef> lefs,
                                    const absl::Span<size_t> rev_lef_rank_buff,
                                    const absl::Span<size_t> fwd_lef_rank_buff,
                                    bool ranks_are_partially_sorted = true,
                                    bool init_buffers = false) {
    Simulation::rank_lefs(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, ranks_are_partially_sorted,
                          init_buffers);
  }

  inline static void test_detect_units_at_chrom_boundaries(
      const Chromosome& chrom, const absl::Span<const Lef> lefs,
      const absl::Span<const size_t> rev_lef_ranks, const absl::Span<const size_t> fwd_lef_ranks,
      const absl::Span<const bp_t> rev_moves, const absl::Span<const bp_t> fwd_moves,
      const absl::Span<collision_t> rev_collisions, const absl::Span<collision_t> fwd_collisions) {
    Simulation::detect_units_at_chrom_boundaries(chrom, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                 rev_moves, fwd_moves, rev_collisions,
                                                 fwd_collisions);
  }

  inline void test_detect_lef_bar_collisions(
      const absl::Span<const Lef> lefs, const absl::Span<const size_t> rev_lef_rank_buff,
      const absl::Span<const size_t> fwd_lef_rank_buff, const absl::Span<const bp_t> rev_move_buff,
      const absl::Span<const bp_t> fwd_move_buff,
      const absl::Span<const ExtrusionBarrier> extr_barriers,
      const boost::dynamic_bitset<>& barrier_mask, const absl::Span<collision_t> rev_collisions,
      const absl::Span<collision_t> fwd_collisions, random::PRNG_t& rand_eng) {
    Simulation::detect_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff,
                                          fwd_move_buff, extr_barriers, barrier_mask,
                                          rev_collisions, fwd_collisions, rand_eng);
  }

  inline static void test_correct_moves_for_lef_bar_collisions(
      const absl::Span<const Lef> lefs, const absl::Span<const ExtrusionBarrier> barriers,
      const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
      const absl::Span<const collision_t> rev_collisions,
      const absl::Span<const collision_t> fwd_collisions) {
    Simulation::correct_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                     rev_collisions, fwd_collisions);
  }

  inline static void test_adjust_moves_for_primary_lef_lef_collisions(
      absl::Span<const Lef> lefs, absl::Span<const ExtrusionBarrier> barriers,
      absl::Span<const size_t> rev_ranks, absl::Span<const size_t> fwd_ranks,
      absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves,
      absl::Span<const collision_t> rev_collisions, absl::Span<const collision_t> fwd_collisions) {
    Simulation::correct_moves_for_primary_lef_lef_collisions(
        lefs, barriers, rev_ranks, fwd_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);
  }

  inline void test_process_collisions(
      const Chromosome& chrom, const absl::Span<const Lef> lefs,
      const absl::Span<const size_t> rev_lef_rank_buff,
      const absl::Span<const size_t> fwd_lef_rank_buff, const absl::Span<bp_t> rev_move_buff,
      const absl::Span<bp_t> fwd_move_buff, const absl::Span<const ExtrusionBarrier> extr_barriers,
      const boost::dynamic_bitset<>& barrier_mask, const absl::Span<collision_t> rev_collisions,
      const absl::Span<collision_t> fwd_collisions, random::PRNG_t& rand_eng) {
    const auto [num_rev_units_at_5prime, num_fwd_units_at_3prime] =
        Simulation::detect_units_at_chrom_boundaries(chrom, lefs, rev_lef_rank_buff,
                                                     fwd_lef_rank_buff, rev_move_buff,
                                                     fwd_move_buff, rev_collisions, fwd_collisions);

    this->detect_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff,
                                    fwd_move_buff, extr_barriers, barrier_mask, rev_collisions,
                                    fwd_collisions, rand_eng, num_rev_units_at_5prime,
                                    num_fwd_units_at_3prime);

    this->detect_primary_lef_lef_collisions(
        lefs, extr_barriers, rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff, fwd_move_buff,
        rev_collisions, fwd_collisions, rand_eng, num_rev_units_at_5prime, num_fwd_units_at_3prime);

    Simulation::correct_moves_for_lef_bar_collisions(lefs, extr_barriers, rev_move_buff,
                                                     fwd_move_buff, rev_collisions, fwd_collisions);

    Simulation::correct_moves_for_primary_lef_lef_collisions(
        lefs, extr_barriers, rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff, fwd_move_buff,
        rev_collisions, fwd_collisions);

    this->process_secondary_lef_lef_collisions(chrom, lefs, extr_barriers.size(), rev_lef_rank_buff,
                                               fwd_lef_rank_buff, rev_move_buff, fwd_move_buff,
                                               rev_collisions, fwd_collisions, rand_eng,
                                               num_rev_units_at_5prime, num_fwd_units_at_3prime);
  }

  inline void test_process_lef_lef_collisions(
      const Chromosome& chrom, absl::Span<const Lef> lefs,
      absl::Span<const ExtrusionBarrier> extr_barriers, absl::Span<const size_t> rev_lef_rank_buff,
      absl::Span<const size_t> fwd_lef_rank_buff, absl::Span<bp_t> rev_move_buff,
      absl::Span<bp_t> fwd_move_buff, absl::Span<collision_t> rev_collisions,
      absl::Span<collision_t> fwd_collisions, random::PRNG_t& rand_eng) {
    Simulation::detect_primary_lef_lef_collisions(lefs, extr_barriers, rev_lef_rank_buff,
                                                  fwd_lef_rank_buff, rev_move_buff, fwd_move_buff,
                                                  rev_collisions, fwd_collisions, rand_eng);

    Simulation::correct_moves_for_primary_lef_lef_collisions(
        lefs, extr_barriers, rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff, fwd_move_buff,
        rev_collisions, fwd_collisions);

    Simulation::process_secondary_lef_lef_collisions(
        chrom, lefs, extr_barriers.size(), rev_lef_rank_buff, fwd_lef_rank_buff, rev_move_buff,
        fwd_move_buff, rev_collisions, fwd_collisions, rand_eng);
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
  inline auto format(const modle::Simulation::Task& t, FormatContext& ctx) -> decltype(ctx.out());
};

template <>
struct fmt::formatter<modle::Simulation::TaskPW> {
  inline constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  inline auto format(const modle::Simulation::TaskPW& t, FormatContext& ctx) -> decltype(ctx.out());
};

template <>
struct fmt::formatter<modle::Simulation::State> {
  inline constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  inline auto format(const modle::Simulation::State& s, FormatContext& ctx) -> decltype(ctx.out());
};

#include "../../simulation_impl.hpp"  // IWYU pragma: export
