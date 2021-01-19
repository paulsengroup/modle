#pragma once

#ifdef USE_XOSHIRO
#include <XoshiroCpp.hpp>
#endif
#include <cstdint>  // for uint*_t, UINT*_MAX
#include <iosfwd>   // for size_t
#include <limits>
#include <memory>       // for unique_ptr
#include <random>       // for mt19937, geometric_distribution
#include <string_view>  // for string_view
#include <utility>      // for pair

#include "modle/common.hpp"

namespace modle {

#ifdef USE_XOSHIRO
using PRNG = XoshiroCpp::Xoshiro256PlusPlus;
using seeder = XoshiroCpp::SplitMix64;
#else
using PRNG = std::mt19937_64;
using seeder = std::seed_seq;
#endif

struct Chromosome;
class ExtrusionBarrier;
class ExtrusionUnit;

/** \brief Class to model a Loop Extrusion Factor (LEF).
 *
 * A Lef consists of two ExtrusionUnit%s: the Lef::_left_unit always moves in the 5' direction and a
 * Lef::_right_unit, which always moves in the 3' direction (assuming Lef is bound to the positive
 * strand of a DNA molecule).
 *
 * A Lef can bind to the DNA through the Lef::bind_at_pos member function. This function calls
 * Lef::_lifetime_generator to generate the initial lifetime of the Lef-DNA binding. The lifetime is
 * decreased at every call to Lef::extrude. Once the lifetime reaches 0, Lef::unload is called, and
 * this Lef instance detaches from the DNA molecule (or DNA::Bin%s to be more precise) and becomes
 * available for the rebind.
 *
 * The DNA binding lifetime can be extended by a certain amount when both ExtrusionUnit%s are
 * stalled because of an ExtrusionBarrier. The magnitude of the stall and lifetime increase depends
 * on the orientation of the two ExtrusionBarrier%s.
 *
 * Most functions of this class are basically just wrappers around member functions of the
 * ExtrusionUnit class, which are the ones who do the actual the heavy-lifting.
 */
class Lef {
 public:
  inline Lef(uint32_t bin_size, uint32_t avg_processivity, double probability_of_extruder_bypass,
             double unloader_strength_coefficient);

  [[nodiscard]] inline std::string_view get_chr_name() const;
  [[nodiscard]] inline uint32_t get_loop_size() const;
  [[nodiscard]] inline uint32_t get_avg_processivity() const;
  [[nodiscard]] inline const DNA::Bin& get_first_bin() const;
  [[nodiscard]] inline const DNA::Bin& get_last_bin() const;
  [[nodiscard]] inline Chromosome* get_ptr_to_chr();
  [[nodiscard]] inline std::pair<uint32_t, uint32_t> get_pos() const;
  [[nodiscard]] inline double get_probability_of_extr_unit_bypass() const;
  [[nodiscard]] inline std::size_t get_bin_size() const;
  [[nodiscard]] inline std::size_t get_nbins() const;
  [[nodiscard]] inline uint64_t get_tot_bp_extruded() const;
  inline void reset_tot_bp_extruded();

  /// Calls extrude on the ExtrusionUnit%s. Returns the number of bp extruded
  inline uint32_t extrude(modle::PRNG& rand_eng);
  /// Register a contact between the DNA::Bin%s associated with the left and right ExtrusionUnit%s
  inline void register_contact();
  [[nodiscard]] inline std::pair<DNA::Bin*, DNA::Bin*> get_ptr_to_bins();
  [[nodiscard]] inline bool is_bound() const;
  inline void bind_chr_at_random_pos(Chromosome* chr, modle::PRNG& rand_eng,
                                     bool register_contact = false);
  inline void assign_to_chr(Chromosome* chr);
  inline void bind_at_pos(Chromosome* chr, uint32_t pos, modle::PRNG& rand_eng,
                          bool register_contact);
  /** Call ExtrusionUnit::check_constraints on the left and right ExtrusionUnit%s, which in turn
   * check whether there last round of extrusion produced a collision between one of the
   * ExtrusionUnit%s and another instance of ExtrusionUnit, or if the current ExtrusionUnit has
   * reached an ExtrusionBoundary.
   *
   * This function also takes care of extending Lef%'s lifetime where appropriate.
   */
  inline void check_constraints(modle::PRNG& rang_eng);
  inline bool try_rebind(modle::PRNG& rand_eng, double prob_of_rebinding, bool register_contact);
  inline bool try_rebind(modle::PRNG& rand_eng);
  inline std::size_t bind_at_random_pos(modle::PRNG& rand_eng, bool register_contact = false);

 private:
  Chromosome* _chr{nullptr};
  uint32_t _lifetime{0};
  uint32_t _avg_processivity;
  double _probability_of_extr_unit_bypass;
  double _unloader_strength_coeff;
  std::geometric_distribution<uint32_t> _lifetime_generator;
  // This will be used to add an offset of -2, -1, +1, +2 when binding to a bin. This is needed to
  // avoid the undesired zebra pattern described in #3
  // std::discrete_distribution<int8_t> _bin_idx_offset_generator{{1, 1, 0, 1, 1}};
  std::uniform_int_distribution<int64_t> _bin_idx_offset_generator{-2, 2};
  std::uniform_real_distribution<> _prob_of_rebind_generator{0.0, 1.0};
  uint64_t _binding_pos{UINT64_MAX};
  /// I am using a unique_ptr instead of storing the extr. units directly to avoid the cyclic
  /// dependency between dna.hpp and lefs.hpp
  std::unique_ptr<ExtrusionUnit> _left_unit;
  std::unique_ptr<ExtrusionUnit> _right_unit;
  uint64_t _tot_bp_extruded{0};

  /// This function resets the state of the Lef and its ExtrusionUnit%s.
  inline void unload();

  /// This function uses computes the probability of unloading given the average LEF processivity,
  /// bin size and the number of active extrusion units.
  [[nodiscard]] inline double compute_prob_of_unloading(uint32_t bin_size,
                                                        uint8_t n_of_active_extr_units = 2) const;
};

class ExtrusionUnit {
  friend class Lef;

 public:
  inline explicit ExtrusionUnit(Lef& lef, double prob_of_extr_unit_bypass);
  [[nodiscard]] inline std::size_t get_pos() const;
  [[nodiscard]] inline dna::Direction get_extr_direction() const;
  [[nodiscard]] inline bool is_stalled() const;
  [[nodiscard]] inline bool is_bound() const;
  inline uint64_t check_constraints(modle::PRNG& rand_eng);
  inline bool try_extrude(modle::PRNG& rand_eng);
  [[nodiscard]] inline double get_prob_of_extr_unit_bypass() const;
  [[nodiscard]] inline std::size_t get_bin_index() const;
  [[nodiscard]] inline std::size_t get_bin_size() const;
  [[nodiscard]] inline DNA::Bin& get_bin();
  [[nodiscard]] inline const DNA::Bin& get_bin() const;

 private:
  Lef& _parent_lef;
  std::size_t _bin_idx{std::numeric_limits<std::size_t>::max()};
  DNA* _dna{nullptr};
  ExtrusionBarrier* _blocking_barrier{nullptr};
  dna::Direction _direction{dna::Direction::none};
  uint32_t _stalls_left{0};
  std::geometric_distribution<uint32_t> _n_stall_generator;

  inline void set_stalls(uint32_t n);
  inline void increment_stalls(uint32_t n = 1);
  inline void decrement_stalls(uint32_t n = 1);
  inline void reset_stalls();
  inline void unload();
  inline void bind(Chromosome* chr, uint32_t pos, dna::Direction direction, modle::PRNG& rand_eng);

  inline uint32_t check_for_extruder_collisions(modle::PRNG& rang_eng);
  [[nodiscard]] inline uint64_t check_for_extrusion_barrier(modle::PRNG& rang_eng);
  inline bool try_moving_to_next_bin();
  inline bool try_moving_to_prev_bin();
  [[nodiscard]] inline bool hard_stall() const;
};

}  // namespace modle

#include "../../lefs_impl.hpp"
