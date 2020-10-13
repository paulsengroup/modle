#pragma once

#include <cstdint>
#include <memory>
#include <random>

#include "modle/contacts.hpp"
#include "modle/dna.hpp"
#include "modle/extr_barrier.hpp"

namespace modle {
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
  Lef(uint32_t bin_size, uint32_t avg_processivity, double probability_of_extruder_bypass);

  [[nodiscard]] std::string_view get_chr_name() const;
  [[nodiscard]] uint32_t get_loop_size() const;
  [[nodiscard]] uint32_t get_avg_processivity() const;
  [[nodiscard]] const DNA::Bin& get_first_bin() const;
  [[nodiscard]] const DNA::Bin& get_last_bin() const;
  [[nodiscard]] Chromosome* get_ptr_to_chr();
  [[nodiscard]] std::pair<uint32_t, uint32_t> get_pos() const;
  [[nodiscard]] double get_probability_of_extr_unit_bypass() const;

  /// Calls extrude on the ExtrusionUnit%s. Returns 0, 1 or 2 depending on how many DNA::Bin%s have
  /// been extruded.
  uint32_t extrude();
  /// Register a contact between the DNA::Bin%s associated with the left and right ExtrusionUnit%s
  void register_contact();
  [[nodiscard]] std::pair<DNA::Bin*, DNA::Bin*> get_ptr_to_bins();
  [[nodiscard]] bool is_bound() const;
  void bind_at_pos(Chromosome& chr, uint32_t pos, std::mt19937& rand_gen);
  /** Call ExtrusionUnit::check_constraints on the left and right ExtrusionUnit%s, which in turn
   * check whether there last round of extrusion produced a collision between one of the
   * ExtrusionUnit%s and another instance of ExtrusionUnit, or if the current ExtrusionUnit has
   * reached an ExtrusionBoundary.
   *
   * This function also takes care of extending Lef%'s lifetime where appropriate.
   */
  void check_constraints(std::mt19937& rand_gen);
  bool try_rebind(Chromosome& chr, std::mt19937& rand_eng, double prob_of_rebinding);

 private:
  Chromosome* _chr{nullptr};
  uint32_t _lifetime{0};
  uint32_t _avg_processivity;
  double _probability_of_extr_unit_bypass;
  std::geometric_distribution<uint32_t> _lifetime_generator;
  /// I am using a unique_ptr instead of storing the extr. units directly to avoid the cyclic
  /// dependency between dna.hpp and lefs.hpp
  std::unique_ptr<ExtrusionUnit> _left_unit;
  std::unique_ptr<ExtrusionUnit> _right_unit;

  /// This function resets the state of the Lef and its ExtrusionUnit%s.
  void unload();

  /// This function uses computes the probability of unloading given the average LEF processivity,
  /// bin size and the number of active extrusion units.
  [[nodiscard]] double compute_prob_of_unloading(uint32_t bin_size,
                                                 uint8_t n_of_active_extr_units = 2) const;
};

class ExtrusionUnit {
  friend class Lef;

 public:
  explicit ExtrusionUnit(Lef& lef, double prob_of_extr_unit_bypass);
  [[nodiscard]] uint32_t get_pos() const;
  [[nodiscard]] DNA::Direction get_extr_direction() const;
  [[nodiscard]] bool is_stalled() const;
  [[nodiscard]] bool is_bound() const;
  [[nodiscard]] uint64_t check_constraints(std::mt19937& rand_gen);
  bool try_extrude();
  [[nodiscard]] double get_prob_of_extr_unit_bypass() const;

 private:
  Lef& _parent_lef;
  DNA::Bin* _bin{nullptr};
  ExtrusionBarrier* _blocking_barrier{nullptr};
  DNA::Direction _direction{};
  uint32_t _stalls_left{0};
  std::geometric_distribution<uint32_t> _n_of_stalls_gen;

  void set_stalls(uint32_t n);
  void increment_stalls(uint32_t n = 1);
  void decrement_stalls(uint32_t n = 1);
  void reset_stalls();
  void unload();
  void bind(Chromosome* chr, uint32_t pos, DNA::Direction direction);

  uint32_t check_for_extruder_collisions(std::mt19937& rand_gen);
  [[nodiscard]] uint64_t check_for_extrusion_barrier(std::mt19937& rand_gen);
  bool try_moving_to_next_bin();
  bool try_moving_to_prev_bin();
  [[nodiscard]] bool hard_stall() const;
};

}  // namespace modle
