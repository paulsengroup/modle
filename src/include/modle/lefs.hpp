#ifndef MODLE_LEFS_HPP
#define MODLE_LEFS_HPP

#include <cstdint>
#include <memory>
#include <random>

#include "modle/contacts.hpp"
#include "modle/dna.hpp"
#include "modle/extr_barrier.hpp"

namespace modle {
class Lef;

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

class Lef {
 public:
  Lef(uint32_t bin_size, uint32_t avg_processivity, double probability_of_extruder_bypass);

  uint32_t extrude();
  void register_contact();
  [[nodiscard]] std::pair<DNA::Bin*, DNA::Bin*> get_ptr_to_bins();
  [[nodiscard]] bool is_bound() const;
  void bind_at_pos(Chromosome& chr, uint32_t pos,
                   std::mt19937& rand_gen);  // Returns false if bin is already occupied
  [[nodiscard]] std::string_view get_chr_name() const;
  void check_constraints(std::mt19937& rand_gen);
  [[nodiscard]] uint32_t get_loop_size() const;
  [[nodiscard]] uint32_t get_avg_processivity() const;
  bool try_rebind(Chromosome& chr, std::mt19937& rand_eng, double prob_of_rebinding);
  [[nodiscard]] const DNA::Bin& get_first_bin() const;
  [[nodiscard]] const DNA::Bin& get_last_bin() const;
  [[nodiscard]] Chromosome* get_ptr_to_chr();
  [[nodiscard]] std::pair<uint32_t, uint32_t> get_pos() const;
  [[nodiscard]] double get_probability_of_extr_unit_bypass() const;

 private:
  Chromosome* _chr{nullptr};
  uint32_t _lifetime{0};
  uint32_t _avg_processivity;
  double _probability_of_extr_unit_bypass;
  std::geometric_distribution<uint32_t> _lifetime_generator;
  std::unique_ptr<ExtrusionUnit> _left_unit;
  std::unique_ptr<ExtrusionUnit> _right_unit;

  void unload();

  [[nodiscard]] double compute_prob_of_unloading(uint32_t bin_size,
                                                 uint8_t n_of_active_extr_units = 2) const;
};

}  // namespace modle

#endif  // MODLE_LEFS_HPP
