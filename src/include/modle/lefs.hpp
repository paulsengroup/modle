#ifndef MODLE_LEFS_HPP
#define MODLE_LEFS_HPP

#include <cstdint>
#include <memory>
#include <random>

#include "modle/dna.hpp"

namespace modle {

class ExtrusionUnit {
  friend class Lef;

 public:
  enum Direction { fwd = true, rev = false };
  explicit ExtrusionUnit() = default;
  ExtrusionUnit(uint32_t pos, DNA* dna, std::shared_ptr<DNA::Bin> bin, uint32_t extrusion_speed,
                Direction direction, bool stall = false);
  [[nodiscard]] uint32_t get_pos() const;

 private:
  uint32_t _pos{UINT32_MAX};
  DNA* _dna{nullptr};
  std::shared_ptr<DNA::Bin> _bin{nullptr};
  uint32_t _extrusion_speed{1};
  Direction _direction{};
  bool _stalled{false};

  uint32_t extrude();
  void set_stall();
  void remove_stall();
  void unload();
  bool bind(DNA* dna, uint32_t pos, std::shared_ptr<DNA::Bin> bin, Direction direction);
  [[nodiscard]] bool is_bound() const;
  [[nodiscard]] bool is_extruding_fwd() const;
  [[nodiscard]] bool is_extruding_rev() const;
  [[nodiscard]] uint32_t get_extrusion_speed() const;
};

class Lef {
 public:
  explicit Lef() = default;
  explicit Lef(uint32_t avg_processivity);
  //  Lef(uint32_t binding_pos, uint32_t avg_processivity, DNA* = nullptr, std::string_view chr =
  //  "");
  //  Lef(uint32_t binding_pos, uint32_t avg_processivity, uint8_t left_extr_speed,
  //      uint8_t right_extr_speed, DNA* = nullptr, std::string_view chr = "");

  uint32_t extrude();
  [[nodiscard]] std::pair<DNA::Bin*, DNA::Bin*> get_ptr_to_bins();
  [[nodiscard]] bool is_bound() const;
  [[nodiscard]] bool bind_at_pos(std::string_view chr_name, DNA& dna,
                                 uint32_t pos);  // Returns false if bin is already occupied
  [[nodiscard]] std::string_view get_chr_name() const;
  void check_constrains();
  [[nodiscard]] uint32_t get_left_extrusion_speed() const;
  [[nodiscard]] uint32_t get_right_extrusion_speed() const;
  void stall_left();
  void stall_right();
  void remove_left_stall();
  void remove_right_stall();
  [[nodiscard]] bool left_is_stalled() const;
  [[nodiscard]] bool right_is_stalled() const;
  [[nodiscard]] uint32_t get_loop_size() const;
  [[nodiscard]] uint32_t get_avg_processivity() const;
  void set_avg_processivity(uint32_t avg_proc);
  bool try_unload(std::default_random_engine& rng);
  bool try_rebind(std::string_view chr_name, DNA& chr, std::default_random_engine& rng,
                  double prob_of_rebinding = 1.0);

 private:
  std::string_view _chr{};
  DNA* _dna{nullptr};

  ExtrusionUnit _left_unit{};
  ExtrusionUnit _right_unit{};
  uint32_t _avg_processivity{0};
  std::bernoulli_distribution _bernoulli_dist;

  void unload();

  [[nodiscard]] double compute_prob_of_unloading() const;
};

}  // namespace modle

#endif  // MODLE_LEFS_HPP
