#ifndef MODLE_LEFS_HPP
#define MODLE_LEFS_HPP

#include <cstdint>
#include <memory>
#include <random>

#include "modle/dna.hpp"

namespace modle {
class Lef {
 public:
  explicit Lef() = default;
  explicit Lef(uint32_t avg_processivity);
  Lef(uint32_t binding_pos, uint32_t avg_processivity, DNA* = nullptr, std::string_view chr = "");
  Lef(uint32_t binding_pos, uint32_t avg_processivity, uint8_t left_extr_speed,
      uint8_t right_extr_speed, DNA* = nullptr, std::string_view chr = "");

  uint32_t extrude();
  [[nodiscard]] std::pair<std::shared_ptr<DNA::Bin>, std::shared_ptr<DNA::Bin>> get_bins();
  [[nodiscard]] bool is_bound() const;
  [[nodiscard]] bool bind_at_pos(std::string_view chr_name, DNA& dna,
                                 uint32_t pos);  // Returns false if bin is already occupied
  [[nodiscard]] std::string_view get_chr_name() const;
  void check_constraints();
  [[nodiscard]] uint32_t get_left_extrusion_speed() const;
  [[nodiscard]] uint32_t get_right_extrusion_speed() const;
  [[nodiscard]] bool left_is_stalled() const;
  [[nodiscard]] bool right_is_stalled() const;
  [[nodiscard]] uint32_t get_loop_size() const;
  [[nodiscard]] uint64_t get_avg_processivity() const;
  void set_avg_processivity(uint32_t avg_proc);
  bool try_unload(std::default_random_engine& rng);
  bool try_rebind(std::string_view chr_name, DNA& chr, std::default_random_engine& rng,
                  double prob_of_rebinding = 1.0);

 private:
  std::string_view _chr{};
  DNA* _dna{nullptr};
  uint32_t _left_pos{UINT32_MAX};
  uint32_t _right_pos{UINT32_MAX};
  std::shared_ptr<DNA::Bin> _left_bin{nullptr};
  std::shared_ptr<DNA::Bin> _right_bin{nullptr};
  uint8_t _left_extrusion_speed{1};
  uint8_t _right_extrusion_speed{1};
  bool _stall_left{false};
  bool _stall_right{false};
  uint32_t _avg_processivity{0};
  std::bernoulli_distribution _bernoulli_dist;

  void unload();
  void stall_left();
  void stall_right();

  void remove_left_stall();
  void remove_right_stall();
  double compute_prob_of_unloading() const;
};

class Lsf {
 public:
  explicit Lsf(uint32_t left_pos, uint32_t right_pos, uint32_t lifetime);
  bool next();

 private:
  uint32_t _left_pos;
  uint32_t _right_pos;
  uint32_t _lifetime;
};

}  // namespace modle

#endif  // MODLE_LEFS_HPP
