#ifndef MODLE_LEFS_HPP
#define MODLE_LEFS_HPP

#include <cstdint>
#include <memory>
#include <random>

#include "modle/dna.hpp"

namespace modle {
class Lef {
 public:
  explicit Lef(uint32_t binding_pos, uint32_t avg_processivity);
  Lef(uint32_t binding_pos, uint32_t avg_processivity, uint8_t left_extr_speed,
      uint8_t right_extr_speed);

  bool extrude(std::default_random_engine &rng);
  [[nodiscard]] std::pair<std::shared_ptr<DNA::Bin>, std::shared_ptr<DNA::Bin>> get_bins();
  [[nodiscard]] bool is_bound() const;
  [[nodiscard]] bool bind_at_pos(std::string_view chr_name, DNA &chr, uint32_t pos); // Returns false if bin is already occupied
  [[nodiscard]] std::string_view get_chr_name() const;

 private:
  std::string_view _chr{};
  uint32_t _left_pos;
  uint32_t _right_pos;
  std::shared_ptr<DNA::Bin> _left_bin{nullptr};
  std::shared_ptr<DNA::Bin> _right_bin{nullptr};
  uint8_t _left_extrusion_speed{1};
  uint8_t _right_extrusion_speed{1};
  uint64_t _avg_processivity;
  std::bernoulli_distribution _dist;
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
