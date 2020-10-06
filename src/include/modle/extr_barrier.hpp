#ifndef MODLE_EXTR_BARRIER_HPP
#define MODLE_EXTR_BARRIER_HPP

#include <cstdint>
#include <random>

namespace modle {
class ExtrusionBarrier {
 public:
  ExtrusionBarrier(uint32_t abs_position, uint32_t bin_size, double prob_of_block, uint64_t seed);
  [[nodiscard]] double get_prob_of_blocking() const;
  [[nodiscard]] bool is_blocking();
  void relocate(uint32_t new_abs_pos, uint32_t bin_size);
  [[nodiscard]] uint32_t get_abs_pos() const;
  [[nodiscard]] uint32_t get_rel_pos() const;

 private:
  uint32_t _abs_position;
  uint32_t _rel_position;  // Possibly useless
  std::default_random_engine _rng;
  std::bernoulli_distribution _block;
  uint64_t _seed;

  static uint32_t abs_to_rel_position(uint32_t abs_pos, uint32_t bin_size);
};
}  // namespace modle

#endif  // MODLE_EXTR_BARRIER_HPP
