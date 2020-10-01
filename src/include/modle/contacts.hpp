#ifndef MODLE_CONTACTS_HPP
#define MODLE_CONTACTS_HPP

#include "absl/container/flat_hash_map.h"

namespace modle {

class ContactMatrix {
 public:
  ContactMatrix(uint32_t width, uint32_t length);
  ContactMatrix(uint32_t width, uint32_t x_len, uint32_t y_len);
  [[nodiscard]] uint32_t get(uint32_t x, uint32_t y) const;
  void increment(uint32_t x, uint32_t y, uint32_t n = 1);
  void decrement(uint32_t i, uint32_t j, uint32_t n = 1);
  void print() const;
  void print_symmetric_matrix() const;
  [[nodiscard]] std::vector<std::vector<uint32_t>> generate_symmetric_matrix() const;

 private:
  uint64_t _width;
  uint64_t _length;
  uint64_t _apparent_x;
  uint64_t _apparent_y;
  std::vector<std::vector<uint32_t>> _matrix;

  [[nodiscard]] std::vector<std::vector<uint32_t>> allocate_matrix() const;
  [[nodiscard]] static uint32_t compute_diagonal_length(uint64_t x, uint64_t y);
};

}  // namespace modle
#endif  // MODLE_CONTACTS_HPP
