#ifndef MODLE_CONTACTS_HPP
#define MODLE_CONTACTS_HPP

#include "absl/container/flat_hash_map.h"

namespace modle {

class ContactMatrix {
 public:
  ContactMatrix(uint32_t width, uint32_t length);
  [[nodiscard]] uint32_t get(uint32_t row, uint32_t col) const;
  void set(uint32_t row, uint32_t col, uint32_t n = 1);
  void increment(uint32_t row, uint32_t col, uint32_t n = 1);
  void decrement(uint32_t row, uint32_t col, uint32_t n = 1);
  void print() const;
  void print_symmetric_matrix() const;
  [[nodiscard]] std::vector<std::vector<uint32_t>> generate_symmetric_matrix() const;
  [[nodiscard]] uint32_t n_cols() const;
  [[nodiscard]] uint32_t n_rows() const;
  [[nodiscard]] uint32_t apparent_n_cols() const;
  [[nodiscard]] uint32_t apparent_n_rows() const;
  void write_full_matrix_to_tsv(const std::string& path_to_file) const;
  void write_to_tsv(const std::string& path_to_file) const;
  // Write method to output table directly, without conversions to "original" coordinates

 private:
  uint64_t _width;
  uint64_t _length;
  std::vector<std::vector<uint32_t>> _matrix;

  [[nodiscard]] std::vector<std::vector<uint32_t>> allocate_matrix() const;
  [[nodiscard]] static uint32_t compute_diagonal_length(uint64_t x, uint64_t y);
};

}  // namespace modle
#endif  // MODLE_CONTACTS_HPP
