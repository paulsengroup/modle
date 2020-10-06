#ifndef MODLE_CONTACTS_HPP
#define MODLE_CONTACTS_HPP

#include "absl/container/flat_hash_map.h"

namespace modle {

class ContactMatrix {
 public:
  ContactMatrix(uint32_t nrows, uint32_t ncols);
  [[nodiscard]] uint32_t get(uint32_t row, uint32_t col) const;
  void set(uint32_t row, uint32_t col, uint32_t n = 1);
  void increment(uint32_t row, uint32_t col, uint32_t n = 1);
  void decrement(uint32_t row, uint32_t col, uint32_t n = 1);
  void print() const;
  void print_symmetric_matrix() const;
  [[nodiscard]] std::vector<std::vector<uint32_t>> generate_symmetric_matrix() const;
  [[nodiscard]] uint32_t n_cols() const;
  [[nodiscard]] uint32_t n_rows() const;
  void write_full_matrix_to_tsv(const std::string& path_to_file) const;
  void write_to_tsv(const std::string& path_to_file) const;

 private:
  uint64_t _nrows;
  uint64_t _ncols;
  std::vector<std::vector<uint32_t>> _matrix;

  [[nodiscard]] std::vector<std::vector<uint32_t>> allocate_matrix() const;
};

}  // namespace modle
#endif  // MODLE_CONTACTS_HPP
