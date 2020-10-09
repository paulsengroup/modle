#ifndef MODLE_CONTACTS_HPP
#define MODLE_CONTACTS_HPP

#include "absl/container/flat_hash_map.h"

namespace modle {

template<typename I>
class ContactMatrix {
 public:
  ContactMatrix(uint64_t nrows, uint64_t ncols);
  [[nodiscard]] I get(uint64_t row, uint64_t col) const;
  void set(uint64_t row, uint64_t col, I n = 1);
  void increment(uint64_t row, uint64_t col, I n = 1);
  void decrement(uint64_t row, uint64_t col, I n = 1);
  void print() const;
  [[nodiscard]] std::vector<std::vector<I>> generate_symmetric_matrix() const;
  [[nodiscard]] uint64_t n_cols() const;
  [[nodiscard]] uint64_t n_rows() const;
  void write_full_matrix_to_tsv(const std::string& path_to_file) const;
  void write_to_tsv(const std::string& path_to_file) const;

 private:
  uint64_t _nrows;
  uint64_t _ncols;
  std::vector<std::vector<I>> _matrix;
  uint64_t _updates_missed{0};

  [[nodiscard]] std::vector<std::vector<uint32_t>> allocate_matrix() const;
};

}  // namespace modle

#include "modle/impl/contacts.hpp"

#endif  // MODLE_CONTACTS_HPP
