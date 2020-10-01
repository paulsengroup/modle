#include "modle/contacts.hpp"

#include <cassert>
#include <cmath>

#include "absl/strings/str_format.h"

namespace modle {

ContactMatrix::ContactMatrix(uint32_t width, uint32_t length)
    : _width(width),
      _length(length),
      _apparent_x(std::sqrt((length * length) / 2)),
      _apparent_y(std::sqrt((length * length) / 2)),
      _matrix(allocate_matrix()) {}

ContactMatrix::ContactMatrix(uint32_t width, uint32_t x_len, uint32_t y_len)
    : _width(width),
      _length(compute_diagonal_length(x_len, y_len)),
      _apparent_x(x_len),
      _apparent_y(y_len),
      _matrix(allocate_matrix()) {}

void ContactMatrix::increment(uint32_t x, uint32_t y, uint32_t n) {
  // TODO: Given that we are only storing half of the symmetric matrix, it may be a good idea to
  // invert the values of i and j when y > x (see print_symmetric for an example)
  uint32_t j = y;
  uint32_t i = j - x;
  assert(i < this->_matrix.size());
  assert(j < this->_matrix[0].size());
  this->_matrix[i][j] += n;
}

void ContactMatrix::decrement(uint32_t x, uint32_t y, uint32_t n) {
  // TODO: Given that we are only storing half of the symmetric matrix, it may be a good idea to
  // invert the values of i and j when y > x (see print_symmetric for an example)
  uint32_t j = y;
  uint32_t i = j - x;
  assert(i < this->_matrix.size());
  assert(j < this->_matrix[0].size());
  assert(n < this->_matrix[i][j]);
  this->_matrix[i][j] -= n;
}

void ContactMatrix::print() const {
  for (const auto &row : this->_matrix) {
    for (const auto &n : row) {
      absl::PrintF("%u\t", n);
    }
    absl::PrintF("\n");
  }
}

std::vector<std::vector<uint32_t>> ContactMatrix::allocate_matrix() const {
  std::vector<std::vector<uint32_t>> m;
  m.reserve(this->_width);

  for (uint32_t i = 0; i < this->_width; ++i) {
    std::vector<uint32_t> v(this->_length + this->_width, 0);
    m.emplace_back(std::move(v));
  }
  assert(this->_matrix.size() == this->_width);
  for (const auto &row : this->_matrix) assert(row.size() == this->_length + this->_width);

  return m;
}

uint32_t ContactMatrix::get(uint32_t x, uint32_t y) const {
  // TODO: Given that we are only storing half of the symmetric matrix, it may be a good idea to
  // invert the values of i and j when y > x (see print_symmetric for an example)
  uint32_t j = y;
  uint32_t i = j - x;
  assert(i < this->_matrix.size());
  assert(j < this->_matrix[0].size());
  return this->_matrix[i][j];
}

uint32_t ContactMatrix::compute_diagonal_length(uint64_t x, uint64_t y) {
  return std::sqrt((x * x) + (y * y));
}

void ContactMatrix::print_symmetric_matrix() const {
  absl::PrintF(" \t");
  for (uint32_t y = 0; y < this->_apparent_y; ++y) {
    absl::PrintF("%u\t", y);
  }
  absl::PrintF("\n");
  for (uint32_t y = 0; y < this->_apparent_y; ++y) {
    absl::PrintF("%u\t", y);
    for (uint32_t x = 0; x < this->_apparent_x; ++x) {
      uint32_t j = x;
      uint32_t i = j - y;
      if (y > x) {
        j = y;
        i = j - x;
      }

      if (i >= this->_width) {
        absl::PrintF("0\t");
      } else {
        absl::PrintF("%d\t", this->_matrix[i][j]);
      }
    }
    absl::PrintF("\n");
  }
}

std::vector<std::vector<uint32_t>> ContactMatrix::generate_symmetric_matrix() const {
  std::vector<std::vector<uint32_t>> m;
  m.reserve(this->_apparent_y);
  for (uint32_t y = 0; y < this->_apparent_y; ++y) {
    std::vector<uint32_t> row(this->_apparent_x, 0);
    for (uint32_t x = 0; x < this->_apparent_x; ++x) {
      uint32_t j = x;
      uint32_t i = j - y;
      if (y > x) {
        j = y;
        i = j - x;
      }

      if (i < this->_width) {
        row[i] = this->_matrix[i][j];
      }
    }
    m.emplace_back(std::move(row));
  }
  return m;
}

}  // namespace modle