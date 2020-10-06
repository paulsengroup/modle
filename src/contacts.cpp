#include "modle/contacts.hpp"

#include <zlib.h>

#include <cassert>
#include <cmath>

#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"

namespace modle {

ContactMatrix::ContactMatrix(uint32_t width, uint32_t length)
    : _width(width), _length(length), _matrix(allocate_matrix()) {}

void ContactMatrix::increment(uint32_t row, uint32_t col, uint32_t n) {
  uint32_t j = col;
  uint32_t i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (i >= this->_matrix.size() || j >= this->_matrix[0].size()) return; // Temporary workaround to avoid segfaults
  assert(i < this->_matrix.size());
  assert(j < this->_matrix[0].size());
  this->_matrix[i][j] += n;
}

void ContactMatrix::set(uint32_t row, uint32_t col, uint32_t n) {
  uint32_t j = col;
  uint32_t i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  assert(i < this->_matrix.size());
  assert(j < this->_matrix[0].size());
  this->_matrix[i][j] = n;
}

void ContactMatrix::decrement(uint32_t row, uint32_t col, uint32_t n) {
  uint32_t j = col;
  uint32_t i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  assert(i < this->_matrix.size());
  assert(j < this->_matrix[0].size());
  assert(n <= this->_matrix[i][j]);
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

uint32_t ContactMatrix::get(uint32_t row, uint32_t col) const {
  uint32_t j = col;
  uint32_t i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  assert(i < this->_matrix.size());
  assert(j < this->_matrix[0].size());
  return this->_matrix[i][j];
}

uint32_t ContactMatrix::compute_diagonal_length(uint64_t x, uint64_t y) {
  return std::sqrt((x * x) + (y * y));
}

void ContactMatrix::print_symmetric_matrix() const {
  std::string buff;
  buff.reserve(64 * 1024 * 1024);
  for (uint32_t y = 0; y < this->_length; ++y) {
    //    absl::FPrintF(stderr, "%u\t", y);
    for (uint32_t x = 0; x < this->_length; ++x) {
      uint32_t j = x;
      uint32_t i = j - y;
      if (y > x) {
        j = y;
        i = j - x;
      }

      if (i >= this->_width) {
        buff += "0\t";
      } else {
        buff += std::to_string(this->_matrix[i][j]) + "\t";
      }
    }
    buff += "\n";
    if (buff.size() >= 64 * 1024 * 1023) {
      absl::FPrintF(stderr, "%s", buff);
      buff.clear();
    }
  }
  if (!buff.empty()) {
    absl::FPrintF(stderr, "%s", buff);
  }
}

std::vector<std::vector<uint32_t>> ContactMatrix::generate_symmetric_matrix() const {
  std::vector<std::vector<uint32_t>> m;
  m.reserve(this->_length);
  for (uint32_t y = 0; y < this->_length; ++y) {
    std::vector<uint32_t> row(this->_length, 0);
    for (uint32_t x = 0; x < this->_length; ++x) {
      uint32_t j = x;
      uint32_t i = j - y;
      if (y > x) {
        j = y;
        i = j - x;
      }

      if (i < this->_width) {
        row[x] = this->_matrix[i][j];
      }
    }
    m.emplace_back(std::move(row));
  }
  return m;
}

uint32_t ContactMatrix::n_rows() const { return this->_matrix.size(); }
uint32_t ContactMatrix::n_cols() const { return this->_matrix[0].size(); }
uint32_t ContactMatrix::apparent_n_rows() const { return this->_length; }
uint32_t ContactMatrix::apparent_n_cols() const { return this->_length; }

void ContactMatrix::write_to_tsv(const std::string &path_to_file) const {
  auto gzf = gzopen(path_to_file.c_str(), "w");
  if (gzf == Z_NULL) {
    throw std::runtime_error(
        absl::StrFormat("Unable to open file '%s' for writing.", path_to_file));
  }
  std::string buff;
  buff.reserve(1024 * 1024);
  for (const auto &row : this->_matrix) {
    buff += absl::StrJoin(row, "\t");
    buff += "\n";
    gzwrite(gzf, buff.c_str(), buff.size());
    buff.clear();
  }
  if (gzclose(gzf) != Z_OK) {
    throw std::runtime_error(
        absl::StrFormat("An error occurred while closing the following file '%s'.", path_to_file));
  }
}

void ContactMatrix::write_full_matrix_to_tsv(const std::string &path_to_file) const {
  auto gzf = gzopen(path_to_file.c_str(), "w");
  if (gzf == Z_NULL) {
    throw std::runtime_error(
        absl::StrFormat("Unable to open file '%s' for writing.", path_to_file));
  }
  std::string buff;
  buff.reserve(1024 * 1024);
  for (uint32_t y = 0; y < this->_length; ++y) {
    //    absl::FPrintF(stderr, "%u\t", y);
    for (uint32_t x = 0; x < this->_length; ++x) {
      uint32_t j = x;
      uint32_t i = j - y;
      if (y > x) {
        j = y;
        i = j - x;
      }

      if (i >= this->_width) {
        buff += "0\t";
      } else {
        buff += absl::StrFormat("%lu\t", this->_matrix[i][j]);
      }
    }
    buff += "\n";
    if (buff.size() >= 1024 * 1023) {
      gzwrite(gzf, buff.c_str(), buff.size());
      buff.clear();
    }
  }
  if (!buff.empty()) {
    gzwrite(gzf, buff.c_str(), buff.size());
  }
  if (gzclose(gzf) != Z_OK) {
    throw std::runtime_error(
        absl::StrFormat("An error occurred while closing the following file '%s'.", path_to_file));
  }
}

}  // namespace modle