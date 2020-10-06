#include "modle/contacts.hpp"

#include <zlib.h>

#include <cassert>
#include <cmath>

#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"

namespace modle {

ContactMatrix::ContactMatrix(uint64_t nrows, uint64_t ncols)
    : _nrows(nrows), _ncols(ncols), _matrix(allocate_matrix()) {}

void ContactMatrix::increment(uint64_t row, uint64_t col, uint32_t n) {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (i >= this->n_rows() || j >= this->n_cols()) {
    ++this->_updates_missed;
    return;
  }
  assert(i < this->n_rows());
  assert(j < this->n_cols());
  this->_matrix[i][j] += n;
}

void ContactMatrix::set(uint64_t row, uint64_t col, uint32_t n) {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (i >= this->n_rows() || j >= this->n_cols()) {
    ++this->_updates_missed;
    return;
  }
  assert(i < this->n_rows());
  assert(j < this->n_cols());
  this->_matrix[i][j] = n;
}

void ContactMatrix::decrement(uint64_t row, uint64_t col, uint32_t n) {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (i >= this->n_rows() || j >= this->n_cols()) {
    ++this->_updates_missed;
    return;
  }
  assert(i < this->n_rows());
  assert(j < this->n_cols());
  assert(n <= this->_matrix[i][j]);
  this->_matrix[i][j] -= n;
}

void ContactMatrix::print() const {
  for (const auto &row : this->_matrix) {
    absl::PrintF("%s\n", absl::StrJoin(row, "\t"));
  }
}

std::vector<std::vector<uint32_t>> ContactMatrix::allocate_matrix() const {
  std::vector<std::vector<uint32_t>> m;
  m.reserve(this->_nrows);

  for (uint64_t i = 0; i < this->_nrows; ++i) {
    std::vector<uint32_t> v(this->_ncols, 0);
    m.emplace_back(std::move(v));
  }
  return m;
}

uint32_t ContactMatrix::get(uint64_t row, uint64_t col) const {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (i >= this->n_rows() || j >= this->n_cols()) return 0;
  return this->_matrix[i][j];
}

std::vector<std::vector<uint32_t>> ContactMatrix::generate_symmetric_matrix() const {
  std::vector<std::vector<uint32_t>> m;
  m.reserve(this->_ncols);
  for (uint32_t y = 0; y < this->_ncols; ++y) {
    std::vector<uint32_t> row(this->_ncols, 0);
    for (uint32_t x = 0; x < this->_ncols; ++x) {
      uint32_t j = x;
      uint32_t i = j - y;
      if (y > x) {
        j = y;
        i = j - x;
      }

      if (i < this->_nrows) {
        row[x] = this->_matrix[i][j];
      }
    }
    m.emplace_back(std::move(row));
  }
  return m;
}

uint64_t ContactMatrix::n_rows() const {
  assert(this->_nrows == this->_matrix.size());
  return this->_nrows;
}
uint64_t ContactMatrix::n_cols() const {
  assert(this->_ncols == this->_matrix.at(0).size());
  return this->_ncols;
}

void ContactMatrix::write_to_tsv(const std::string &path_to_file) const {
  if (this->_updates_missed > 0) {
    absl::FPrintF(stderr, "WARNING: There were %lu missed updates!\n", this->_updates_missed);
  }
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
  if (this->_updates_missed > 0) {
    absl::FPrintF(stderr, "WARNING: There were %lu missed updates!\n", this->_updates_missed);
  }
  auto gzf = gzopen(path_to_file.c_str(), "w");
  if (gzf == Z_NULL) {
    throw std::runtime_error(
        absl::StrFormat("Unable to open file '%s' for writing.", path_to_file));
  }
  std::vector<uint32_t> row(this->_ncols, 0);
  std::string buff;
  for (uint32_t y = 0; y < this->_ncols; ++y) {
    for (uint32_t x = 0; x < this->_ncols; ++x) {
      uint32_t j = x;
      uint32_t i = j - y;
      if (y > x) {
        j = y;
        i = j - x;
      }
      if (i < this->_nrows)
        row[j] = this->_matrix[i][j];
      else
        row[j] = 0;
    }
    buff = absl::StrJoin(row, "\t") + "\n";
    gzwrite(gzf, buff.c_str(), buff.size());
  }
  if (gzclose(gzf) != Z_OK) {
    throw std::runtime_error(
        absl::StrFormat("An error occurred while closing the following file '%s'.", path_to_file));
  }
}

}  // namespace modle