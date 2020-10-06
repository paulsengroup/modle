#include "modle/contacts.hpp"

#include <zlib.h>

#include <cassert>
#include <cmath>

#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"

namespace modle {

ContactMatrix::ContactMatrix(uint32_t nrows, uint32_t ncols)
    : _nrows(nrows), _ncols(ncols), _matrix(allocate_matrix()) {}

void ContactMatrix::increment(uint32_t row, uint32_t col, uint32_t n) {
  uint32_t j = col;
  uint32_t i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (i >= this->_matrix.size() || j >= this->_matrix[0].size()) {
    absl::FPrintF(
        stderr,
        "WARNING: ContactMatrix::increment: ignoring one increment for matrix[%lu][%lu] += %lu! "
        "Reason: out of range (nrows=%lu; ncols=%lu;)\n",
        i, j, n, this->_nrows, this->_ncols);
    return;  // Temporary workaround to avoid segfaults
  }
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
  if (i >= this->_matrix.size() || j >= this->_matrix[0].size()) {
    absl::FPrintF(stderr,
                  "WARNING: ContactMatrix::set: ignoring one update for matrix[%lu][%lu] = %lu! "
                  "Reason: out of range (nrows=%lu; ncols=%lu;)\n",
                  i, j, n, this->_nrows, this->_ncols);
    return;  // Temporary workaround to avoid segfaults
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
  if (i >= this->_matrix.size() || j >= this->_matrix[0].size()) {
    absl::FPrintF(
        stderr,
        "WARNING: ContactMatrix::increment: ignoring one decrement for matrix[%lu][%lu] -= %lu! "
        "Reason: out of range (nrows=%lu; ncols=%lu;)\n",
        i, j, n, this->_nrows, this->_ncols);
    return;  // Temporary workaround to avoid segfaults
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
  m.reserve(this->_nrows);

  for (uint32_t i = 0; i < this->_nrows; ++i) {
    std::vector<uint32_t> v(this->_ncols, 0);
    m.emplace_back(std::move(v));
  }
  return m;
}

uint32_t ContactMatrix::get(uint32_t row, uint32_t col) const {
  uint32_t j = col;
  uint32_t i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (i >= this->_matrix.size() || j >= this->_matrix[0].size()) return 0;
  return this->_matrix[i][j];
}

void ContactMatrix::print_symmetric_matrix() const {
  std::string buff;
  buff.reserve(64 * 1024 * 1024);
  for (uint32_t y = 0; y < this->_ncols; ++y) {
    //    absl::FPrintF(stderr, "%u\t", y);
    for (uint32_t x = 0; x < this->_ncols; ++x) {
      uint32_t j = x;
      uint32_t i = j - y;
      if (y > x) {
        j = y;
        i = j - x;
      }

      if (i >= this->_nrows) {
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

uint32_t ContactMatrix::n_rows() const {
  assert(this->_nrows == this->_matrix.size());
  return this->_matrix.size();
}
uint32_t ContactMatrix::n_cols() const {
  assert(this->_ncols == this->_matrix.at(0).size());
  return this->_matrix[0].size();
}

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
  std::vector<uint32_t> row(this->_ncols, 0);
  std::string buff;
  for (uint32_t y = 0; y < this->_ncols; ++y) {
    //    absl::FPrintF(stderr, "%u\t", y);
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
    buff = absl::StrFormat("%s\n", absl::StrJoin(row, "\t"));
    gzwrite(gzf, buff.c_str(), buff.size());
  }
  if (gzclose(gzf) != Z_OK) {
    throw std::runtime_error(
        absl::StrFormat("An error occurred while closing the following file '%s'.", path_to_file));
  }
}

}  // namespace modle