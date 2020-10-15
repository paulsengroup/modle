// #include <zlib.h>
#include <bzlib.h>

#include <cassert>
#include <cmath>

#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "modle/contacts.hpp"

namespace modle {
template <typename I>
ContactMatrix<I>::ContactMatrix(uint64_t nrows, uint64_t ncols)
    : _nrows(nrows), _ncols(ncols), _matrix(allocate_matrix()) {}

template <typename I>
void ContactMatrix<I>::increment(uint64_t row, uint64_t col, I n) {
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

template <typename I>
void ContactMatrix<I>::set(uint64_t row, uint64_t col, I n) {
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

template <typename I>
void ContactMatrix<I>::decrement(uint64_t row, uint64_t col, I n) {
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

template <typename I>
void ContactMatrix<I>::print(bool full) const {
  if (full) {
    std::vector<I> row(this->_ncols, 0);
    for (uint64_t y = 0; y < this->_ncols; ++y) {
      std::fill(row.begin(), row.end(), 0);
      for (uint64_t x = 0; x < this->_ncols; ++x) {
        auto j = x;
        auto i = j - y;
        if (y > x) {
          j = y;
          i = j - x;
        }
        if (i >= this->_nrows)
          row[x] = 0;
        else
          row[x] = this->_matrix[i][j];
      }
      absl::PrintF("%s\n", absl::StrJoin(row, "\t"));
    }
  } else {
    for (const auto &row : this->_matrix) {
      absl::PrintF("%s\n", absl::StrJoin(row, "\t"));
    }
  }
}

template <typename I>
std::vector<std::vector<uint32_t>> ContactMatrix<I>::allocate_matrix() const {
  std::vector<std::vector<uint32_t>> m;
  m.reserve(this->_nrows);

  for (uint64_t i = 0; i < this->_nrows; ++i) {
    std::vector<uint32_t> v(this->_ncols, 0);
    m.emplace_back(std::move(v));
  }
  return m;
}

template <typename I>
I ContactMatrix<I>::get(uint64_t row, uint64_t col) const {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (i >= this->n_rows() || j >= this->n_cols()) return 0;
  return this->_matrix[i][j];
}

template <typename I>
std::vector<std::vector<I>> ContactMatrix<I>::generate_symmetric_matrix() const {
  std::vector<std::vector<I>> m;
  m.reserve(this->_ncols);
  for (uint64_t y = 0; y < this->_ncols; ++y) {
    std::vector<I> row(this->_ncols, 0);
    for (I x = 0; x < this->_ncols; ++x) {
      auto j = x;
      auto i = j - y;
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

template <typename I>
uint64_t ContactMatrix<I>::n_rows() const {
  assert(this->_nrows == this->_matrix.size());
  return this->_nrows;
}

template <typename I>
uint64_t ContactMatrix<I>::n_cols() const {
  assert(this->_ncols == this->_matrix.at(0).size());
  return this->_ncols;
}

template <typename I>
std::pair<uint32_t, uint32_t> ContactMatrix<I>::write_to_tsv(
    const std::string &path_to_file) const {
  auto fp = fopen(path_to_file.c_str(), "wb");
  int bz_status;
  if (!fp)
    throw std::runtime_error(
        absl::StrFormat("An error occurred while opening file '%s' for writing.", path_to_file));
  if (this->_updates_missed > 0) {
    absl::FPrintF(stderr, "WARNING: There were %lu missed updates!\n", this->_updates_missed);
  }

  auto bzf = BZ2_bzWriteOpen(&bz_status, fp, 9 /* block size */, 0 /* verbosity */,
                             /* work factor, 0 == default == 30 */ 0);
  if (bz_status != BZ_OK) {
    BZ2_bzWriteClose(nullptr, bzf, 0, nullptr, nullptr);
    fclose(fp);
    throw std::runtime_error(
        absl::StrFormat("Unable to open file '%s' for writing.", path_to_file));
  }
  std::string buff;
  for (const auto &row : this->_matrix) {
    buff = absl::StrJoin(row, "\t");
    buff += "\n";
    BZ2_bzWrite(&bz_status, bzf, buff.data(), buff.size());
    if (bz_status != BZ_OK) {
      BZ2_bzWriteClose(nullptr, bzf, 0, nullptr, nullptr);
      fclose(fp);
      throw std::runtime_error(
          absl::StrFormat("An error occurred while writing to file '%s'.", path_to_file));
    }
  }
  uint32_t bytes_in, bytes_out;
  BZ2_bzWriteClose(&bz_status, bzf, 0, &bytes_in, &bytes_out);
  if (fclose(fp) != 0 || bz_status != BZ_OK) {
    throw std::runtime_error(
        absl::StrFormat("An error occurred while closing file '%s'.", path_to_file));
  }
  return std::make_pair(bytes_in, bytes_out);
}

template <typename I>
std::pair<uint32_t, uint32_t> ContactMatrix<I>::write_full_matrix_to_tsv(
    const std::string &path_to_file) const {
  auto fp = fopen(path_to_file.c_str(), "wb");
  int bz_status;
  if (!fp)
    throw std::runtime_error(
        absl::StrFormat("An error occurred while opening file '%s' for writing.", path_to_file));
  if (this->_updates_missed > 0) {
    absl::FPrintF(stderr, "WARNING: There were %lu missed updates!\n", this->_updates_missed);
  }

  auto bzf = BZ2_bzWriteOpen(&bz_status, fp, 9 /* block size */, 0 /* verbosity */,
                             /* work factor, 0 == default == 30 */ 0);
  if (bz_status != BZ_OK) {
    BZ2_bzWriteClose(nullptr, bzf, 0, nullptr, nullptr);
    fclose(fp);
    throw std::runtime_error(
        absl::StrFormat("Unable to open file '%s' for writing.", path_to_file));
  }
  std::vector<I> row(this->_ncols, 0);
  std::string buff;
  for (uint64_t y = 0; y < this->_ncols; ++y) {
    //    std::fill(row.begin(), row.end(), 0);

    for (uint64_t x = 0; x < this->_ncols; ++x) {
      auto j = x;
      auto i = j - y;
      if (y > x) {
        j = y;
        i = j - x;
      }
      if (i >= this->_nrows)
        row[x] = 0;
      else
        row[x] = this->_matrix[i][j];
    }
    buff = absl::StrJoin(row, "\t") + "\n";
    BZ2_bzWrite(&bz_status, bzf, buff.data(), buff.size());
    if (bz_status != BZ_OK) {
      BZ2_bzWriteClose(nullptr, bzf, 0, nullptr, nullptr);
      fclose(fp);
      throw std::runtime_error(
          absl::StrFormat("An error occurred while writing to file '%s'.", path_to_file));
    }
  }
  uint32_t bytes_in, bytes_out;
  BZ2_bzWriteClose(&bz_status, bzf, 0, &bytes_in, &bytes_out);
  if (fclose(fp) != 0 || bz_status != BZ_OK) {
    throw std::runtime_error(
        absl::StrFormat("An error occurred while closing file '%s'.", path_to_file));
  }
  return std::make_pair(bytes_in, bytes_out);
}

template <typename I>
void ContactMatrix<I>::handle_bzip_write_errors(int status, std::string_view path) {
  if (status != BZ_OK)
    throw std::runtime_error(
        absl::StrFormat("An error occurred while compressing file '%s'.", path));
}

}  // namespace modle
