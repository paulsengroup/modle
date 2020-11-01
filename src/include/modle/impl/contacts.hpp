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
    : _nrows(nrows), _ncols(ncols), _matrix(_nrows * _ncols, 0) {}

template <typename I>
I &ContactMatrix<I>::at(uint64_t i, uint64_t j) {
  assert((j * this->_nrows) + i < this->_matrix.size());
  return this->_matrix[(j * this->_nrows) + i];
}

template <typename I>
const I &ContactMatrix<I>::at(uint64_t i, uint64_t j) const {
  assert((j * this->_nrows) + i < this->_matrix.size());
  return this->_matrix[(j * this->_nrows) + i];
}

template <typename I>
void ContactMatrix<I>::increment(uint64_t row, uint64_t col, I n) {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (j >= this->n_cols())
    throw std::logic_error(absl::StrFormat("ContactMatrix::increment(%d, %d, %d): j=%d > ncols=%d.",
                                           row, col, n, j, this->n_cols()));
  if (this->_updates_missed += i >= this->n_rows()) return;

  assert(i < this->n_rows());
  assert(j < this->n_cols());
  this->at(i, j) += n;
}

template <typename I>
void ContactMatrix<I>::set(uint64_t row, uint64_t col, I n) {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (j >= this->n_cols())
    throw std::logic_error(absl::StrFormat("ContactMatrix::set(%d, %d, %d): j=%d > ncols=%d.", row,
                                           col, n, j, this->n_cols()));
  if (this->_updates_missed += i >= this->n_rows()) return;

  assert(i < this->n_rows());
  assert(j < this->n_cols());
  this->at(i, j) = n;
}

template <typename I>
void ContactMatrix<I>::decrement(uint64_t row, uint64_t col, I n) {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (j >= this->n_cols())
    throw std::logic_error(absl::StrFormat("ContactMatrix::decrement(%d, %d, %d): j=%d > ncols=%d.",
                                           row, col, n, j, this->n_cols()));
  if (this->_updates_missed += i >= this->n_rows()) return;

  assert(i < this->n_rows());
  assert(j < this->n_cols());
  assert(n <= this->_matrix[(j * this->_nrows) + i]);
  this->at(i, j) -= n;
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
          row[x] = this->at(i, j);
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
I ContactMatrix<I>::get(uint64_t row, uint64_t col) const {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  if (i >= this->n_rows() || j >= this->n_cols()) return 0;
  return this->at(i, j);
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
        row[x] = this->at(i, j);
      }
    }
    m.emplace_back(std::move(row));
  }
  return m;
}

template <typename I>
uint64_t ContactMatrix<I>::n_rows() const {
  return this->_nrows;
}

template <typename I>
uint64_t ContactMatrix<I>::n_cols() const {
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
  std::vector<I> row(this->_ncols, 0);
  std::string buff;

  for (auto i = 0UL; i < this->_nrows; ++i) {
    for (auto j = 0UL; j < this->_ncols; ++j) {
      assert(i * j < this->_matrix.size());
      row[j] = this->at(i, j);
    }

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
        row[x] = this->at(i, j);
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

}  // namespace modle
