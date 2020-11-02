#include <bzlib.h>

#include <array>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <cassert>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string_view>

#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "modle/contacts.hpp"
#include "modle/utils.hpp"

namespace modle {
template <class I>
ContactMatrix<I>::ContactMatrix(uint64_t nrows, uint64_t ncols, bool fill_with_random_numbers)
    : _nrows(nrows), _ncols(ncols), _matrix(_nrows * _ncols, 0) {
  if (fill_with_random_numbers) {
    std::random_device rand_dev;
    std::mt19937_64 rand_eng{rand_dev()};
    std::uniform_int_distribution<I> dist;
    std::generate(this->_matrix.begin(), this->_matrix.end(), [&]() { return dist(rand_eng); });
  }
}

template <class I>
ContactMatrix<I>::ContactMatrix(const std::string &path_to_file, char sep) : _nrows(0), _ncols(0) {
  std::ifstream fp(path_to_file, std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_istream in;
  std::string line;
  std::array<char, 8'192> buff{};
  fp.rdbuf()->pubsetbuf(buff.data(), buff.size());
  try {
    in.push(boost::iostreams::bzip2_decompressor());
    in.push(fp);
    for (this->_nrows = 0UL; std::getline(in, line); ++this->_nrows) {
      std::vector<std::string_view> toks = absl::StrSplit(line, sep);
      if (this->_ncols == 0) this->_ncols = toks.size();
      if (this->_ncols != toks.size()) {
        throw std::runtime_error(
            absl::StrFormat("Expected %lu columns, got %lu", this->_ncols, toks.size()));
      }
    }
    this->_matrix = std::vector<I>(this->_ncols * this->_nrows, 0);

    if (!in.eof() || in.bad() || fp.bad()) throw std::runtime_error("General IO error");

    fp.clear();
    fp.seekg(0);
    in.reset();
    in.push(boost::iostreams::bzip2_decompressor());
    in.push(fp);
    uint64_t j, n;
    for (auto i = 0UL; std::getline(in, line); ++i) {
      j = 0;
      for (const auto &tok : absl::StrSplit(line, sep)) {
        utils::parse_numeric_or_throw(tok, n);
        this->at(i, j++) = n;
      }
    }
    if (!in.eof() || in.bad() || fp.bad()) throw std::runtime_error("General IO error");

  } catch (const boost::iostreams::bzip2_error &err) {
    throw std::runtime_error(absl::StrFormat("An error occurred while parsing file '%s': %s.",
                                             path_to_file, err.what()));
  }
}

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
  std::ofstream fp(path_to_file, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream out;
  std::string buff;
  if (this->_updates_missed > 0) {
    absl::FPrintF(stderr, "WARNING: There were %lu missed updates!\n", this->_updates_missed);
  }
  try {
    out.push(boost::iostreams::bzip2_compressor(boost::iostreams::bzip2_params(9)));
    out.push(fp);
    std::vector<I> row(this->_ncols, 0);
    for (auto i = 0UL; i < this->_nrows; ++i) {
      for (auto j = 0UL; j < this->_ncols; ++j) {
        assert(i * j < this->_matrix.size());
        row[j] = this->at(i, j);
      }
      buff = absl::StrJoin(row, "\t") + "\n";
      out.write(buff.data(), buff.size());
    }
    if (out.bad() || fp.bad()) throw std::runtime_error("General IO error");
  } catch (const boost::iostreams::bzip2_error &err) {
    throw std::runtime_error(absl::StrFormat("An error occurred while writing file '%s': %s.",
                                             path_to_file, err.what()));
  }
  out.tellp();
  fp.tellp();
  return std::make_pair(fp.tellp(), out.tellp());
}

template <typename I>
void ContactMatrix<I>::clear_missed_updates_counter() {
  this->_updates_missed = 0;
}

template <typename I>
std::pair<uint32_t, uint32_t> ContactMatrix<I>::write_full_matrix_to_tsv(
    const std::string &path_to_file) const {
  std::ofstream fp(path_to_file, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream out;
  std::string buff;
  if (this->_updates_missed > 0) {
    absl::FPrintF(stderr, "WARNING: There were %lu missed updates!\n", this->_updates_missed);
  }
  try {
    out.push(boost::iostreams::bzip2_compressor(boost::iostreams::bzip2_params(9)));
    out.push(fp);
    std::vector<I> row(this->_ncols, 0);
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
      out.write(buff.data(), buff.size());
    }
    if (out.bad() || fp.bad()) throw std::runtime_error("General IO error");
  } catch (const boost::iostreams::bzip2_error &err) {
    throw std::runtime_error(absl::StrFormat("An error occurred while writing file '%s': %s.",
                                             path_to_file, err.what()));
  }
  out.tellp();
  fp.tellp();
  return std::make_pair(fp.tellp(), out.tellp());
}

}  // namespace modle
