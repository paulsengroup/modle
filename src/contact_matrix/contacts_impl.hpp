#pragma once

#include <absl/strings/match.h>
#include <absl/strings/str_format.h>
#include <absl/strings/str_join.h>
#include <absl/strings/str_split.h>
#include <absl/strings/strip.h>
#include <absl/types/span.h>
#include <bzlib.h>
#include <fmt/printf.h>

#include <array>
#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/process.hpp>
#include <cassert>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string_view>

#include "modle/contacts.hpp"
#include "modle/suppress_compiler_warnings.hpp"
#include "modle/utils.hpp"

namespace modle {
template <class I>
ContactMatrix<I>::ContactMatrix(std::size_t nrows, std::size_t ncols, bool fill_with_random_numbers)
    : _nrows(nrows), _ncols(ncols), _contacts(_nrows * _ncols, 0) {
  if (fill_with_random_numbers) {
    std::random_device rand_dev;
    std::mt19937_64 rand_eng{rand_dev()};
    std::uniform_int_distribution<I> dist{0, std::numeric_limits<I>::max()};
    for (auto i = 0UL; i < _ncols; ++i) {
      for (auto j = i; j < i + _nrows && j < _ncols; ++j) {
        this->set(i, j, dist(rand_eng));
      }
    }
  }
}

template <class I>
ContactMatrix<I>::ContactMatrix(const std::string &path_to_file,
                                std::normal_distribution<double> *noise_generator, uint64_t seed,
                                char sep)
    : _nrows(0), _ncols(0) {
  std::unique_ptr<std::mt19937_64> rand_eng{nullptr};
  if (noise_generator) {
    rand_eng = std::make_unique<std::mt19937_64>(seed);
  }

  std::ifstream fp(path_to_file, std::ios_base::binary);
  boost::iostreams::filtering_istream in;
  std::string line;

  try {
    in.push(boost::iostreams::bzip2_decompressor());
    in.push(fp);

    const auto header = this->parse_header(path_to_file, in, false);
    this->_nrows = header.nrows;
    this->_ncols = header.ncols;
    this->_contacts = std::vector<I>(this->_ncols * this->_nrows, 0);

    std::size_t j{};
    std::size_t i{};
    I n{};
    const auto bin_size = static_cast<double>(header.bin_size);

    for (i = 0; std::getline(in, line); ++i) {
      if (i > this->_nrows) {
        continue;
      }
      j = 0;
      for (const auto &tok : absl::StrSplit(line, sep)) {
        utils::parse_numeric_or_throw(tok, n);
        if (noise_generator != nullptr) {
          for (auto k = n; k > 0; --k) {
            const auto noise1 = (*noise_generator)(*rand_eng) / bin_size;
            const auto noise2 = (*noise_generator)(*rand_eng) / bin_size;
            DISABLE_WARNING_PUSH
            DISABLE_WARNING_CONVERSION
            const auto pos1 = static_cast<std::size_t>(
                std::clamp(std::round(noise1 + i), 0.0, this->_nrows - 1.0));
            const auto pos2 = static_cast<std::size_t>(
                std::clamp(std::round(noise2 + j), 0.0, this->_ncols - 1.0));
            DISABLE_WARNING_POP
            this->at(pos1, pos2)++;
          }
        } else {
          this->at(i, j) = n;
        }
        ++j;
      }
    }

    if (i != this->_nrows) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected {} rows, got {})"), this->_nrows, i));
    }

    if ((!in && !in.eof()) || (!fp && !fp.eof())) {
      throw fmt::system_error(errno, "IO error while reading file '{}'", path_to_file);
    }

  } catch (const boost::iostreams::bzip2_error &err) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while decompressing file '{}': {})"),
                    path_to_file, err.what()));
  } catch (const std::runtime_error &err) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("An error occurred while parsing file '{}': {})"), path_to_file, err.what()));
  }
}

template <class I>
ContactMatrix<I>::ContactMatrix(std::string_view path_to_file,
                                std::normal_distribution<double> *noise_generator, uint64_t seed,
                                char sep)
    : ContactMatrix<I>::ContactMatrix(std::string{path_to_file}, noise_generator, seed, sep) {}

// TODO Rewrite this function
template <class I>
ContactMatrix<I>::ContactMatrix(const std::string &path_to_file, std::size_t nrows,
                                std::size_t ncols, char sep)
    : _nrows(nrows), _ncols(ncols), _contacts(_nrows * _ncols, 0) {
  std::ifstream fp(path_to_file);
  std::string line;
  std::vector<std::string_view> toks;
  auto read_and_tokenize = [&]() {
    if (std::getline(fp, line)) {
      if (toks = absl::StrSplit(line, sep); toks.size() != 3) {
        throw std::runtime_error(fmt::format(
            "Malformed file: expected 3 fields, got {}: line that triggered the error: '{}'",
            toks.size(), line));
      }
      return true;
    }
    if (!fp.eof() || fp.bad()) {
      throw fmt::system_error(errno, "Unable to read from file '{}'", path_to_file);
    }
    return false;
  };

  try {
    read_and_tokenize();
    uint64_t offset{};
    uint64_t bin_size{};
    std::size_t i{};
    std::size_t j{};
    double contacts{};
    if (toks[0] != toks[1]) {
      throw std::runtime_error(
          "The first entry in this file is expected to report counts for the same bin!");
    }
    utils::parse_numeric_or_throw(toks[0], offset);
    utils::parse_real_or_throw(toks[2], contacts);
    this->set(0, 0, std::round<uint64_t>(contacts));

    read_and_tokenize();
    utils::parse_numeric_or_throw(toks[1], bin_size);
    assert(bin_size > offset);
    bin_size -= offset;

    do {
      utils::parse_numeric_or_throw(toks[0], i);
      utils::parse_numeric_or_throw(toks[1], j);
      utils::parse_real_or_throw(toks[2], contacts);
      i = (i - offset) / bin_size;
      j = (j - offset) / bin_size;
      this->set(i, j, std::round<uint64_t>(contacts));
    } while (read_and_tokenize());

    if (!fp && !fp.eof()) {
      throw fmt::system_error(errno, "IO error while reading file '{}'", path_to_file);
    }

  } catch (const boost::iostreams::bzip2_error &err) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while decompressing file '{}': {})"),
                    path_to_file, err.what()));
  } catch (const std::runtime_error &err) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("An error occurred while parsing file '{}': {})"), path_to_file, err.what()));
  }
}

template <typename I>
I &ContactMatrix<I>::at(std::size_t i, std::size_t j) {
#ifndef NDEBUG
  if ((j * this->_nrows) + i > this->_contacts.size()) {
    throw std::runtime_error(fmt::format(
        "ContactMatrix::at tried to access element m[{}][{}] of a matrix of shape [{}][{}]! "
        "({} >= {})",
        i, j, this->nrows(), this->ncols(), (j * this->_nrows) + i, this->_contacts.size()));
  }
#endif
  return this->_contacts[(j * this->_nrows) + i];
}

template <typename I>
const I &ContactMatrix<I>::at(std::size_t i, std::size_t j) const {
#ifndef NDEBUG
  if ((j * this->_nrows) + i > this->_contacts.size()) {
    throw std::runtime_error(fmt::format(
        "ContactMatrix::at tried to access element m[{}][{}] of a matrix of shape [{}][{}]! "
        "({} >= {})",
        i, j, this->nrows(), this->ncols(), (j * this->_nrows) + i, this->_contacts.size()));
  }
#endif
  return this->_contacts[(j * this->_nrows) + i];
}

template <typename I>
template <typename I2>
void ContactMatrix<I>::set(std::size_t row, std::size_t col, I2 n) {
  static_assert(std::is_integral<I2>::value,
                "ContactMatrix<I>::set expects the parameter n to be an integer type.");
#ifndef NDEBUG
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_COMPARE
  DISABLE_WARNING_BOOL_COMPARE
  if (n < std::numeric_limits<I>::min() || n > std::numeric_limits<I>::max()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("ContactMatrix<I>::set(row={}, col={}, n={}): Overflow detected: "
                               "n={} is outside the range of representable numbers ({}-{})"),
                    row, col, n, n, std::numeric_limits<I>::min(), std::numeric_limits<I>::max()));
  }
  DISABLE_WARNING_POP
#endif

  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  try {
    if (j > this->ncols()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("attempt to access an element past the end of the contact matrix: "
                                 "j={}; ncols={}: j > ncols"),
                      j, this->ncols()));
    }
    if (i > this->nrows()) {
      ++this->_updates_missed;
      return;
    }

    auto &m = this->at(i, j);
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_CONVERSION
    DISABLE_WARNING_SIGN_CONVERSION
    DISABLE_WARNING_SIGN_COMPARE
#ifndef NDEBUG
    assert(this->_tot_contacts >= m);
    if constexpr (std::is_signed_v<I2>) {
      if (n < 0) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Setting counts to a negative value (n={}) is not allowed"), n));
      }
      assert(this->_tot_contacts - m <
             std::numeric_limits<decltype(this->_tot_contacts)>::max() - n);
    }
#endif
    this->_tot_contacts -= m;
    this->_tot_contacts += n;
    m = n;
    DISABLE_WARNING_POP
  } catch (const std::runtime_error &err) {
    throw std::logic_error(
        fmt::format(FMT_STRING("ContactMatrix::set({}, {}, {}): {})"), row, col, n, err.what()));
  }
}

template <typename I>
template <typename I2>
void ContactMatrix<I>::increment(std::size_t row, std::size_t col, I2 n) {
  static_assert(std::is_integral<I2>::value,
                "ContactMatrix<I>::increment expects the parameter n to be an integer type.");

  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  try {
    if (j > this->ncols()) {
      throw std::runtime_error(fmt::format(FMT_STRING("j={} > ncols={})"), j, this->ncols()));
    }
    if (i > this->nrows()) {
      ++this->_updates_missed;
      return;
    }

    DISABLE_WARNING_PUSH
    DISABLE_WARNING_SIGN_COMPARE
    DISABLE_WARNING_SIGN_CONVERSION
#ifndef NDEBUG
    auto &m = this->at(i, j);
    if constexpr (std::is_signed<I2>::value) {
      if (n < 0) {
        throw std::logic_error(fmt::format(
            FMT_STRING("ContactMatrix<I>::increment(row={}, col={}, n={}): consider using "
                       "ContactMatrix<I>::decrement instead of incrementing by a negative number."),
            row, col, n));
      }
    }
    if (std::numeric_limits<I>::max() - n < m) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Overflow detected: "
                     "incrementing m={} by n={} would result in a number outside of range {}-{}"),
          m, n, std::numeric_limits<I>::min(), std::numeric_limits<I>::max()));
    }
    if (std::numeric_limits<I>::max() - n < this->_tot_contacts) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Overflow detected: "
                     "incrementing _tot_contacts={} by n={} would result in a number outside of "
                     "range {}-{}"),
          this->_tot_contacts, n, std::numeric_limits<I>::min(), std::numeric_limits<I>::max()));
    }
#endif

    this->at(i, j) += n;
    this->_tot_contacts += n;
    DISABLE_WARNING_POP
  } catch (const std::runtime_error &err) {
    throw std::logic_error(
        fmt::format(FMT_STRING("ContactMatrix<I>::increment(row={}, col={}, n={}): {}"), row, col,
                    n, err.what()));
  }
}

template <typename I>
template <typename I2>
void ContactMatrix<I>::decrement(std::size_t row, std::size_t col, I2 n) {
  static_assert(std::is_integral<I2>::value,
                "ContactMatrix<I>::decrement expects the parameter n to be an integer type.");
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  try {
    if (j > this->ncols()) {
      throw std::runtime_error(fmt::format(FMT_STRING("j={} > ncols={})"), j, this->ncols()));
    }
    if (i > this->nrows()) {
      ++this->_updates_missed;
      return;
    }

    DISABLE_WARNING_PUSH
    DISABLE_WARNING_SIGN_COMPARE
    DISABLE_WARNING_SIGN_CONVERSION
#ifndef NDEBUG
    auto &m = this->at(i, j);
    if constexpr (std::is_signed<I2>::value) {
      if (n < 0) {
        throw std::logic_error(fmt::format(
            FMT_STRING("ContactMatrix<I>::decrement(row={}, col={}, n={}): consider using "
                       "ContactMatrix<I>::increment instead of decrementing by a negative number."),
            row, col, n));
      }
    }

    if (std::numeric_limits<I>::min() + n > m) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Overflow detected: "
                     "decrementing m={} by n={} would result in a number outside of range {}-{}"),
          m, n, std::numeric_limits<I>::min(), std::numeric_limits<I>::max()));
    }
    if (std::numeric_limits<I>::min() + n > this->_tot_contacts) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Overflow detected: "
                     "decrementing _tot_contacts={} by n={} would result in a number outside of "
                     "range {}-{}"),
          this->_tot_contacts, n, std::numeric_limits<I>::min(), std::numeric_limits<I>::max()));
    }
#endif

    this->at(i, j) -= n;
    this->_tot_contacts -= n;
    DISABLE_WARNING_POP
  } catch (const std::runtime_error &err) {
    throw std::logic_error(fmt::format(FMT_STRING("ContactMatrix::decrement({}, {}, {}): {})"), row,
                                       col, n, err.what()));
  }
}

template <typename I>
void ContactMatrix<I>::increment(std::size_t row, std::size_t col) {
  this->increment(row, col, static_cast<uint8_t>(1));
}

template <typename I>
void ContactMatrix<I>::decrement(std::size_t row, std::size_t col) {
  this->decrement(row, col, static_cast<uint8_t>(1));
}

template <typename I>
void ContactMatrix<I>::print(bool full) const {
  if (full) {
    std::vector<I> row(this->_ncols, 0);
    for (std::size_t y = 0; y < this->_ncols; ++y) {
      std::fill(row.begin(), row.end(), 0);
      for (std::size_t x = 0; x < this->_ncols; ++x) {
        auto j = x;
        auto i = j - y;
        if (y > x) {
          j = y;
          i = j - x;
        }
        if (i >= this->_nrows) {
          row[x] = 0;
        } else {
          row[x] = this->at(i, j);
        }
      }
      fmt::print(FMT_STRING("{}\n"), absl::StrJoin(row, "\t"));
    }
  } else {
    std::vector<I> row(this->ncols());
    for (auto i = 0UL; i < this->nrows(); ++i) {
      for (auto j = 0UL; j < this->ncols(); ++j) {
        row[j] = this->at(i, j);
      }
      fmt::print(FMT_STRING("{}\n"), absl::StrJoin(row, "\t"));
    }
  }
}

template <typename I>
I ContactMatrix<I>::get(std::size_t row, std::size_t col) const {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
#ifndef NDEBUG
  if (i >= this->ncols() || j >= this->ncols()) {
    throw std::logic_error(fmt::format(
        FMT_STRING("ContactMatrix<I>::get(row={}, col={}) tried to access an element outside of "
                   "the space {}x{}, which is what this contact matrix is supposed to represent"),
        row, col, col, col));
  }
#endif

  if (i >= this->nrows()) {
    return 0;
  }

  return this->at(i, j);
}

template <typename I>
std::vector<std::vector<I>> ContactMatrix<I>::generate_symmetric_matrix() const {
  std::vector<std::vector<I>> m;
  m.reserve(this->_ncols);
  for (std::size_t y = 0; y < this->_ncols; ++y) {
    std::vector<I> row(this->_ncols, 0);
    for (std::size_t x = 0; x < this->_ncols; ++x) {
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
std::size_t ContactMatrix<I>::nrows() const {
  return this->_nrows;
}

template <typename I>
std::size_t ContactMatrix<I>::ncols() const {
  return this->_ncols;
}

template <typename I>
std::size_t ContactMatrix<I>::npixels() const {
  return this->_nrows * this->_ncols;
}

template <typename I>
std::size_t ContactMatrix<I>::npixels_after_masking() const {
  auto npixels = this->npixels();
  const auto mask = this->generate_mask_for_bins_without_contacts();
  if (mask.all()) {
    return npixels;
  }
  if (mask.none()) {
    return 0;
  }

  auto count_zeros = [&mask](std::size_t start, std::size_t end) {
    assert(start <= end);
    std::size_t n = 0;
    while (start < end) {
      n += !mask[start++];
    }
    return n;
  };

  assert(this->nrows() <= this->ncols());
  for (auto i = 0UL; i < this->ncols(); ++i) {
    if (!mask[i]) {
      if (i < this->nrows()) {  // Upper left corner of cmatrix
        npixels -= this->nrows() - count_zeros(0UL, i);
        npixels -= i;
        assert(npixels <= this->npixels());
      } else if (i > this->ncols() - this->nrows()) {  // Lower right corner of cmatrix
        npixels -= this->nrows() - count_zeros(i - this->nrows(), i);
        npixels -= this->ncols() - i;
        assert(npixels <= this->npixels());
      } else {  // Middle of cmatrix
        npixels -= (2 * this->nrows()) - 1 - count_zeros(i - this->nrows(), i);
        assert(npixels <= this->npixels());
      }
    }
  }
  return npixels;
}

template <typename I>
std::pair<uint32_t, uint32_t> ContactMatrix<I>::write_to_tsv(const std::string &path_to_file,
                                                             std::string_view header,
                                                             int bzip2_block_size) const {
  std::ofstream fp(path_to_file, std::ios_base::binary);
  boost::iostreams::filtering_ostream out;
  std::string buff;
  uint64_t raw_size = 0;
  if (this->_updates_missed > 0) {
    fmt::print(stderr, "WARNING: There were {} missed updates!\n", this->_updates_missed);
  }
  try {
    out.push(boost::iostreams::bzip2_compressor(boost::iostreams::bzip2_params(bzip2_block_size)));
    out.push(fp);
    std::vector<I> row(this->_ncols, 0);
    if (!header.empty()) {
      if (absl::EndsWith(header, "\n")) {
        fmt::print(out, header);
      } else {
        fmt::print(out, "{}\n", header);
      }
    }
    for (auto i = 0UL; i < this->_nrows; ++i) {
      for (auto j = 0UL; j < this->_ncols; ++j) {
        assert(i * j < this->_contacts.size());
        row[j] = this->at(i, j);
      }
      // TODO: Figure out how to get the raw size when using fmt::print
      // fmt::print(out, "{}\n", absl::StrJoin(row, "\t"));
      buff = absl::StrJoin(row, "\t");
      fmt::print(out, "{}\n", buff);
      raw_size += buff.size();
    }
    if (!out || !fp) {
      throw fmt::system_error(errno, "IO error while writing to file '{}'", path_to_file);
    }

  } catch (const boost::iostreams::bzip2_error &err) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while compressing file '{}': {})"), path_to_file,
                    err.what()));
  }

  out.reset();
  return std::make_pair(raw_size, std::filesystem::file_size(path_to_file));
}

template <typename I>
void ContactMatrix<I>::clear_missed_updates_counter() {
  this->_updates_missed = 0;
}

template <typename I>
absl::Span<const I> ContactMatrix<I>::get_raw_count_vector() const {
  return absl::MakeConstSpan(this->_contacts);
}

template <typename I>
uint64_t ContactMatrix<I>::get_tot_contacts() const {
  return this->_tot_contacts;
}

template <typename I>
double ContactMatrix<I>::get_avg_contact_density() const {
  return static_cast<double>(this->get_tot_contacts()) /
         static_cast<double>(this->_contacts.size());
}

template <typename I>
uint64_t ContactMatrix<I>::get_n_of_missed_updates() const {
  return this->_updates_missed;
}

template <typename I>
typename ContactMatrix<I>::Header ContactMatrix<I>::parse_header(std::string_view path_to_file) {
  std::ifstream fp(path_to_file.data(), std::ios_base::binary);
  if (!fp) {
    throw fmt::system_error(errno, "Unable to open file '{}' for reading", path_to_file);
  }
  boost::iostreams::filtering_istream in;
  in.push(boost::iostreams::bzip2_decompressor());
  in.push(fp);
  if (!in) {
    throw fmt::system_error(errno, "Unable to decompress file '{}'", path_to_file);
  }
  return parse_header(path_to_file, in, false);
}

template <typename I>
typename ContactMatrix<I>::Header ContactMatrix<I>::parse_header(
    std::string_view path_to_file, boost::iostreams::filtering_istream &in, bool rewind_file) {
  ContactMatrix<I>::Header header;
  std::string buff;
  if (!std::getline(in, buff)) {
    throw fmt::system_error(errno, "IO error while reading file '{}'", path_to_file);
  }
  if (rewind_file) {
    in.seekg(0);
    if (!in) {
      throw fmt::system_error(errno, "IO error while rewinding file '{}'", path_to_file);
    }
  }

  std::vector<std::string_view> toks = absl::StrSplit(buff, '\t');
  if (toks.size() != Header::N_OF_EXPECTED_TOKENS || !absl::StartsWith(buff, "#")) {
    throw std::runtime_error(
        fmt::format("Malformed header: header should have the following structure: "
                    "#chr_name\\tbin_size\\tstart\\tend\\tdiagonal_width: got '{}'",
                    buff));
  }

  header.chr_name = absl::StripPrefix(toks[0], "#");
  modle::utils::parse_numeric_or_throw(toks[1], header.bin_size);
  modle::utils::parse_numeric_or_throw(toks[2], header.start);
  modle::utils::parse_numeric_or_throw(toks[3], header.end);
  std::string err;
  if (header.end < header.start) {
    err = fmt::format(FMT_STRING("end position < start position in header '{}')"), buff);
  }
  if (header.bin_size > header.end - header.start) {
    absl::StrAppendFormat(&err, "%sbin_size > end - start in header '%s'", err.empty() ? "" : "\n",
                          buff);
  }
  modle::utils::parse_numeric_or_throw(toks[4], header.diagonal_width);
  if (header.bin_size > header.diagonal_width) {
    absl::StrAppendFormat(&err, "%sbin_size > diagonal_width in header '%s'",
                          err.empty() ? "" : "\n", buff);
  }

  if (!err.empty()) {
    throw std::runtime_error(fmt::format(FMT_STRING("Malformed header: {})"), err));
  }
  header.ncols = static_cast<std::size_t>(std::ceil(static_cast<double>(header.end - header.start) /
                                                    static_cast<double>(header.bin_size)));
  header.nrows = static_cast<std::size_t>(
      std::ceil(static_cast<double>(header.diagonal_width) / static_cast<double>(header.bin_size)));
  return header;
}

template <typename I>
uint64_t ContactMatrix<I>::get_matrix_size_in_bytes() const {
  return this->_contacts.size() * sizeof(I);
}

template <typename I>
double ContactMatrix<I>::get_matrix_size_in_mb() const {
  return static_cast<double>(this->get_matrix_size_in_bytes()) / 1.0e6;
}

template <typename I>
std::pair<uint32_t, uint32_t> ContactMatrix<I>::write_full_matrix_to_tsv(
    const std::string &path_to_file, int bzip2_block_size) const {
  std::ofstream fp(path_to_file, std::ios_base::binary);
  boost::iostreams::filtering_ostream out;
  uint64_t raw_size = 0;
  if (this->_updates_missed > 0) {
    fmt::fprintf(stderr, "WARNING: There were {} missed updates!\n", this->_updates_missed);
  }
  try {
    out.push(boost::iostreams::bzip2_compressor(boost::iostreams::bzip2_params(bzip2_block_size)));
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
        if (i >= this->_nrows) {
          row[x] = 0;
        } else {
          row[x] = this->at(i, j);
        }
      }
      fmt::print(out, "{}\n", absl::StrJoin(row, "\t"));
    }
    if (!out || !fp) {
      throw fmt::system_error(errno, "IO error while writing to file '{}'", path_to_file);
    }
  } catch (const boost::iostreams::bzip2_error &err) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("An error occurred while writing file '{}': {})"), path_to_file, err.what()));
  }
  out.reset();
  return std::make_pair(raw_size, std::filesystem::file_size(path_to_file));
}

template <typename I>
boost::dynamic_bitset<> ContactMatrix<I>::generate_mask_for_bins_without_contacts() const {
  boost::dynamic_bitset<> mask{};
  this->generate_mask_for_bins_without_contacts(mask);
  return mask;
}
template <typename I>
void ContactMatrix<I>::generate_mask_for_bins_without_contacts(
    boost::dynamic_bitset<> &mask) const {
  mask.resize(this->ncols());
  mask.reset();
  for (auto i = 0UL; i < this->ncols(); ++i) {
    // Set bitmask to 1 if row contains at least one non-zero value
    for (auto j = i; j < (i + this->nrows()) && j < this->ncols(); ++j) {
      if ((mask[i] |= this->get(i, j))) {
        break;
      }
    }
    // Set bitmask to 1 if column "above" the current bin contains at least one non-zero value
    for (auto j = i; j > 0 && j > (i - this->nrows()); --j) {
      if ((mask[i] |= this->get(i, j))) {
        break;
      }
    }
  }
}

template <typename I>
void ContactMatrix<I>::reset() {
  std::fill(this->_contacts.begin(), this->_contacts.end(), 0);
  this->_tot_contacts = 0;
  this->_updates_missed = 0;
}

template <typename I>
bool ContactMatrix<I>::empty() const {
  return this->_tot_contacts == 0;
}

}  // namespace modle
