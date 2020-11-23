#include <bzlib.h>

#include <array>
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

#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "fmt/printf.h"
#include "modle/contacts.hpp"
#include "modle/utils.hpp"

namespace modle {
template <class I>
ContactMatrix<I>::ContactMatrix(std::size_t nrows, std::size_t ncols, bool fill_with_random_numbers)
    : _nrows(nrows), _ncols(ncols), _contacts(_nrows * _ncols, 0) {
  if (fill_with_random_numbers) {
    std::random_device rand_dev;
    std::mt19937_64 rand_eng{rand_dev()};
    std::uniform_int_distribution<I> dist;
    std::generate(this->_contacts.begin(), this->_contacts.end(), [&]() { return dist(rand_eng); });
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
    uint64_t n{};
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
            // TODO: Figure out how to improve readability by removing some static_casts
            const auto noise1 = (*noise_generator)(*rand_eng) / bin_size;
            const auto noise2 = (*noise_generator)(*rand_eng) / bin_size;
            const auto pos1 =
                static_cast<std::size_t>(std::clamp(std::round(noise1 + static_cast<double>(i)),
                                                    0.0, static_cast<double>(this->_nrows) - 1.0));
            const auto pos2 =
                static_cast<std::size_t>(std::clamp(std::round(noise2 + static_cast<double>(j)),
                                                    0.0, static_cast<double>(this->_ncols) - 1.0));
            this->at(pos1, pos2) += static_cast<I>(n);
          }
        } else {
          this->at(i, j) = static_cast<I>(n);
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
  if ((j * this->_nrows) + i > this->_contacts.size()) {
    throw std::runtime_error(fmt::format(
        "ContactMatrix::at tried to access element m[{}][{}] of a matrix of shape [{}][{}]! "
        "({} >= {})",
        i, j, this->n_rows(), this->n_cols(), (j * this->_nrows) + i, this->_contacts.size()));
  }
  return this->_contacts[(j * this->_nrows) + i];
}

template <typename I>
const I &ContactMatrix<I>::at(std::size_t i, std::size_t j) const {
  if ((j * this->_nrows) + i > this->_contacts.size()) {
    throw std::runtime_error(fmt::format(
        "ContactMatrix::at tried to access element m[{}][{}] of a matrix of shape [{}][{}]! "
        "({} >= {})",
        i, j, this->n_rows(), this->n_cols(), (j * this->_nrows) + i, this->_contacts.size()));
  }
  return this->_contacts[(j * this->_nrows) + i];
}

template <typename I>
template <typename I2>
void ContactMatrix<I>::set(std::size_t row, std::size_t col, I2 n) {
  static_assert(std::is_integral<I2>::value,
                "ContactMatrix<I>::set expects the parameter n to be an integer type.");
  if constexpr (std::is_signed<I2>::value && std::is_unsigned<I>::value) {
    if (n < 0) {
      throw std::logic_error(
          fmt::format("ContactMatrix<unsigned-like>::set was called with n < 0: n = {}", n));
    }
  }

  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  try {
    if (j > this->n_cols()) {
      throw std::runtime_error(fmt::format(FMT_STRING("j={} > ncols={})"), j, this->n_cols()));
    }
    if (i > this->n_rows()) {
      ++this->_updates_missed;
      return;
    }

    assert(i <= this->n_rows());
    assert(j <= this->n_cols());
    assert(n >= std::numeric_limits<I>::min());
    assert(n <= std::numeric_limits<I>::max());

    this->at(i, j) = static_cast<I>(n);
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
  if constexpr (std::is_signed<I2>::value) {
    if (n < 0) {
      throw std::logic_error(
          fmt::format("ContactMatrix<I>::increment was called with n < 0. Consider using "
                      "::decrement(..., -n) instead: n = {}",
                      n));
    }
  }
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  try {
    if (j > this->n_cols()) {
      throw std::runtime_error(fmt::format(FMT_STRING("j={} > ncols={})"), j, this->n_cols()));
    }
    if (i > this->n_rows()) {
      ++this->_updates_missed;
      return;
    }

    assert(i <= this->n_rows());
    assert(j <= this->n_cols());
    assert(this->at(i, j) <= std::numeric_limits<I>::max() - n);

    this->at(i, j) += static_cast<I>(n);
  } catch (const std::runtime_error &err) {
    throw std::logic_error(fmt::format(FMT_STRING("ContactMatrix::increment({}, {}, {}): {})"), row,
                                       col, n, err.what()));
  }
}

template <typename I>
template <typename I2>
void ContactMatrix<I>::decrement(std::size_t row, std::size_t col, I2 n) {
  static_assert(std::is_integral<I2>::value,
                "ContactMatrix<I>::decrement expects the parameter n to be an integer type.");
  if constexpr (std::is_signed<I2>::value) {
    if (n < 0) {
      throw std::logic_error(
          fmt::format("ContactMatrix<I>::decrement was called with n < 0. Consider using "
                      "::increment(..., -n) instead: n = {}",
                      n));
    }
  }
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  try {
    if (j > this->n_cols()) {
      throw std::runtime_error(fmt::format(FMT_STRING("j={} > ncols={})"), j, this->n_cols()));
    }
    if (i > this->n_rows()) {
      ++this->_updates_missed;
      return;
    }

    assert(i <= this->n_rows());
    assert(j <= this->n_cols());
    assert(n <= this->at(i, j));

    this->at(i, j) -= static_cast<I>(n);
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
      fmt::print(FMT_STRING("%s\n)"), absl::StrJoin(row, "\t"));
    }
  } else {
    std::vector<I> row(this->n_cols());
    for (auto i = 0UL; i < this->n_rows(); ++i) {
      for (auto j = 0UL; j < this->n_cols(); ++j) {
        row[j] = this->at(i, j);
      }
      fmt::print(FMT_STRING("%s\n)"), absl::StrJoin(row, "\t"));
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
  if (i >= this->n_rows() || j >= this->n_cols()) {
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
std::size_t ContactMatrix<I>::n_rows() const {
  return this->_nrows;
}

template <typename I>
std::size_t ContactMatrix<I>::n_cols() const {
  return this->_ncols;
}

template <typename I>
std::pair<uint32_t, uint32_t> ContactMatrix<I>::write_to_tsv(const std::string &path_to_file,
                                                             std::string_view header,
                                                             int bzip2_block_size) const {
  std::ofstream fp(path_to_file, std::ios_base::binary);
  assert(header.empty() || header.back() == '\n');
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
      out.write(header.data(), static_cast<long>(header.size()));
    }
    for (auto i = 0UL; i < this->_nrows; ++i) {
      for (auto j = 0UL; j < this->_ncols; ++j) {
        assert(i * j < this->_contacts.size());
        row[j] = this->at(i, j);
      }
      fmt::print(out, "{}\n", absl::StrJoin(row, "\t"));
      buff = absl::StrCat(absl::StrJoin(row, "\t"), "\n");
      out.write(buff.data(), static_cast<long>(buff.size()));
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
const std::vector<I> &ContactMatrix<I>::get_raw_count_vector() const {
  return this->_contacts;
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
    throw fmt::system_error(errno, "IO error while reading file '%s'", path_to_file);
  }
  if (rewind_file) {
    in.seekg(0);
    if (!in) {
      throw fmt::system_error(errno, "IO error while rewinding file '%s'", path_to_file);
    }
  }

  std::vector<std::string_view> toks = absl::StrSplit(buff, '\t');
  if (toks.size() != Header::N_OF_EXPECTED_TOKENS || buff.front() != '#') {
    throw std::runtime_error(
        fmt::format("Malformed header: header should have the following structure: "
                    "#chr_name\\tbin_size\\tstart\\tend\\tdiagonal_width: got '{}'",
                    buff));
  }
  header.chr_name = toks[0].substr(1);  // Remove the leading #
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
std::pair<uint32_t, uint32_t> ContactMatrix<I>::write_full_matrix_to_tsv(
    const std::string &path_to_file, int bzip2_block_size) const {
  std::ofstream fp(path_to_file, std::ios_base::binary);
  boost::iostreams::filtering_ostream out;
  std::string buff;
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
      buff = absl::StrJoin(row, "\t") + "\n";
      out.write(buff.data(), static_cast<long>(buff.size()));
      raw_size += buff.size();
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

}  // namespace modle
