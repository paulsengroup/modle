#include <bzlib.h>

#include <array>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
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
#include "modle/contacts.hpp"
#include "modle/utils.hpp"

namespace modle {
template <class I>
ContactMatrix<I>::ContactMatrix(uint64_t nrows, uint64_t ncols, bool fill_with_random_numbers)
    : _nrows(nrows), _ncols(ncols), _contacts(_nrows * _ncols, 0) {
  if (fill_with_random_numbers) {
    std::random_device rand_dev;
    std::mt19937_64 rand_eng{rand_dev()};
    std::uniform_int_distribution<I> dist;
    std::generate(this->_contacts.begin(), this->_contacts.end(), [&]() { return dist(rand_eng); });
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
    this->_contacts = std::vector<I>(this->_ncols * this->_nrows, 0);

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
    throw std::runtime_error(absl::StrFormat("An error occurred while decompressing file '%s': %s.",
                                             path_to_file, err.what()));
  } catch (const std::runtime_error &err) {
    throw std::runtime_error(absl::StrFormat("An error occurred while parsing file '%s': %s.",
                                             path_to_file, err.what()));
  }
}

template <class I>
ContactMatrix<I>::ContactMatrix(const std::string &path_to_file, uint64_t nrows, uint64_t ncols,
                                char sep)
    : _nrows(nrows), _ncols(ncols), _contacts(_nrows * _ncols, 0) {
  std::ifstream fp(path_to_file);
  std::string line;
  std::array<char, 8'192> buff{};
  std::vector<std::string_view> toks;
  fp.rdbuf()->pubsetbuf(buff.data(), buff.size());

  auto read_and_tokenize = [&]() {
    if (std::getline(fp, line)) {
      if (toks = absl::StrSplit(line, sep); toks.size() != 3) {
        throw std::runtime_error(absl::StrFormat(
            "Malformed file: expected 3 fields, got %lu: line that triggered the error: '%s'",
            toks.size(), line));
      }
      return true;
    }
    if (!fp.eof() || fp.bad()) {
      throw std::runtime_error("Unable to read from file.");
    }
    return false;
  };

  try {
    read_and_tokenize();
    uint64_t offset, bin_size, i, j;
    double contacts;
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

    if (!fp.eof() || fp.bad()) throw std::runtime_error("General IO error");

  } catch (const boost::iostreams::bzip2_error &err) {
    throw std::runtime_error(absl::StrFormat("An error occurred while decompressing file '%s': %s.",
                                             path_to_file, err.what()));
  } catch (const std::runtime_error &err) {
    throw std::runtime_error(absl::StrFormat("An error occurred while parsing file '%s': %s.",
                                             path_to_file, err.what()));
  }
}

template <typename I>
I &ContactMatrix<I>::at(uint64_t i, uint64_t j) {
  if ((j * this->_nrows) + i > this->_contacts.size()) {
    throw std::runtime_error(absl::StrFormat(
        "ContactMatrix::at tried to access element m[%lu][%lu] of a matrix of shape [%lu][%lu]! "
        "(%lu >= %lu).",
        i, j, this->n_rows(), this->n_cols(), (j * this->_nrows) + i, this->_contacts.size()));
  }
  return this->_contacts[(j * this->_nrows) + i];
}

template <typename I>
const I &ContactMatrix<I>::at(uint64_t i, uint64_t j) const {
  if ((j * this->_nrows) + i > this->_contacts.size()) {
    throw std::runtime_error(absl::StrFormat(
        "ContactMatrix::at tried to access element m[%lu][%lu] of a matrix of shape [%lu][%lu]! "
        "(%lu >= %lu).",
        i, j, this->n_rows(), this->n_cols(), (j * this->_nrows) + i, this->_contacts.size()));
  }
  return this->_contacts[(j * this->_nrows) + i];
}

template <typename I>
void ContactMatrix<I>::increment(uint64_t row, uint64_t col, I n) {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  try {
    if (j > this->n_cols())
      throw std::runtime_error(absl::StrFormat("j=%d > ncols=%d.", j, this->n_cols()));
    if (i > this->n_rows()) {
      ++this->_updates_missed;
      return;
    }

    assert(i <= this->n_rows());
    assert(j <= this->n_cols());

    this->at(i, j) += n;
  } catch (const std::runtime_error &err) {
    throw std::logic_error(
        absl::StrFormat("ContactMatrix::increment(%d, %d, %d): %s", row, col, n, err.what()));
  }
}

template <typename I>
void ContactMatrix<I>::set(uint64_t row, uint64_t col, I n) {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  try {
    if (j > this->n_cols())
      throw std::runtime_error(absl::StrFormat("j=%d > ncols=%d.", j, this->n_cols()));
    if (i > this->n_rows()) {
      ++this->_updates_missed;
      return;
    }

    assert(i <= this->n_rows());
    assert(j <= this->n_cols());
    this->at(i, j) = n;
  } catch (const std::runtime_error &err) {
    throw std::logic_error(
        absl::StrFormat("ContactMatrix::set(%d, %d, %d): %s", row, col, n, err.what()));
  }
}

template <typename I>
void ContactMatrix<I>::decrement(uint64_t row, uint64_t col, I n) {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  try {
    if (j > this->n_cols())
      throw std::runtime_error(absl::StrFormat("j=%d > ncols=%d.", j, this->n_cols()));
    if (i > this->n_rows()) {
      ++this->_updates_missed;
      return;
    }

    assert(i <= this->n_rows());
    assert(j <= this->n_cols());
    assert(n <= this->at(i, j));
    this->at(i, j) -= n;
  } catch (const std::runtime_error &err) {
    throw std::logic_error(
        absl::StrFormat("ContactMatrix::decrement(%d, %d, %d): %s", row, col, n, err.what()));
  }
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
    std::vector<I> row(this->n_cols());
    for (auto i = 0UL; i < this->n_rows(); ++i) {
      for (auto j = 0UL; j < this->n_cols(); ++j) {
        row[j] = this->at(i, j);
      }
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

      if (i < this->_nrows) row[x] = this->at(i, j);
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
  uint64_t raw_size = 0;
  if (this->_updates_missed > 0) {
    absl::FPrintF(stderr, "WARNING: There were %lu missed updates!\n", this->_updates_missed);
  }
  try {
    out.push(boost::iostreams::bzip2_compressor(boost::iostreams::bzip2_params(9)));
    out.push(fp);
    std::vector<I> row(this->_ncols, 0);
    for (auto i = 0UL; i < this->_nrows; ++i) {
      for (auto j = 0UL; j < this->_ncols; ++j) {
        assert(i * j < this->_contacts.size());
        row[j] = this->at(i, j);
      }
      buff = absl::StrJoin(row, "\t") + "\n";
      out.write(buff.data(), buff.size());
      raw_size += buff.size();
    }
    if (out.bad() || fp.bad()) throw std::runtime_error("General IO error");
  } catch (const boost::iostreams::bzip2_error &err) {
    throw std::runtime_error(absl::StrFormat("An error occurred while writing file '%s': %s.",
                                             path_to_file, err.what()));
  }
  out.flush();
  return std::make_pair(raw_size, std::filesystem::file_size(path_to_file));
}

template <typename I>
std::pair<uint32_t, uint32_t> ContactMatrix<I>::write_to_hic(const std::string &path_to_juicer,
                                                             const std::string &chr_name,
                                                             uint64_t bin_size,
                                                             const std::string &path_to_chr_sizes,
                                                             const std::string &path_to_file,
                                                             const std::string &tmp_dir) const {
  std::string tmp_file = absl::StrFormat("%s/%s.juicer.tmp.gz", tmp_dir, chr_name);
  std::ofstream fp(tmp_file, std::ios_base::out | std::ios_base::binary);
  boost::process::ipstream juicer_stderr;
  std::string buff;
  uint64_t raw_size = 0;
  if (this->_updates_missed > 0) {
    absl::FPrintF(stderr, "WARNING: There were %lu missed updates!\n", this->_updates_missed);
  }
  try {
    boost::iostreams::filtering_ostream juicer_in;
    juicer_in.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(1)));
    juicer_in.push(fp);
    absl::FPrintF(stderr, "Writing contacts to feed to Juicer pre...");
    // TODO: The -1 here is a workaround to discard the last bin (which usually does not have the
    // same size as the others). There's probably a better way to deal with this
    for (auto i = 0UL; i < this->_ncols - 1; ++i) {
      for (auto j = 0UL; j < this->_ncols - 1; ++j) {
        if (auto n = this->get(i, j); n != 0) {
          // strand1 chr1 pos1 frag1 strand2 chr2 pos2 frag2
          uint64_t pos1 = ((2 * i * bin_size) + bin_size) / 2;
          uint64_t pos2 = ((2 * j * bin_size) + bin_size) / 2;
          buff = absl::StrFormat("1 %s %lu 0 0 %s %lu 1\n", chr_name, pos1, chr_name, pos2);
          for (auto k = this->get(i, j); k > 0; --k) {
            juicer_in.write(buff.data(), buff.size());
            raw_size += buff.size();
            if (juicer_in.bad() || fp.bad()) throw std::runtime_error("General IO error");
          }
        }
      }
    }
    if (juicer_in.bad() || fp.bad()) throw std::runtime_error("General IO error");
    juicer_in.flush();
  } catch (const boost::iostreams::gzip_error &err) {
    std::filesystem::remove(tmp_file);
    throw std::runtime_error(absl::StrFormat("An error occurred while writing file '%s': %s.",
                                             path_to_file, err.what()));
  } catch (const std::runtime_error &err) {
    std::filesystem::remove(tmp_file);
    throw std::runtime_error(absl::StrFormat("An error occurred while writing file '%s': %s.",
                                             path_to_file, err.what()));
  }

  absl::FPrintF(stderr, " DONE!\n");
  std::vector<std::string> args{boost::process::search_path("java").string(),
                                "-Xms512m",
                                "-Xmx2048m",
                                "-jar",
                                path_to_juicer,
                                "pre",
                                tmp_file,
                                path_to_file,
                                path_to_chr_sizes};
  absl::FPrintF(stderr, "Running: %s\n", absl::StrJoin(args, " "));
  boost::process::child juicer(args, boost::process::std_in.close(),
                               boost::process::std_out > boost::process::null,
                               boost::process::std_err > juicer_stderr);

  juicer.wait();
  std::filesystem::remove(tmp_file);
  if (auto ec = juicer.exit_code(); ec != 0) {
    for (; std::getline(juicer_stderr, buff);) absl::FPrintF(stderr, "%s", buff);
    juicer_stderr.close();
    throw std::runtime_error(absl::StrFormat("Juicer exited with error %lu.", ec));
  }
  juicer_stderr.close();
  assert(!std::filesystem::exist(tmp_file));
  return {raw_size, std::filesystem::file_size(path_to_file.data())};
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
std::pair<uint32_t, uint32_t> ContactMatrix<I>::write_full_matrix_to_tsv(
    const std::string &path_to_file) const {
  std::ofstream fp(path_to_file, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream out;
  std::string buff;
  uint64_t raw_size = 0;
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
      raw_size += buff.size();
    }
    if (out.bad() || fp.bad()) throw std::runtime_error("General IO error");
  } catch (const boost::iostreams::bzip2_error &err) {
    throw std::runtime_error(absl::StrFormat("An error occurred while writing file '%s': %s.",
                                             path_to_file, err.what()));
  }
  out.flush();
  return std::make_pair(raw_size, std::filesystem::file_size(path_to_file));
}

}  // namespace modle
