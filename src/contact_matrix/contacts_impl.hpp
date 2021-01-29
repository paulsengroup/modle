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

template <typename I>
I ContactMatrix<I>::get(std::size_t row, std::size_t col) const {
  const auto [i, j] = transpose_coords(row, col);
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
  const auto [i, j] = transpose_coords(row, col);

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
void ContactMatrix<I>::add(std::size_t row, std::size_t col, I2 n) {
  static_assert(std::is_integral<I2>::value,
                "ContactMatrix<I>::add expects the parameter n to be an integer type.");

  const auto [i, j] = transpose_coords(row, col);

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
    throw std::logic_error(fmt::format(
        FMT_STRING("ContactMatrix<I>::add(row={}, col={}, n={}): {}"), row, col, n, err.what()));
  }
}

template <typename I>
template <typename I2>
void ContactMatrix<I>::subtract(std::size_t row, std::size_t col, I2 n) {
  static_assert(std::is_integral<I2>::value,
                "ContactMatrix<I>::subtract expects the parameter n to be an integer type.");

  const auto [i, j] = transpose_coords(row, col);

  try {
    if (j > this->ncols()) {
      throw std::runtime_error(fmt::format(FMT_STRING("j={} > ncols={})"), j, this->ncols()));
    }
    if (i > this->nrows()) {
      ++this->_updates_missed;
      return;
    }

    DISABLE_WARNING_PUSH
    DISABLE_WARNING_CONVERSION
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
    throw std::logic_error(fmt::format(FMT_STRING("ContactMatrix::subtract({}, {}, {}): {})"), row,
                                       col, n, err.what()));
  }
}

template <typename I>
void ContactMatrix<I>::increment(std::size_t row, std::size_t col) {
  this->add(row, col, 1U);
}

template <typename I>
void ContactMatrix<I>::decrement(std::size_t row, std::size_t col) {
  this->subtract(row, col, 1U);
}

template <typename I>
std::size_t ContactMatrix<I>::ncols() const {
  return this->_ncols;
}

template <typename I>
std::size_t ContactMatrix<I>::nrows() const {
  return this->_nrows;
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
      // We are processing pixels in the upper left corner of cmatrix
      if (i < this->nrows()) {
        npixels -= this->nrows() - count_zeros(0UL, i);
        npixels -= i;
        assert(npixels <= this->npixels());
        // We are processing pixels in the lower right corner of cmatrix
      } else if (i > this->ncols() - this->nrows()) {
        npixels -= this->nrows() - count_zeros(i - this->nrows(), i);
        npixels -= this->ncols() - i;
        assert(npixels <= this->npixels());
      } else {  // We are processing pixels in the center of the cmatrix
        npixels -= (2 * this->nrows()) - 1 - count_zeros(i - this->nrows(), i);
        assert(npixels <= this->npixels());
      }
    }
  }
  return npixels;
}

template <typename I>
uint64_t ContactMatrix<I>::get_n_of_missed_updates() const {
  return this->_updates_missed;
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
uint64_t ContactMatrix<I>::get_matrix_size_in_bytes() const {
  return this->_contacts.size() * sizeof(I);
}

template <typename I>
double ContactMatrix<I>::get_matrix_size_in_mb() const {
  return static_cast<double>(this->get_matrix_size_in_bytes()) / 1.0e6;
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
boost::dynamic_bitset<> ContactMatrix<I>::generate_mask_for_bins_without_contacts() const {
  boost::dynamic_bitset<> mask{};
  this->generate_mask_for_bins_without_contacts(mask);
  return mask;
}

template <typename I>
void ContactMatrix<I>::clear_missed_updates_counter() {
  this->_updates_missed = 0;
}

template <typename I>
void ContactMatrix<I>::reset() {
  std::fill(this->_contacts.begin(), this->_contacts.end(), 0);
  this->_tot_contacts = 0;
  this->_updates_missed = 0;
}

template <typename I>
bool ContactMatrix<I>::empty() const {
  assert(std::all_of(this->_contacts.begin(), this->_contacts.end(),
                     [](const auto n) { return n == 0; }));
  return this->_tot_contacts == 0;
}

template <typename I>
absl::Span<const I> ContactMatrix<I>::get_raw_count_vector() const {
  return absl::MakeConstSpan(this->_contacts);
}

template <typename I>
void ContactMatrix<I>::compute_row_wise_contact_histogram(std::vector<uint64_t> &buff) const {
  buff.resize(this->nrows());
  std::fill(buff.begin(), buff.end(), 0);

  for (auto i = 0UL; i < this->ncols(); ++i) {
    for (auto j = i; j < i + this->nrows() && j < this->ncols(); ++j) {
      // j - i corresponds to the distance from the diagonal
      buff[j - i] += this->get(j, i);
    }
  }
}

template <typename I>
std::vector<uint64_t> ContactMatrix<I>::compute_row_wise_contact_histogram() const {
  std::vector<uint64_t> buff(this->nrows());
  this->compute_row_wise_contact_histogram(buff);
  return buff;
}

template <typename I>
void ContactMatrix<I>::deplete_contacts(double depletion_multiplier) {
  const auto hist = this->compute_row_wise_contact_histogram();
  const auto effective_nbins = this->generate_mask_for_bins_without_contacts().count();
  // This histogram contains the average contact number (instead of the total)
  std::vector<uint64_t> row_wise_avg_contacts(hist.size());
  std::transform(hist.begin(), hist.end(), row_wise_avg_contacts.begin(), [&](const auto n) {
    return static_cast<uint64_t>(std::round((depletion_multiplier * static_cast<double>(n)) /
                                            static_cast<double>(effective_nbins)));
  });

  for (auto i = 0UL; i < this->ncols(); ++i) {
    for (auto j = i; j < i + this->nrows() && j < this->ncols(); ++j) {
      // j - i corresponds to the distance from the diagonal
      if (this->get(j, i) > row_wise_avg_contacts[j - i]) {
        this->subtract(j, i, row_wise_avg_contacts[j - i]);
      } else {
        this->set(j, i, 0);
      }
    }
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
std::pair<std::size_t, std::size_t> ContactMatrix<I>::transpose_coords(std::size_t row,
                                                                       std::size_t col) {
  auto j = col;
  auto i = j - row;
  if (row > col) {
    j = row;
    i = j - col;
  }
  return std::make_pair(i, j);
}

}  // namespace modle
