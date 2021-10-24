// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/strings/str_join.h>        // for StrJoin
#include <absl/types/span.h>              // for Span, MakeConstSpan
#include <cpp-sort/sorter_facade.h>       // for sorter_facade
#include <cpp-sort/sorters/ska_sorter.h>  // for ska_sort, ska_sorter
#include <fmt/format.h>                   // for FMT_STRING, join
#include <fmt/ostream.h>

#include <algorithm>                                // for clamp, min, fill, equal_range
#include <atomic>                                   // for memory_order_relaxed
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cassert>                                  // for assert
#include <cmath>                                    // for round
#include <iostream>                                 // for cout
#include <limits>                                   // for numeric_limits
#include <mutex>                                    // for mutex
#include <numeric>                                  // for accumulate
#include <stdexcept>                                // for runtime_error, logic_error
#include <type_traits>                              // for is_integral, is_signed, remove_re...
#include <utility>                                  // for pair, make_pair, pair<>::second
#include <vector>                                   // for vector, allocator

#include "modle/common/common.hpp"                      // for u64, i64
#include "modle/common/random.hpp"                      // for PRNG, uniform_int_distribution
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_PUSH, DISABLE_WAR...
#include "modle/common/utils.hpp"                       // for ndebug_defined, ndebug_not_defined

namespace modle {

#if defined(__clang__) && __clang_major__ < 7
template <class N>
ContactMatrix<N>::ContactMatrix(ContactMatrix<N> &&other) noexcept
    : _nrows(other.nrows()),
      _ncols(other.ncols()),
      _contacts(std::move(other._contacts)),
      _mtxes(_ncols),
      _tot_contacts(other.get_tot_contacts()),
      _updates_missed(other.get_n_of_missed_updates()) {}

#endif

template <class N>
ContactMatrix<N>::ContactMatrix(const ContactMatrix<N> &other)
    : _nrows(other.nrows()),
      _ncols(other.ncols()),
      _contacts(other._contacts),
      _mtxes(_ncols),
      _tot_contacts(other._tot_contacts.load()),
      _updates_missed(other._updates_missed.load()) {}

template <class N>
ContactMatrix<N>::ContactMatrix(const usize nrows, const usize ncols,
                                const bool fill_with_random_numbers, const u64 seed)
    : _nrows(std::min(nrows, ncols)),
      _ncols(ncols),
      _contacts(_nrows * _ncols + 1, N(0)),
      _mtxes(_ncols + 1) {
  if (fill_with_random_numbers) {
    auto rand_eng = random::PRNG(seed);

    auto dist = []() {
      if constexpr (std::is_floating_point_v<N>) {
        return random::uniform_real_distribution<N>{0, 65553};
      }
      uint64_t max_ = std::min(u64(65553), static_cast<u64>((std::numeric_limits<N>::max)()));
      return random::uniform_int_distribution<N>{0, static_cast<N>(max_)};
    }();
    for (usize i = 0; i < _ncols; ++i) {
      for (usize j = i; j < i + _nrows && j < _ncols; ++j) {
        this->set(i, j, dist(rand_eng));
      }
    }
  }
}

template <class N>
ContactMatrix<N>::ContactMatrix(const bp_t length, const bp_t diagonal_width, const bp_t bin_size,
                                const bool fill_with_random_numbers)
    : ContactMatrix((diagonal_width + bin_size - 1) / bin_size, (length + bin_size - 1) / bin_size,
                    fill_with_random_numbers) {}

template <class N>
ContactMatrix<N>::ContactMatrix(const absl::Span<const N> contacts, const usize nrows,
                                const usize ncols, const usize tot_contacts,
                                const usize updates_missed)
    : _nrows(nrows),
      _ncols(ncols),
      _contacts(contacts.begin(), contacts.end()),
      _mtxes(_ncols),
      _updates_missed(static_cast<i64>(updates_missed)),
      _tot_contacts(static_cast<i64>(tot_contacts)) {
  if (this->_tot_contacts == 0 && !this->_contacts.empty()) {
    this->_tot_contacts = std::accumulate(this->_contacts.begin(), this->_contacts.end(), usize(0));
  }
}

#if defined(__clang__) && __clang_major__ < 7
template <class N>
ContactMatrix<N> &ContactMatrix<N>::operator=(ContactMatrix<N> &&other) noexcept {
  if (this == &other) {
    return *this;
  }
  _nrows = other.nrows();
  _ncols = other.ncols();
  _contacts = std::move(other._contacts);
  _mtxes.resize(_ncols);
  _tot_contacts = other._tot_contacts.load();
  _updates_missed = other._updates_missed.load();

  return *this;
}
#endif

template <class N>
ContactMatrix<N> &ContactMatrix<N>::operator=(const ContactMatrix<N> &other) {
  if (this == &other) {
    return *this;
  }

  _nrows = other.nrows();
  _ncols = other.ncols();
  if (!_contacts.empty()) {
    _contacts.resize(other._contacts.size());
    std::copy(other._contacts.begin(), other._contacts.end(), _contacts.begin());
  } else {
    _contacts = other._contacts;
  }
  _mtxes.resize(_ncols);
  _tot_contacts = other._tot_contacts.load();
  _updates_missed = other._updates_missed.load();

  return *this;
}

template <class N>
N ContactMatrix<N>::unsafe_get(const usize row, const usize col) const
    noexcept(utils::ndebug_defined()) {
  const auto [i, j] = transpose_coords(row, col);
  if constexpr (utils::ndebug_not_defined()) {
    if (i >= this->ncols() || j >= this->ncols()) {
      throw std::logic_error(fmt::format(
          FMT_STRING("ContactMatrix<N>::get(row={}, col={}) tried to access an element outside of "
                     "the space {}x{}, which is what this contact matrix is supposed to represent"),
          row, col, col, col));
    }
  }

  if (i >= this->nrows()) {
    return 0;
  }

  return this->at(i, j);
}

// NOTE the pixels returned by this function go from the diagonal towards the periphery
// Example: given the following matrix
//          1  2  3
//          2  4  5
//          3  5  6
// Fetching col #3 would yield 6 5 3
template <class N>
void ContactMatrix<N>::unsafe_get_column(const usize col, std::vector<N> &buff,
                                         const usize row_offset) const
    noexcept(utils::ndebug_defined()) {
  assert(row_offset <= col);  // NOLINT
  const auto [rowt, colt] = transpose_coords(col - row_offset, col);

  const auto first_idx = (colt * this->nrows()) + rowt;
  const auto last_idx = std::min(first_idx + this->nrows() - row_offset, this->_contacts.size());

  buff.resize(last_idx - first_idx);

  const auto first = this->_contacts.begin() + static_cast<isize>(first_idx);
  const auto last = this->_contacts.begin() + static_cast<isize>(last_idx);
  std::copy(first, last, buff.begin());
}

// NOTE the pixels returned by this function go from the diagonal towards the periphery
// Example: given the following matrix
//          1  2  3
//          2  4  5
//          3  5  6
// Fetching col #1 would yield 1 2 3
template <class N>
void ContactMatrix<N>::unsafe_get_row(const usize row, std::vector<N> &buff,
                                      const usize col_offset) const
    noexcept(utils::ndebug_defined()) {
  assert(row >= col_offset);  // NOLINT
  buff.resize(std::min(row + 1, this->nrows()) - col_offset);

  for (usize i = 0; i < buff.size(); ++i) {
    assert(i < row + 1);  // NOLINT
    buff[i] = this->unsafe_get(row, row - col_offset - i);
  }
}

template <class N>
N ContactMatrix<N>::unsafe_get(const usize row, const usize col, const usize block_size) const
    noexcept(utils::ndebug_defined()) {
  assert(block_size > 0);              // NOLINT
  assert(block_size < this->nrows());  // NOLINT
  // For now we only support blocks with an odd size
  assert(block_size % 2 != 0);  // NOLINT
  if (block_size == 1) {
    return this->unsafe_get(row, col);
  }

  if constexpr (utils::ndebug_not_defined()) {
    const auto [i, j] = transpose_coords(row, col);
    if (i >= this->ncols() || j >= this->ncols()) {
      throw std::logic_error(fmt::format(
          FMT_STRING("ContactMatrix<N>::get(row={}, col={}) tried to access an element outside of "
                     "the space {}x{}, which is what this contact matrix is supposed to represent"),
          row, col, col, col));
    }
  }

  // Edges are handled like shown here: https://en.wikipedia.org/wiki/File:Extend_Edge-Handling.png
  const auto bs = static_cast<i64>(block_size);
  const auto first_row = static_cast<i64>(row) - (bs / 2);
  const auto first_col = static_cast<i64>(col) - (bs / 2);
  N n{0};

  for (auto i = first_row; i < first_row + bs; ++i) {
    for (auto j = first_col; j < first_col + bs; ++j) {
      const auto ii = static_cast<usize>(std::clamp(i, i64(0), static_cast<i64>(this->_ncols - 1)));
      const auto jj = static_cast<usize>(std::clamp(j, i64(0), static_cast<i64>(this->_ncols - 1)));
      n += this->unsafe_get(ii, jj);
    }
  }
  return n;
}

template <class N>
void ContactMatrix<N>::unsafe_set(const usize row, const usize col,
                                  const N n) noexcept(utils::ndebug_defined()) {
  this->internal_set(row, col, n, nullptr);
}

template <class N>
void ContactMatrix<N>::set(const usize row, const usize col,
                           const N n) noexcept(utils::ndebug_defined()) {
  const auto [_, j] = transpose_coords(row, col);
  assert(j < this->_mtxes.size());  // NOLINT
  this->internal_set(row, col, n, &this->_mtxes[j]);
}

template <class N>
void ContactMatrix<N>::internal_set(const usize row, const usize col, const N n,
                                    std::mutex *mtx) noexcept(utils::ndebug_defined()) {
  const auto [i, j] = this->transpose_coords(row, col);
  if constexpr (utils::ndebug_not_defined()) {
    try {
      this->bound_check_column(j);
    } catch (const std::runtime_error &err) {
      throw std::logic_error(
          fmt::format(FMT_STRING("ContactMatrix::set({}, {}, {}): {})"), i, j, n, err.what()));
    }
  }

  if (i > this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, i64(1), std::memory_order_relaxed);
    return;
  }

  [[maybe_unused]] auto lck = [&]() {
    if (mtx) {
      return std::unique_lock(*mtx);
    }
    return std::unique_lock<std::mutex>{};
  }();
  auto &m = this->at(i, j);
  if (n > m) {
    this->_tot_contacts += n - m;
  } else {
    this->_tot_contacts -= m - n;
  }
  m = n;
}

template <class N>
void ContactMatrix<N>::add(const usize row, const usize col,
                           const N n) noexcept(utils::ndebug_defined()) {
  assert(n > 0);  // NOLINT Use subtract to add a negative number
  const auto [i, j] = transpose_coords(row, col);

  if constexpr (utils::ndebug_not_defined()) {
    try {
      this->bound_check_column(j);
      this->check_for_overflow_on_add(i, j, n);
    } catch (const std::runtime_error &err) {
      throw std::logic_error(
          fmt::format(FMT_STRING("ContactMatrix::add({}, {}, {}): {})"), row, col, n, err.what()));
    }
  }
  if (i > this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, i64(1), std::memory_order_relaxed);
    return;
  }

  {
    std::scoped_lock<std::mutex> lck(this->_mtxes[j]);
    this->at(i, j) += n;
  }
  this->_tot_contacts += n;
}

template <class N>
void ContactMatrix<N>::subtract(const usize row, const usize col,
                                const N n) noexcept(utils::ndebug_defined()) {
  assert(n >= 0);  // NOLINT Use add to subtract a negative number
  const auto [i, j] = transpose_coords(row, col);

  if constexpr (utils::ndebug_not_defined()) {
    try {
      this->bound_check_column(j);
      this->check_for_overflow_on_subtract(i, j, n);
    } catch (const std::runtime_error &err) {
      throw std::logic_error(fmt::format(FMT_STRING("ContactMatrix::subtract({}, {}, {}): {})"),
                                         row, col, n, err.what()));
    }
  }

  if (i > this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, i64(1), std::memory_order_relaxed);
    return;
  }

  {
    std::scoped_lock<std::mutex> lck(this->_mtxes[j]);
    this->at(i, j) -= n;
  }
  this->_tot_contacts -= n;
}

template <class N>
void ContactMatrix<N>::increment(usize row, usize col) noexcept(utils::ndebug_defined()) {
  this->add(row, col, N(1));
}

template <class N>
void ContactMatrix<N>::decrement(usize row, usize col) noexcept(utils::ndebug_defined()) {
  this->subtract(row, col, N(1));
}

template <class N>
constexpr usize ContactMatrix<N>::ncols() const noexcept(utils::ndebug_defined()) {
  return this->_ncols;
}

template <class N>
constexpr usize ContactMatrix<N>::nrows() const noexcept(utils::ndebug_defined()) {
  return this->_nrows;
}

template <class N>
constexpr usize ContactMatrix<N>::npixels() const noexcept(utils::ndebug_defined()) {
  return this->_nrows * this->_ncols;
}

template <class N>
usize ContactMatrix<N>::unsafe_npixels_after_masking() const {
  auto npixels = this->npixels();
  const auto mask = this->unsafe_generate_mask_for_bins_without_contacts();
  if (mask.all()) {
    return npixels;
  }
  if (mask.none()) {
    return 0;
  }

  auto count_zeros = [&mask](usize start, usize end) {
    assert(start <= end);
    usize n = 0;
    while (start < end) {
      n += !mask[start++];
    }
    return n;
  };

  assert(this->nrows() <= this->ncols());
  for (usize i = 0; i < this->ncols(); ++i) {
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

template <class N>
constexpr usize ContactMatrix<N>::get_n_of_missed_updates() const noexcept {
  return static_cast<usize>(this->_updates_missed.load());
}

template <class N>
constexpr double ContactMatrix<N>::unsafe_get_fraction_of_missed_updates() const noexcept {
  if (this->empty()) {
    return 0.0;
  }
  return static_cast<double>(this->get_n_of_missed_updates()) /
         static_cast<double>(this->get_tot_contacts() + this->get_n_of_missed_updates());
}

template <class N>
constexpr usize ContactMatrix<N>::get_tot_contacts() const noexcept(utils::ndebug_defined()) {
  return static_cast<usize>(this->_tot_contacts.load());
}

template <class N>
double ContactMatrix<N>::get_avg_contact_density() const noexcept(utils::ndebug_defined()) {
  return static_cast<double>(this->get_tot_contacts()) / static_cast<double>(this->npixels());
}

template <class N>
constexpr usize ContactMatrix<N>::get_matrix_size_in_bytes() const
    noexcept(utils::ndebug_defined()) {
  return this->npixels() * sizeof(N);
}

template <class N>
constexpr double ContactMatrix<N>::get_matrix_size_in_mb() const noexcept(utils::ndebug_defined()) {
  const double mb = 1.0e6;
  return static_cast<double>(this->get_matrix_size_in_bytes()) / mb;
}

template <class N>
void ContactMatrix<N>::unsafe_print(std::ostream &out_stream, bool full) const {
  if (full) {
    std::vector<N> row(this->_ncols, 0);
    for (usize y = 0; y < this->_ncols; ++y) {
      std::fill(row.begin(), row.end(), 0);
      for (usize x = 0; x < this->_ncols; ++x) {
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
      fmt::print(out_stream, FMT_STRING("{}\n"), fmt::join(row, "\t"));
    }
  } else {
    std::vector<N> row(this->ncols());
    for (usize i = 0; i < this->nrows(); ++i) {
      for (auto j = i; j < this->ncols(); ++j) {
        row[j] = this->at(i, j);
      }
      fmt::print(out_stream, FMT_STRING("{}\n"), fmt::join(row, "\t"));
    }
  }
  out_stream << std::flush;
}

template <class N>
void ContactMatrix<N>::unsafe_print(bool full) const {
  this->unsafe_print(std::cout, full);
}

template <class N>
std::vector<std::vector<N>> ContactMatrix<N>::unsafe_generate_symmetric_matrix() const {
  std::vector<std::vector<N>> m;
  m.reserve(this->_ncols);
  for (usize y = 0; y < this->_ncols; ++y) {
    std::vector<N> row(this->_ncols, 0);
    for (usize x = 0; x < this->_ncols; ++x) {
      const auto [j, i] = [&]() {
        if (y > x) {
          return std::make_pair(y, y - x);
        }
        return std::make_pair(x, x - y);
      }();

      if (i < this->_nrows) {
        row[x] = this->at(i, j);
      }
    }
    m.emplace_back(std::move(row));
  }
  return m;
}

template <class N>
void ContactMatrix<N>::unsafe_import_from_txt(const boost::filesystem::path &path, const char sep) {
  assert(boost::filesystem::exists(path));  // NOLINT
  std::ifstream fp(path.string());

  std::string buff;
  std::vector<std::string_view> toks;
  usize i = 0;
  while (std::getline(fp, buff)) {
    toks = absl::StrSplit(buff, sep);
    if (i == 0) {
      this->unsafe_resize(toks.size(), toks.size());
      this->unsafe_reset();
    }
    for (usize j = i; j < this->ncols(); ++j) {
      this->set(i, j, utils::parse_numeric_or_throw<N>(toks[j]));
    }
    ++i;
  }
  assert(i != 0);  // NOLINT guard against empty files
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_ATTRIBUTES
template <class N>
__attribute__((no_sanitize("integer"))) void
ContactMatrix<N>::unsafe_generate_mask_for_bins_without_contacts(
    boost::dynamic_bitset<> &mask) const {
  DISABLE_WARNING_POP
  mask.resize(this->ncols());
  mask.reset();

  for (usize i = 0; i < this->ncols(); ++i) {
    // Set bitmask to 1 if row contains at least one non-zero value
    for (auto j = i; j < (i + this->nrows()) && j < this->ncols(); ++j) {
      if ((mask[i] |= this->unsafe_get(i, j))) {
        break;
      }
    }
    // Set bitmask to 1 if column "above" the current bin contains at least one non-zero value
    // NOTE i - this->nrows() can overflow. This is intended
    for (auto j = i; j > 0 && j > (i - this->nrows()); --j) {
      if ((mask[i] |= this->unsafe_get(i, j))) {
        break;
      }
    }
  }
}

template <class N>
boost::dynamic_bitset<> ContactMatrix<N>::unsafe_generate_mask_for_bins_without_contacts() const {
  boost::dynamic_bitset<> mask{};
  this->unsafe_generate_mask_for_bins_without_contacts(mask);
  return mask;
}

template <class N>
void ContactMatrix<N>::clear_missed_updates_counter() {
  this->_updates_missed = 0;
}

template <class N>
void ContactMatrix<N>::unsafe_reset() {
  std::fill(this->_contacts.begin(), this->_contacts.end(), 0);
  this->_tot_contacts = 0;
  this->clear_missed_updates_counter();
}

template <class N>
void ContactMatrix<N>::unsafe_resize(const usize nrows, const usize ncols) {
  if (nrows == this->_nrows && ncols == this->_ncols) {
    return;
  }

  if (ncols != this->_ncols) {
    std::vector<std::mutex> mtxes(ncols);
    std::swap(this->_mtxes, mtxes);
  }

  this->_nrows = std::min(nrows, ncols);
  this->_ncols = ncols;
  this->_contacts.resize(this->npixels(), N(0));
}

template <class N>
void ContactMatrix<N>::unsafe_resize(const bp_t length, const bp_t diagonal_width,
                                     const bp_t bin_size) {
  const auto nrows = (diagonal_width + bin_size - 1) / bin_size;
  const auto ncols = (length + bin_size - 1) / bin_size;
  this->unsafe_resize(nrows, ncols);
}

template <class N>
bool ContactMatrix<N>::empty() const {
  return this->_tot_contacts == 0;
}

template <class N>
absl::Span<const N> ContactMatrix<N>::get_raw_count_vector() const {
  return absl::MakeConstSpan(this->_contacts);
}

template <class N>
absl::Span<N> ContactMatrix<N>::get_raw_count_vector() {
  return absl::MakeSpan(this->_contacts);
}

template <class N>
void ContactMatrix<N>::unsafe_compute_row_wise_contact_histogram(std::vector<u64> &buff) const {
  buff.resize(this->nrows());
  std::fill(buff.begin(), buff.end(), 0);

  for (usize i = 0; i < this->ncols(); ++i) {
    for (auto j = i; j < i + this->nrows() && j < this->ncols(); ++j) {
      // j - i corresponds to the distance from the diagonal
      buff[j - i] += this->unsafe_get(j, i);
    }
  }
}

template <class N>
std::vector<u64> ContactMatrix<N>::unsafe_compute_row_wise_contact_histogram() const {
  std::vector<u64> buff(this->nrows());
  this->unsafe_compute_row_wise_contact_histogram(buff);
  return buff;
}

template <class N>
void ContactMatrix<N>::unsafe_deplete_contacts(double depletion_multiplier) {
  const auto hist = this->unsafe_compute_row_wise_contact_histogram();
  const auto effective_nbins = this->unsafe_generate_mask_for_bins_without_contacts().count();
  // This histogram contains the average contact number (instead of the total)
  std::vector<u64> row_wise_avg_contacts(hist.size());
  std::transform(hist.begin(), hist.end(), row_wise_avg_contacts.begin(), [&](const auto n) {
    return static_cast<u64>(std::round((depletion_multiplier * static_cast<double>(n)) /
                                       static_cast<double>(effective_nbins)));
  });

  for (usize i = 0; i < this->ncols(); ++i) {
    for (usize j = i; j < i + this->nrows() && j < this->ncols(); ++j) {
      // j - i corresponds to the distance from the diagonal
      if (this->unsafe_get(j, i) > row_wise_avg_contacts[j - i]) {
        this->subtract(j, i, static_cast<N>(row_wise_avg_contacts[j - i]));
      } else {
        this->set(j, i, N(0));
      }
    }
  }
}

template <class N>
N &ContactMatrix<N>::at(const usize i, const usize j) noexcept(utils::ndebug_defined()) {
  if constexpr (utils::ndebug_not_defined()) {
    if ((j * this->_nrows) + i > this->_contacts.size()) {
      throw std::runtime_error(fmt::format(FMT_STRING("ContactMatrix::at tried to access element "
                                                      "m[{}][{}] of a matrix of shape [{}][{}]! "
                                                      "({} >= {})"),
                                           i, j, this->nrows(), this->ncols(),
                                           (j * this->_nrows) + i, this->_contacts.size()));
    }
  }
  return this->_contacts[(j * this->_nrows) + i];
}

template <class N>
const N &ContactMatrix<N>::at(const usize i, const usize j) const
    noexcept(utils::ndebug_defined()) {
  if constexpr (utils::ndebug_not_defined()) {
    if ((j * this->_nrows) + i > this->_contacts.size()) {
      throw std::runtime_error(fmt::format(FMT_STRING("ContactMatrix::at tried to access element "
                                                      "m[{}][{}] of a matrix of shape [{}][{}]! "
                                                      "({} >= {})"),
                                           i, j, this->nrows(), this->ncols(),
                                           (j * this->_nrows) + i, this->_contacts.size()));
    }
  }
  return this->_contacts[(j * this->_nrows) + i];
}

template <class N>
std::pair<usize, usize> ContactMatrix<N>::transpose_coords(
    const usize row, const usize col) noexcept(utils::ndebug_defined()) {
  if (row > col) {
    return std::make_pair(row - col, row);
  }
  return std::make_pair(col - row, col);
}

template <class N>
void ContactMatrix<N>::bound_check_column(const usize col) const {
  if (col > this->ncols()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("caught an attempt to access element past the end of the contact matrix: "
                   "j={}; ncols={}: j > ncols"),
        col, this->ncols()));
  }
}

template <class N>
void ContactMatrix<N>::check_for_overflow_on_add(const usize row, const usize col,
                                                 const N n) const {
  assert(n >= 0);
  const auto lo = (std::numeric_limits<N>::min)();
  const auto hi = (std::numeric_limits<N>::max)();

  const auto m = this->at(row, col);
  if (hi - n < m) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Overflow detected: incrementing m={} by n={} would result in a "
                               "number outside of range {}-{}"),
                    m, n, lo, hi));
  }
  if (hi - n < this->_tot_contacts) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Overflow detected: incrementing _tot_contacts={} by n={} would result in a "
                   "number outside of range {}-{}"),
        this->_tot_contacts, n, lo, hi));
  }
}

template <class N>
void ContactMatrix<N>::check_for_overflow_on_subtract(usize row, usize col, const N n) const {
  assert(n >= 0);
  const auto lo = (std::numeric_limits<N>::min)();
  const auto hi = (std::numeric_limits<N>::max)();

  const auto m = this->at(row, col);
  if (lo + n > m) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Overflow detected: decrementing m={} by n={} would result in a "
                               "number outside of range {}-{}"),
                    m, n, lo, hi));
  }
  if (lo + n > this->_tot_contacts) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Overflow detected: decrementing _tot_contacts={} by n={} would result in a "
                   "number outside of range {}-{}"),
        this->_tot_contacts, n, lo, hi));
  }
}

}  // namespace modle

// IWYU pragma: private, include "modle/contacts.hpp"
// IWYU pragma: no_include <boost/core/checked_delete.hpp>
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <boost/exception/exception.hpp>
// IWYU pragma: no_include <boost/move/utility_core.hpp>
// IWYU pragma: no_include <H5Public.h>
