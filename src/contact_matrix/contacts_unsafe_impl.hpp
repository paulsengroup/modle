// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/strings/str_split.h>  // for StrSplit, Splitter
#include <fmt/format.h>              // for FMT_STRING, join

#include <algorithm>                                // for clamp, fill, max, min
#include <atomic>                                   // for atomic_fetch_add_explicit
#include <boost/config.hpp>                         // IWYU pragma: keep for BOOST_UNLIKELY
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/filesystem/operations.hpp>          // for exists
#include <boost/filesystem/path.hpp>                // for path
#include <cassert>                                  // for assert
#include <cmath>                                    // for sqrt, round
#include <fstream>                                  // for flush, ostream
#include <iostream>                                 // for cout
#include <numeric>                                  // for accumulate
#include <shared_mutex>                             // for shared_mutex
#include <string>                                   // for allocator, string
#include <string_view>                              // for string_view
#include <utility>                                  // for make_pair
#include <vector>                                   // for vector

#include "modle/common/common.hpp"                // for usize, i64, u64, bp_t, isize
#include "modle/common/utils.hpp"                 // for convolve, parse_numeric_or_throw
#include "modle/compressed_io/compressed_io.hpp"  // for Reader
#include "modle/interval_tree.hpp"                // for IITree

namespace modle {

template <class N>
N ContactMatrix<N>::unsafe_get(const usize row, const usize col) const {
  const auto [i, j] = transpose_coords(row, col);
  this->bound_check_coords(i, j);

  if (i >= this->nrows()) {
    return 0;
  }

  return this->unsafe_at(i, j);
}

// NOTE the pixels returned by this function go from the diagonal towards the periphery
// Example: given the following matrix
//          1  2  3
//          2  4  5
//          3  5  6
// Fetching col #3 would yield 6 5 3
template <class N>
void ContactMatrix<N>::unsafe_get_column(const usize col, std::vector<N> &buff,
                                         const usize row_offset) const {
  assert(row_offset <= col);  // NOLINT
  const auto [rowt, colt] = transpose_coords(col - row_offset, col);
  this->bound_check_coords(rowt, colt);

  const auto first_idx = (colt * this->nrows()) + row_offset;
  const auto last_idx =
      first_idx + std::min(this->ncols() - colt - row_offset, this->nrows() - row_offset);

  buff.resize(last_idx - first_idx);
  assert(buff.size() <= this->nrows());  // NOLINT

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
                                      const usize col_offset) const {
  assert(row >= col_offset);  // NOLINT
  buff.resize(std::clamp(this->ncols() - row - col_offset, usize(0), this->nrows() - col_offset));

  for (usize i = 0; i < buff.size(); ++i) {
    buff[i] = this->unsafe_get(row, row + col_offset + i);
  }
}

template <class N>
N ContactMatrix<N>::unsafe_get_block(const usize row, const usize col,
                                     const usize block_size) const {
  assert(block_size > 0);              // NOLINT
  assert(block_size < this->nrows());  // NOLINT
  // For now we only support blocks with an odd size
  assert(block_size % 2 != 0);  // NOLINT
  if (block_size == 1) {
    return this->get(row, col);
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
void ContactMatrix<N>::unsafe_get_block(const usize row, const usize col, const usize block_size,
                                        std::vector<N> &buff) const {
  assert(block_size > 0);              // NOLINT
  assert(block_size < this->nrows());  // NOLINT
  // For now we only support blocks with an odd size
  assert(block_size % 2 != 0);  // NOLINT
  if (BOOST_UNLIKELY(block_size == 1)) {
    buff.resize(1);
    buff.front() = this->get(row, col);
    return;
  }

  // Edges are handled like shown here: https://en.wikipedia.org/wiki/File:Extend_Edge-Handling.png
  const auto bs = static_cast<i64>(block_size);
  const auto first_row = static_cast<i64>(row) - (bs / 2);
  const auto first_col = static_cast<i64>(col) - (bs / 2);
  buff.resize(block_size * block_size);

  usize k = 0;
  for (auto i = first_row; i < first_row + bs; ++i) {
    for (auto j = first_col; j < first_col + bs; ++j) {
      const auto ii = static_cast<usize>(std::clamp(i, i64(0), static_cast<i64>(this->_ncols - 1)));
      const auto jj = static_cast<usize>(std::clamp(j, i64(0), static_cast<i64>(this->_ncols - 1)));
      buff[k++] = this->unsafe_get(ii, jj);
    }
  }
}

template <class N>
void ContactMatrix<N>::unsafe_set(const usize row, const usize col, const N n) {
  const auto [i, j] = this->transpose_coords(row, col);
  this->bound_check_coords(i, j);

  if (i > this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, i64(1), std::memory_order_relaxed);
    return;
  }

  this->unsafe_at(i, j) = n;
  this->_tot_contacts_outdated = true;
}

template <class N>
void ContactMatrix<N>::unsafe_add(const usize row, const usize col, const N n) {
  assert(n > 0);  // NOLINT Use subtract to add a negative number
  const auto [i, j] = transpose_coords(row, col);
  this->bound_check_coords(i, j);

  if (i > this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, i64(1), std::memory_order_relaxed);
    return;
  }

  this->unsafe_at(i, j) += n;
  this->_tot_contacts_outdated = true;
}

template <class N>
void ContactMatrix<N>::unsafe_subtract(const usize row, const usize col, const N n) {
  assert(n >= 0);  // NOLINT Use add to subtract a negative number
  const auto [i, j] = transpose_coords(row, col);
  this->bound_check_coords(i, j);

  if (i > this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, i64(1), std::memory_order_relaxed);
    return;
  }

  this->unsafe_at(i, j) -= n;
  this->_tot_contacts_outdated = true;
}

template <class N>
void ContactMatrix<N>::unsafe_increment(usize row, usize col) {
  this->unsafe_add(row, col, N(1));
}

template <class N>
void ContactMatrix<N>::unsafe_decrement(usize row, usize col) {
  this->unsafe_subtract(row, col, N(1));
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
constexpr double ContactMatrix<N>::unsafe_get_fraction_of_missed_updates() const noexcept {
  if (this->empty() || this->get_n_of_missed_updates() == 0) {
    return 0.0;
  }
  const auto missed_updates = static_cast<double>(this->get_n_of_missed_updates());
  return missed_updates / (static_cast<double>(this->unsafe_get_tot_contacts()) + missed_updates);
}

template <class N>
usize ContactMatrix<N>::unsafe_get_tot_contacts() const noexcept {
  if (this->_tot_contacts_outdated) {
    this->_tot_contacts = std::accumulate(this->_contacts.begin(), this->_contacts.end(), N(0));
    this->_tot_contacts_outdated = false;
  }
  return static_cast<usize>(this->_tot_contacts.load());
}

template <class N>
double ContactMatrix<N>::unsafe_get_avg_contact_density() const {
  return static_cast<double>(this->unsafe_get_tot_contacts()) /
         static_cast<double>(this->npixels());
}

template <class N>
N ContactMatrix<N>::unsafe_get_min_count() const noexcept {
  if (this->unsafe_get_tot_contacts() == 0) {
    return 0;
  }
  return *std::min_element(this->_contacts.begin(), this->_contacts.end());
}

template <class N>
N ContactMatrix<N>::unsafe_get_max_count() const noexcept {
  if (this->unsafe_get_tot_contacts() == 0) {
    return 0;
  }
  return *std::max_element(this->_contacts.begin(), this->_contacts.end());
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
        row[x] = this->unsafe_at(i, j);
      }
    }
    m.emplace_back(std::move(row));
  }
  return m;
}

template <class N>
void ContactMatrix<N>::unsafe_import_from_txt(const boost::filesystem::path &path, const char sep) {
  assert(boost::filesystem::exists(path));  // NOLINT
  compressed_io::Reader r(path);

  std::string buff;
  std::vector<std::string_view> toks;
  usize i = 0;
  while (r.getline(buff)) {
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

template <class N>
void ContactMatrix<N>::unsafe_generate_mask_for_bins_without_contacts(
    boost::dynamic_bitset<> &mask) const {
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
    std::vector<std::shared_mutex> mtxes(compute_number_of_mutexes(nrows, ncols));
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
N &ContactMatrix<N>::unsafe_at(const usize i, const usize j) {
  // this->bound_check_coords(i, j);
  return this->_contacts[(j * this->_nrows) + i];
}

template <class N>
const N &ContactMatrix<N>::unsafe_at(const usize i, const usize j) const {
  // this->bound_check_coords(i, j);
  return this->_contacts[(j * this->_nrows) + i];
}

template <class N>
ContactMatrix<double> ContactMatrix<N>::unsafe_gaussian_diff(const double sigma1,
                                                             const double sigma2,
                                                             const double min_value,
                                                             const double max_value) const {
  assert(sigma1 <= sigma2);  // NOLINT
  ContactMatrix<double> bmatrix(this->nrows(), this->ncols());

  const auto gauss_kernel1 = stats::compute_gauss_kernel(sigma1);
  const auto gauss_kernel2 = stats::compute_gauss_kernel(/*sqrt(kernel1.size(), */ sigma2);
  std::vector<N> pixels(std::max(gauss_kernel1.size(), gauss_kernel2.size()));

  [[maybe_unused]] const auto block_size1 =
      static_cast<usize>(std::sqrt(static_cast<double>(gauss_kernel1.size())));
  [[maybe_unused]] const auto block_size2 =
      static_cast<usize>(std::sqrt(static_cast<double>(gauss_kernel2.size())));
  assert(block_size1 * block_size1 <= gauss_kernel1.size());  // NOLINT
  assert(block_size2 * block_size2 <= gauss_kernel2.size());  // NOLINT

  for (usize i = 0; i < this->ncols(); ++i) {
    for (usize j = i; j < std::min(i + this->nrows(), this->ncols()); ++j) {
      this->unsafe_get_block(i, j, block_size1, pixels);
      const auto n1 = utils::convolve(gauss_kernel1, pixels);
      this->unsafe_get_block(i, j, block_size2, pixels);
      const auto n2 = utils::convolve(gauss_kernel2, pixels);

      bmatrix.set(i, j, std::clamp(n1 - n2, min_value, max_value));
    }
  }
  return bmatrix;
}

template <class N>
template <class M, class>
void ContactMatrix<N>::unsafe_normalize(const ContactMatrix<N> &input_matrix,
                                        ContactMatrix<M> &output_matrix, M lb, M ub) noexcept {
  assert(ub >= lb);  // NOLINT
  // The unsafe resize takes care of the case where &input_matrix == &output_matrix
  output_matrix.unsafe_resize(input_matrix.nrows(), input_matrix.ncols());

  const auto min_count = input_matrix.unsafe_get_min_count();
  const auto max_count = input_matrix.unsafe_get_max_count();
  const auto scaling_factor = ub - lb;

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  std::transform(input_matrix._contacts.begin(), input_matrix._contacts.end(),
                 output_matrix._contacts.begin(), [&](const auto count) {
                   // https://stats.stackexchange.com/a/281164
                   const auto n =
                       static_cast<M>(count - min_count) / static_cast<M>(max_count - min_count);
                   return (n * scaling_factor) + static_cast<M>(lb);
                 });
  DISABLE_WARNING_POP

  output_matrix._updates_missed = input_matrix._updates_missed.load();
  output_matrix._tot_contacts_outdated = true;
}

template <class N>
template <class FP, class>
ContactMatrix<FP> ContactMatrix<N>::unsafe_normalize(const double lb, const double ub) const {
  ContactMatrix<FP> m(this->nrows(), this->ncols());
  ContactMatrix<N>::unsafe_normalize(*this, m, lb, ub);
  return m;
}

template <class N>
inline void ContactMatrix<N>::unsafe_normalize_inplace(const N lb, const N ub) noexcept {
  ContactMatrix<N>::unsafe_normalize(*this, *this, lb, ub);
}

template <class N>
void ContactMatrix<N>::unsafe_clamp(const ContactMatrix<N> &input_matrix,
                                    ContactMatrix<N> &output_matrix, const N lb,
                                    const N ub) noexcept {
  assert(lb <= ub);  // NOLINT
  // The unsafe resize takes care of the case where &input_matrix == &output_matrix
  output_matrix.unsafe_resize(input_matrix.nrows(), input_matrix.ncols());

  std::transform(input_matrix._contacts.begin(), input_matrix._contacts.end(),
                 output_matrix._contacts.begin(),
                 [&](const auto count) { return std::clamp(count, lb, ub); });
  output_matrix._updates_missed = input_matrix._updates_missed.load();
  output_matrix._tot_contacts_outdated = true;
}

template <class N>
ContactMatrix<N> ContactMatrix<N>::unsafe_clamp(const N lb, const N ub) const {
  ContactMatrix<N> m(this->nrows(), this->ncols());
  ContactMatrix<N>::unsafe_clamp(*this, m, lb, ub);
  return m;
}

template <class N>
inline void ContactMatrix<N>::unsafe_clamp_inplace(const N lb, const N ub) noexcept {
  ContactMatrix<N>::unsafe_clamp(*this, *this, lb, ub);
}

template <class N>
template <class N1, class N2>
void ContactMatrix<N>::unsafe_discretize(const ContactMatrix<N> &input_matrix,
                                         ContactMatrix<N1> &output_matrix,
                                         const IITree<N2, N1> &mappings) noexcept {
  output_matrix.unsafe_resize(input_matrix.nrows(), input_matrix.ncols());
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  std::transform(input_matrix._contacts.begin(), input_matrix._contacts.end(),
                 output_matrix._contacts.begin(), [&](const auto n) {
                   if (auto it = mappings.find_overlaps(static_cast<N>(n), static_cast<N>(n));
                       it.first != it.second) {
                     return *it.first;
                   }
                   return static_cast<N1>(n);
                 });
  DISABLE_WARNING_POP
  output_matrix._tot_contacts_outdated = true;
  output_matrix._updates_missed = input_matrix._updates_missed.load();
}

template <class N>
template <class N1, class N2>
ContactMatrix<N1> ContactMatrix<N>::unsafe_discretize(const IITree<N2, N1> &mappings) const {
  ContactMatrix<N1> m(this->nrows(), this->ncols());
  ContactMatrix<N>::unsafe_discretize(*this, m, mappings);
  return m;
}

template <class N>
template <class M>
void ContactMatrix<N>::unsafe_discretize_inplace(const IITree<M, N> &mappings) noexcept {
  ContactMatrix<N>::unsafe_discretize(*this, *this, mappings);
}

template <class N>
template <class M, class>
ContactMatrix<M> ContactMatrix<N>::unsafe_as() const {
  ContactMatrix<M> m(this->nrows(), this->ncols());
  std::transform(this->_contacts.begin(), this->_contacts.end(), m._contacts.begin(),
                 [](const auto n) {
                   if constexpr (std::is_floating_point_v<N> && !std::is_floating_point_v<M>) {
                     return static_cast<M>(std::round(n));
                   } else {
                     return static_cast<M>(n);
                   }
                 });
  m._tot_contats_outdated = true;
  m._updates_missed = this->_updates_missed.load();
}

}  // namespace modle

// IWYU pragma: private, include "modle/contacts.hpp"
// IWYU pragma: no_include <boost/core/checked_delete.hpp>
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <boost/exception/exception.hpp>
// IWYU pragma: no_include <boost/move/utility_core.hpp>
// IWYU pragma: no_include <H5Public.h>
