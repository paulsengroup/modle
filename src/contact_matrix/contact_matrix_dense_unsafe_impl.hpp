// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <atomic>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <shared_mutex>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/compressed_io/compressed_io.hpp"
#include "modle/interval_tree.hpp"

namespace modle {

template <class N>
N ContactMatrixDense<N>::unsafe_get(const std::size_t row, const std::size_t col) const {
  const auto [i, j] = internal::transpose_coords(row, col);
  bound_check_coords(i, j);

  if (i >= nrows()) {
    return 0;
  }

  return unsafe_at(i, j);
}

// NOTE the pixels returned by this function go from the diagonal towards the periphery
// Example: given the following matrix
//          1  2  3
//          2  4  5
//          3  5  6
// Fetching col #3 would yield 6 5 3
template <class N>
void ContactMatrixDense<N>::unsafe_get_column(const std::size_t col, std::vector<N> &buff,
                                              const std::size_t row_offset) const {
  assert(row_offset <= col);
  const auto [rowt, colt] = internal::transpose_coords(col - row_offset, col);
  bound_check_coords(rowt, colt);

  const auto first_idx = (colt * nrows()) + row_offset;
  const auto last_idx = first_idx + std::min(ncols() - colt - row_offset, nrows() - row_offset);

  buff.resize(last_idx - first_idx);
  assert(buff.size() <= nrows());

  const auto first = _contacts.begin() + static_cast<std::ptrdiff_t>(first_idx);
  const auto last = _contacts.begin() + static_cast<std::ptrdiff_t>(last_idx);

  std::copy(first, last, buff.begin());
}

// NOTE the pixels returned by this function go from the diagonal towards the periphery
// Example: given the following matrix
//          1  2  3
//          2  4  5
//          3  5  6
// Fetching col #1 would yield 1 2 3
template <class N>
void ContactMatrixDense<N>::unsafe_get_row(const std::size_t row, std::vector<N> &buff,
                                           const std::size_t col_offset) const {
  assert(row >= col_offset);
  buff.resize(std::clamp(ncols() - row - col_offset, std::size_t(0), nrows() - col_offset));

  for (std::size_t i = 0; i < buff.size(); ++i) {
    buff[i] = unsafe_get(row, row + col_offset + i);
  }
}

template <class N>
N ContactMatrixDense<N>::unsafe_get_block(const std::size_t row, const std::size_t col,
                                          const std::size_t block_size) const {
  assert(block_size > 0);
  assert(block_size < nrows());
  // For now we only support blocks with an odd size
  assert(block_size % 2 != 0);
  if (block_size == 1) {
    return unsafe_get(row, col);
  }

  // Edges are handled like shown here: https://en.wikipedia.org/wiki/File:Extend_Edge-Handling.png
  const auto bs = static_cast<std::int64_t>(block_size);
  const auto first_row = static_cast<std::int64_t>(row) - (bs / 2);
  const auto first_col = static_cast<std::int64_t>(col) - (bs / 2);

  N n{0};
  for (auto i = first_row; i < first_row + bs; ++i) {
    for (auto j = first_col; j < first_col + bs; ++j) {
      const auto ii = static_cast<std::size_t>(
          std::clamp(i, std::int64_t(0), static_cast<std::int64_t>(_ncols - 1)));
      const auto jj = static_cast<std::size_t>(
          std::clamp(j, std::int64_t(0), static_cast<std::int64_t>(_ncols - 1)));
      n += unsafe_get(ii, jj);
    }
  }
  return n;
}

template <class N>
void ContactMatrixDense<N>::unsafe_get_block(const std::size_t row, const std::size_t col,
                                             const std::size_t block_size,
                                             std::vector<N> &buff) const {
  assert(block_size > 0);
  assert(block_size < nrows());
  // For now we only support blocks with an odd size
  assert(block_size % 2 != 0);
  if (block_size == 1) [[unlikely]] {
    buff.resize(1);
    buff.front() = unsafe_get(row, col);
    return;
  }

  // Edges are handled like shown here: https://en.wikipedia.org/wiki/File:Extend_Edge-Handling.png
  const auto bs = static_cast<std::int64_t>(block_size);
  const auto first_row = static_cast<std::int64_t>(row) - (bs / 2);
  const auto first_col = static_cast<std::int64_t>(col) - (bs / 2);
  buff.resize(block_size * block_size);

  std::size_t k = 0;
  for (auto i = first_row; i < first_row + bs; ++i) {
    for (auto j = first_col; j < first_col + bs; ++j) {
      const auto ii = static_cast<std::size_t>(
          std::clamp(i, std::int64_t(0), static_cast<std::int64_t>(_ncols - 1)));
      const auto jj = static_cast<std::size_t>(
          std::clamp(j, std::int64_t(0), static_cast<std::int64_t>(_ncols - 1)));
      buff[k++] = unsafe_get(ii, jj);
    }
  }
}

template <class N>
void ContactMatrixDense<N>::unsafe_set(const std::size_t row, const std::size_t col, const N n) {
  const auto [i, j] = internal::transpose_coords(row, col);
  bound_check_coords(i, j);

  if (i >= nrows()) {
    std::atomic_fetch_add_explicit(&_updates_missed, std::size_t(1), std::memory_order_relaxed);
    return;
  }

  unsafe_at(i, j) = n;
  _global_stats_outdated = true;
}

template <class N>
void ContactMatrixDense<N>::unsafe_add(const std::size_t row, const std::size_t col, const N n) {
  assert(n > 0);
  const auto [i, j] = internal::transpose_coords(row, col);
  bound_check_coords(i, j);

  if (i >= nrows()) {
    std::atomic_fetch_add_explicit(&_updates_missed, std::size_t(1), std::memory_order_relaxed);
    return;
  }

  unsafe_at(i, j) += n;
  _global_stats_outdated = true;
}

template <class N>
void ContactMatrixDense<N>::unsafe_subtract(const std::size_t row, const std::size_t col,
                                            const N n) {
  assert(n >= 0);
  const auto [i, j] = internal::transpose_coords(row, col);
  bound_check_coords(i, j);

  if (i >= nrows()) {
    std::atomic_fetch_add_explicit(&_updates_missed, std::size_t(1), std::memory_order_relaxed);
    return;
  }

  unsafe_at(i, j) -= n;
  _global_stats_outdated = true;
}

template <class N>
void ContactMatrixDense<N>::unsafe_increment(std::size_t row, std::size_t col) {
  unsafe_add(row, col, N(1));
}

template <class N>
void ContactMatrixDense<N>::unsafe_decrement(std::size_t row, std::size_t col) {
  unsafe_subtract(row, col, N(1));
}

template <class N>
constexpr double ContactMatrixDense<N>::unsafe_get_fraction_of_missed_updates() const noexcept {
  if (empty() || get_n_of_missed_updates() == 0) {
    return 0.0;
  }
  const auto missed_updates = static_cast<double>(get_n_of_missed_updates());
  return missed_updates / (static_cast<double>(unsafe_get_tot_contacts()) + missed_updates);
}

template <class N>
auto ContactMatrixDense<N>::unsafe_get_tot_contacts() const noexcept -> SumT {
  if (_global_stats_outdated) {
    unsafe_update_global_stats();
  }
  return _tot_contacts.load();
}

template <class N>
std::size_t ContactMatrixDense<N>::unsafe_get_nnz() const noexcept {
  if (_global_stats_outdated) {
    unsafe_update_global_stats();
  }
  return _nnz.load();
}

template <class N>
double ContactMatrixDense<N>::unsafe_get_avg_contact_density() const {
  return static_cast<double>(unsafe_get_tot_contacts()) / static_cast<double>(npixels());
}

template <class N>
N ContactMatrixDense<N>::unsafe_get_min_count() const noexcept {
  if (unsafe_get_tot_contacts() == 0) {
    return 0;
  }
  return *std::min_element(_contacts.begin(), _contacts.end());
}

template <class N>
N ContactMatrixDense<N>::unsafe_get_max_count() const noexcept {
  if (unsafe_get_tot_contacts() == 0) {
    return 0;
  }
  return *std::max_element(_contacts.begin(), _contacts.end());
}

template <class N>
void ContactMatrixDense<N>::unsafe_print(std::ostream &out_stream, bool full) const {
  if (full) {
    std::vector<N> row(_ncols, 0);
    for (std::size_t y = 0; y < _ncols; ++y) {
      std::fill(row.begin(), row.end(), 0);
      for (std::size_t x = 0; x < _ncols; ++x) {
        auto j = x;
        auto i = j - y;
        if (y > x) {
          j = y;
          i = j - x;
        }
        if (i >= _nrows) {
          row[x] = 0;
        } else {
          row[x] = unsafe_at(i, j);
        }
      }
      fmt::print(out_stream, "{}\n", fmt::join(row, "\t"));
    }
  } else {
    std::vector<N> row(ncols());
    for (std::size_t i = 0; i < nrows(); ++i) {
      for (auto j = i; j < ncols(); ++j) {
        row[j] = unsafe_at(i, j);
      }
      fmt::print(out_stream, "{}\n", fmt::join(row, "\t"));
    }
  }
  out_stream << std::flush;
}

template <class N>
void ContactMatrixDense<N>::unsafe_print(bool full) const {
  unsafe_print(std::cout, full);
}

template <class N>
std::vector<std::vector<N>> ContactMatrixDense<N>::unsafe_generate_symmetric_matrix() const {
  std::vector<std::vector<N>> m;
  m.reserve(_ncols);
  for (std::size_t y = 0; y < _ncols; ++y) {
    std::vector<N> row(_ncols, 0);
    for (std::size_t x = 0; x < _ncols; ++x) {
      const auto [j, i] = [&]() {
        if (y > x) {
          return std::make_pair(y, y - x);
        }
        return std::make_pair(x, x - y);
      }();

      if (i < _nrows) {
        row[x] = unsafe_at(i, j);
      }
    }
    m.emplace_back(std::move(row));
  }
  return m;
}

template <class N>
void ContactMatrixDense<N>::unsafe_reset() {
  std::fill(_contacts.begin(), _contacts.end(), 0);
  _tot_contacts = 0;
  clear_missed_updates_counter();
}

template <class N>
void ContactMatrixDense<N>::unsafe_resize(const std::size_t nrows, const std::size_t ncols) {
  if (nrows == _nrows && ncols == _ncols) {
    return;
  }

  if (const auto num_mtxes = compute_number_of_mutexes(nrows, ncols); num_mtxes > _mtxes.size()) {
    std::vector<mutex_t> mtxes(num_mtxes);
    std::swap(_mtxes, mtxes);
  }

  _nrows = std::min(nrows, ncols);
  _ncols = ncols;
  if (npixels() >= _contacts.size()) {
    _contacts.resize(npixels() + 1, N(0));
  }
}

template <class N>
void ContactMatrixDense<N>::unsafe_resize(const bp_t length, const bp_t diagonal_width,
                                          const bp_t bin_size) {
  const auto nrows = (diagonal_width + bin_size - 1) / bin_size;
  const auto ncols = (length + bin_size - 1) / bin_size;
  unsafe_resize(nrows, ncols);
}

template <class N>
N &ContactMatrixDense<N>::unsafe_at(const std::size_t i, const std::size_t j) {
  return _contacts[internal::encode_idx(i, j, _nrows)];
}

template <class N>
const N &ContactMatrixDense<N>::unsafe_at(const std::size_t i, const std::size_t j) const {
  return _contacts[internal::encode_idx(i, j, _nrows)];
}

template <class N>
template <class M, class>
void ContactMatrixDense<N>::unsafe_normalize(const ContactMatrixDense<N> &input_matrix,
                                             ContactMatrixDense<M> &output_matrix, M lb,
                                             M ub) noexcept {
  assert(ub >= lb);
  // The unsafe resize takes care of the case where &input_matrix == &output_matrix
  output_matrix.unsafe_resize(input_matrix.nrows(), input_matrix.ncols());

  if (input_matrix.unsafe_empty()) {
    if (&input_matrix == &output_matrix) {
      return;
    }
    std::fill(output_matrix._contacts.begin(), output_matrix._contacts.end(), M(0));
    output_matrix._updates_missed = input_matrix._updates_missed.load();
    output_matrix._tot_contacts = 0;
    output_matrix._global_stats_outdated = false;
    return;
  }

  const auto min_count = input_matrix.unsafe_get_min_count();
  const auto max_count = input_matrix.unsafe_get_max_count();
  const auto scaling_factor = ub - lb;

  std::transform(input_matrix._contacts.begin(), input_matrix._contacts.end(),
                 output_matrix._contacts.begin(), [&](const auto count) {
                   // https://stats.stackexchange.com/a/281164
                   const auto n = utils::conditional_static_cast<M>(count - min_count) /
                                  utils::conditional_static_cast<M>(max_count - min_count);
                   return (n * scaling_factor) + utils::conditional_static_cast<M>(lb);
                 });

  output_matrix._updates_missed = input_matrix._updates_missed.load();
  output_matrix._global_stats_outdated = true;
}

template <class N>
template <class FP, class>
ContactMatrixDense<FP> ContactMatrixDense<N>::unsafe_normalize(const double lb,
                                                               const double ub) const {
  ContactMatrixDense<FP> m(nrows(), ncols());
  ContactMatrixDense<N>::unsafe_normalize(*this, m, lb, ub);
  return m;
}

template <class N>
inline void ContactMatrixDense<N>::unsafe_normalize_inplace(const N lb, const N ub) noexcept {
  ContactMatrixDense<N>::unsafe_normalize(*this, *this, lb, ub);
}

template <class N>
void ContactMatrixDense<N>::unsafe_clamp(const ContactMatrixDense<N> &input_matrix,
                                         ContactMatrixDense<N> &output_matrix, const N lb,
                                         const N ub) noexcept {
  assert(lb <= ub);
  // The unsafe resize takes care of the case where &input_matrix == &output_matrix
  output_matrix.unsafe_resize(input_matrix.nrows(), input_matrix.ncols());

  std::transform(input_matrix._contacts.begin(), input_matrix._contacts.end(),
                 output_matrix._contacts.begin(),
                 [&](const auto count) { return std::clamp(count, lb, ub); });
  output_matrix._updates_missed = input_matrix._updates_missed.load();
  output_matrix._global_stats_outdated = true;
}

template <class N>
ContactMatrixDense<N> ContactMatrixDense<N>::unsafe_clamp(const N lb, const N ub) const {
  ContactMatrixDense<N> m(nrows(), ncols());
  ContactMatrixDense<N>::unsafe_clamp(*this, m, lb, ub);
  return m;
}

template <class N>
inline void ContactMatrixDense<N>::unsafe_clamp_inplace(const N lb, const N ub) noexcept {
  ContactMatrixDense<N>::unsafe_clamp(*this, *this, lb, ub);
}

template <class N>
template <class N1, class N2>
void ContactMatrixDense<N>::unsafe_discretize(const ContactMatrixDense<N> &input_matrix,
                                              ContactMatrixDense<N1> &output_matrix,
                                              const IITree<N2, N1> &mappings) noexcept {
  output_matrix.unsafe_resize(input_matrix.nrows(), input_matrix.ncols());
  std::transform(input_matrix._contacts.begin(), input_matrix._contacts.end(),
                 output_matrix._contacts.begin(), [&](const auto n) {
                   if (auto it = mappings.find_overlaps(utils::conditional_static_cast<N>(n),
                                                        utils::conditional_static_cast<N>(n));
                       it.first != it.second) {
                     return *it.first;
                   }
                   return utils::conditional_static_cast<N1>(n);
                 });
  output_matrix._global_stats_outdated = true;
  output_matrix._updates_missed = input_matrix._updates_missed.load();
}

template <class N>
template <class N1, class N2>
ContactMatrixDense<N1> ContactMatrixDense<N>::unsafe_discretize(
    const IITree<N2, N1> &mappings) const {
  ContactMatrixDense<N1> m(nrows(), ncols());
  ContactMatrixDense<N>::unsafe_discretize(*this, m, mappings);
  return m;
}

template <class N>
template <class M>
void ContactMatrixDense<N>::unsafe_discretize_inplace(const IITree<M, N> &mappings) noexcept {
  ContactMatrixDense<N>::unsafe_discretize(*this, *this, mappings);
}

template <class N>
template <class M, class>
ContactMatrixDense<M> ContactMatrixDense<N>::unsafe_as() const {
  ContactMatrixDense<M> m(nrows(), ncols());
  std::transform(_contacts.begin(), _contacts.end(), m._contacts.begin(), [](const auto n) {
    if constexpr (std::is_floating_point_v<N> && !std::is_floating_point_v<M>) {
      return static_cast<M>(std::round(n));
    } else {
      return static_cast<M>(n);
    }
  });
  m._global_stats_outdated = true;
  m._updates_missed = _updates_missed.load();

  return m;
}

template <class N>
bool ContactMatrixDense<N>::unsafe_empty() const {
  return unsafe_get_tot_contacts() == 0;
}

template <class N>
void ContactMatrixDense<N>::unsafe_update_global_stats() const noexcept {
  assert(_global_stats_outdated);
  _nnz = std::size_t(
      std::count_if(_contacts.begin(), _contacts.end(), [&](const auto n) { return n != N(0); }));
  assert(_nnz <= npixels());

  _tot_contacts = std::accumulate(_contacts.begin(), _contacts.end(), SumT(0),
                                  [&](const auto accumulator, const auto n) {
                                    return accumulator + utils::conditional_static_cast<SumT>(n);
                                  });

  _global_stats_outdated = false;
}

}  // namespace modle

// IWYU pragma: private, include "modle/contact_matrix_dense.hpp"
// IWYU pragma: no_include <boost/core/checked_delete.hpp>
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <boost/exception/exception.hpp>
// IWYU pragma: no_include <boost/move/utility_core.hpp>
// IWYU pragma: no_include <H5Public.h>
