// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>
#include <fmt/format.h>
#include <xxhash.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <mutex>
#include <shared_mutex>
#include <stdexcept>
#include <utility>

#include "modle/common/common.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"
#include "modle/common/utils.hpp"
#include "modle/internal/contact_matrix_internal.hpp"

namespace modle {

template <class N>
ContactMatrixDense<N>::ContactMatrixDense(const ContactMatrixDense<N> &other)
    : _nrows(other.nrows()),
      _ncols(other.ncols()),
      _contacts(other._contacts),
      _mtxes(compute_number_of_mutexes(nrows(), ncols())),
      _tot_contacts(other._tot_contacts.load()),
      _nnz(other._nnz.load()),
      _global_stats_outdated(other._global_stats_outdated.load()),
      _updates_missed(other._updates_missed.load()) {}

template <class N>
ContactMatrixDense<N>::ContactMatrixDense(const std::size_t nrows_, const std::size_t ncols_)
    : _nrows(std::min(nrows_, ncols_)),
      _ncols(ncols_),
      _contacts(_nrows * _ncols + 1, N(0)),
      _mtxes(compute_number_of_mutexes(nrows(), ncols())) {}

template <class N>
ContactMatrixDense<N>::ContactMatrixDense(const bp_t length, const bp_t diagonal_width,
                                          const bp_t bin_size)
    : ContactMatrixDense((diagonal_width + bin_size - 1) / bin_size,
                         (length + bin_size - 1) / bin_size) {}

template <class N>
ContactMatrixDense<N>::ContactMatrixDense(const absl::Span<const N> contacts,
                                          const std::size_t nrows_, const std::size_t ncols_,
                                          const std::size_t tot_contacts,
                                          const std::size_t updates_missed)
    : _nrows(nrows_),
      _ncols(ncols_),
      _contacts(contacts.begin(), contacts.end()),
      _mtxes(compute_number_of_mutexes(nrows(), ncols())),
      _tot_contacts(static_cast<std::int64_t>(tot_contacts)),
      _global_stats_outdated(true),
      _updates_missed(updates_missed) {
  assert(_contacts.size() == _nrows * _ncols + 1);
  if (tot_contacts == 0) {
    unsafe_update_global_stats();
  }
}

template <class N>
ContactMatrixDense<N> &ContactMatrixDense<N>::operator=(const ContactMatrixDense<N> &other) {
  if (this == &other) {
    return *this;
  }

  const auto lck2 = other.lock();
  _nrows = other.nrows();
  _ncols = other.ncols();
  if (!_contacts.empty()) {
    _contacts.resize(other._contacts.size());
    std::copy(other._contacts.begin(), other._contacts.end(), _contacts.begin());
  } else {
    _contacts = other._contacts;
  }
  _mtxes = std::vector<mutex_t>(other._mtxes.size());
  _tot_contacts = other._tot_contacts.load();
  _nnz = other._nnz.load();
  _global_stats_outdated = other._global_stats_outdated.load();
  _updates_missed = other._updates_missed.load();

  return *this;
}

template <class N>
std::unique_lock<typename ContactMatrixDense<N>::mutex_t> ContactMatrixDense<N>::lock_pixel(
    std::size_t row, std::size_t col) const {
  return std::unique_lock<ContactMatrixDense<N>::mutex_t>(_mtxes[get_pixel_mutex_idx(row, col)]);
}

template <class N>

utils::LockRangeExclusive<typename ContactMatrixDense<N>::mutex_t> ContactMatrixDense<N>::lock()
    const {
  return utils::LockRangeExclusive<mutex_t>(_mtxes);
}

template <class N>
constexpr std::size_t ContactMatrixDense<N>::ncols() const {
  return _ncols;
}

template <class N>
constexpr std::size_t ContactMatrixDense<N>::nrows() const {
  return _nrows;
}

template <class N>
constexpr std::size_t ContactMatrixDense<N>::npixels() const {
  return _nrows * _ncols;
}

template <class N>
constexpr std::size_t ContactMatrixDense<N>::get_n_of_missed_updates() const noexcept {
  return _updates_missed.load();
}

template <class N>
constexpr std::size_t ContactMatrixDense<N>::get_matrix_size_in_bytes() const {
  return ((npixels() + 1) * sizeof(N)) + (_mtxes.size() * sizeof(mutex_t));
}

template <class N>
void ContactMatrixDense<N>::clear_missed_updates_counter() {
  _updates_missed = 0;
}

template <class N>
absl::Span<const N> ContactMatrixDense<N>::get_raw_count_vector() const {
  return absl::MakeConstSpan(_contacts);
}

template <class N>
absl::Span<N> ContactMatrixDense<N>::get_raw_count_vector() {
  return absl::MakeSpan(_contacts);
}

template <class N>
void ContactMatrixDense<N>::bound_check_coords([[maybe_unused]] const std::size_t row,
                                               [[maybe_unused]] const std::size_t col) const {
  if constexpr (utils::ndebug_not_defined()) {
    internal::bound_check_coords(*this, row, col);
  }
}

template <class N>
std::size_t ContactMatrixDense<N>::hash_coordinates(const std::size_t i,
                                                    const std::size_t j) noexcept {
  const std::array<std::size_t, 2> buff{i, j};
  return utils::conditional_static_cast<std::size_t>(
      XXH3_64bits(buff.data(), sizeof(std::size_t) * 2));
}

template <class N>
constexpr std::size_t ContactMatrixDense<N>::compute_number_of_mutexes(
    const std::size_t rows, const std::size_t cols) noexcept {
  if (rows + cols == 0) {
    return 0;
  }
  const auto nthreads = static_cast<std::size_t>(std::thread::hardware_concurrency());
  const auto max_dim = std::max(rows, cols);
  // Clamping to 2-n is needed because get_pixel_mutex_idx expects the number of
  // mutexes to be a multiple of 2
  return utils::next_pow2(std::clamp(max_dim, std::size_t(2), 1000 * nthreads));
}

template <class N>
std::size_t ContactMatrixDense<N>::get_pixel_mutex_idx(const std::size_t row,
                                                       const std::size_t col) const noexcept {
  assert(!_mtxes.empty());
  assert(_mtxes.size() % 2 == 0);
  // equivalent to hash_coordinates(row, col) % _mtxes.size() when _mtxes.size() % 2 == 0
  return (hash_coordinates(row, col) & (_mtxes.size() - 1));
}

}  // namespace modle

// IWYU pragma: private, include "modle/contact_matrix_dense.hpp"
// IWYU pragma: no_include <boost/core/checked_delete.hpp>
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <boost/exception/exception.hpp>
// IWYU pragma: no_include <boost/move/utility_core.hpp>
// IWYU pragma: no_include <H5Public.h>
