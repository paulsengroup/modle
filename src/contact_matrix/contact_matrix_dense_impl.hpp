// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span, MakeConstSpan, MakeSpan
#include <fmt/format.h>       // for FMT_STRING
#include <xxhash.h>           // for XXH_INLINE_XXH3_64bits, XXH3_64bits

#include <algorithm>  // for min, clamp
#include <array>      // for array
#include <cassert>    // for assert
#include <cmath>      // for sqrt
#include <fstream>
#include <limits>        // for numeric_limits
#include <mutex>         // for unique_lock
#include <shared_mutex>  // for shared_lock, shared_mutex
#include <stdexcept>     // for runtime_error, logic_error
#include <utility>       // for make_pair, pair

#include "modle/common/common.hpp"                      // for usize, u64, bp_t, i64
#include "modle/common/random.hpp"                      // for PRNG, uniform_int_distribution
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...
#include "modle/common/utils.hpp"                       // for ndebug_not_defined, next_pow2
#include "modle/internal/contact_matrix_internal.hpp"

namespace modle {

template <class N>
ContactMatrixDense<N>::ContactMatrixDense(const ContactMatrixDense<N> &other)
    : _nrows(other.nrows()),
      _ncols(other.ncols()),
      _contacts(other._contacts),
      _mtxes(compute_number_of_mutexes(this->nrows(), this->ncols())),
      _tot_contacts(other._tot_contacts.load()),
      _nnz(other._nnz.load()),
      _global_stats_outdated(other._global_stats_outdated.load()),
      _updates_missed(other._updates_missed.load()) {}

template <class N>
ContactMatrixDense<N>::ContactMatrixDense(const usize nrows, const usize ncols)
    : _nrows(std::min(nrows, ncols)),
      _ncols(ncols),
      _contacts(_nrows * _ncols + 1, N(0)),
      _mtxes(compute_number_of_mutexes(this->nrows(), this->ncols())) {}

template <class N>
ContactMatrixDense<N>::ContactMatrixDense(const bp_t length, const bp_t diagonal_width,
                                          const bp_t bin_size)
    : ContactMatrixDense((diagonal_width + bin_size - 1) / bin_size,
                         (length + bin_size - 1) / bin_size) {}

template <class N>
ContactMatrixDense<N>::ContactMatrixDense(const absl::Span<const N> contacts, const usize nrows,
                                          const usize ncols, const usize tot_contacts,
                                          const usize updates_missed)
    : _nrows(nrows),
      _ncols(ncols),
      _contacts(contacts.begin(), contacts.end()),
      _mtxes(compute_number_of_mutexes(this->nrows(), this->ncols())),
      _tot_contacts(static_cast<i64>(tot_contacts)),
      _global_stats_outdated(true),
      _updates_missed(static_cast<i64>(updates_missed)) {
  assert(_contacts.size() == _nrows * _ncols + 1);
  if (tot_contacts == 0) {
    this->unsafe_update_global_stats();
  }
}

template <class N>
ContactMatrixDense<N> ContactMatrixDense<N>::create_random_matrix(usize nrows, usize ncols,
                                                                  u64 seed) {
  ContactMatrixDense<N> m(nrows, ncols);
  auto rand_eng = random::PRNG(seed);

  auto dist = [&rand_eng]() {
    if constexpr (std::is_floating_point_v<N>) {
      return random::uniform_real_distribution<N>{0, 65553}(rand_eng);
    } else {
      u64 max_ = std::min(u64(65553), static_cast<u64>((std::numeric_limits<N>::max)()));
      return random::uniform_int_distribution<N>{0, static_cast<N>(max_)}(rand_eng);
    }
  };

  for (usize i = 0; i < m._ncols; ++i) {
    for (usize j = i; j < i + m._nrows && j < m._ncols; ++j) {
      m.set(i, j, dist());
    }
  }
  return m;
}

template <class N>
ContactMatrixDense<N> ContactMatrixDense<N>::create_random_matrix(bp_t length, bp_t diagonal_width,
                                                                  bp_t bin_size, u64 seed) {
  return ContactMatrixDense<N>::create_random_matrix((diagonal_width + bin_size - 1) / bin_size,
                                                     (length + bin_size - 1) / bin_size, seed);
}


template <class N>
ContactMatrixDense<N> &ContactMatrixDense<N>::operator=(const ContactMatrixDense<N> &other) {
  if (this == &other) {
    return *this;
  }

  const auto lck1 = this->lock();
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
    usize row, usize col) const {
  return std::unique_lock<ContactMatrixDense<N>::mutex_t>(
      this->_mtxes[this->get_pixel_mutex_idx(row, col)]);
}

template <class N>

utils::LockRangeExclusive<typename ContactMatrixDense<N>::mutex_t> ContactMatrixDense<N>::lock()
    const {
  return utils::LockRangeExclusive<mutex_t>(this->_mtxes);
}

template <class N>
constexpr usize ContactMatrixDense<N>::ncols() const {
  return this->_ncols;
}

template <class N>
constexpr usize ContactMatrixDense<N>::nrows() const {
  return this->_nrows;
}

template <class N>
constexpr usize ContactMatrixDense<N>::npixels() const {
  return this->_nrows * this->_ncols;
}

template <class N>
constexpr usize ContactMatrixDense<N>::get_n_of_missed_updates() const noexcept {
  return static_cast<usize>(this->_updates_missed.load());
}

template <class N>
constexpr usize ContactMatrixDense<N>::get_matrix_size_in_bytes() const {
  return this->npixels() * sizeof(N);
}

template <class N>
void ContactMatrixDense<N>::clear_missed_updates_counter() {
  this->_updates_missed = 0;
}

template <class N>
absl::Span<const N> ContactMatrixDense<N>::get_raw_count_vector() const {
  return absl::MakeConstSpan(this->_contacts);
}

template <class N>
absl::Span<N> ContactMatrixDense<N>::get_raw_count_vector() {
  return absl::MakeSpan(this->_contacts);
}

template <class N>
void ContactMatrixDense<N>::bound_check_coords([[maybe_unused]] const usize row,
                                               [[maybe_unused]] const usize col) const {
  if constexpr (utils::ndebug_not_defined()) {
    internal::bound_check_coords(*this, row, col);
  }
}

template <class N>
void ContactMatrixDense<N>::check_for_overflow_on_add(const usize row, const usize col,
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
void ContactMatrixDense<N>::check_for_overflow_on_subtract(usize row, usize col, const N n) const {
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

template <class N>
usize ContactMatrixDense<N>::hash_coordinates(const usize i, const usize j) noexcept {
  const std::array<usize, 2> buff{i, j};
  return utils::conditional_static_cast<usize>(XXH3_64bits(buff.data(), sizeof(usize) * 2));
}

template <class N>
constexpr usize ContactMatrixDense<N>::compute_number_of_mutexes(const usize rows,
                                                                 const usize cols) noexcept {
  if (rows + cols == 0) {
    return 0;
  }
  const auto nthreads = static_cast<usize>(std::thread::hardware_concurrency());
  const auto max_dim = std::max(rows, cols);
  // Clamping to 2-n is needed because get_pixel_mutex_idx expects the number of
  // mutexes to be a multiple of 2
  return utils::next_pow2(std::clamp(max_dim, usize(2), 1000 * nthreads));
}

template <class N>
usize ContactMatrixDense<N>::get_pixel_mutex_idx(const usize row, const usize col) const noexcept {
  assert(!this->_mtxes.empty());
  assert(this->_mtxes.size() % 2 == 0);
  // equivalent to hash_coordinates(row, col) % this->_mtxes.size() when _mtxes.size() % 2 == 0
  return (hash_coordinates(row, col) & (this->_mtxes.size() - 1));
}

}  // namespace modle

// IWYU pragma: private, include "modle/contact_matrix_dense.hpp"
// IWYU pragma: no_include <boost/core/checked_delete.hpp>
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <boost/exception/exception.hpp>
// IWYU pragma: no_include <boost/move/utility_core.hpp>
// IWYU pragma: no_include <H5Public.h>
