// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once
#include <absl/strings/str_split.h>  // for StrSplit, Splitter

#include <BS_thread_pool.hpp>  // for BS::thread_pool, parallelize_loop
#include <atomic>              // for atomic_fetch_add_explicit
#include <cassert>             // for assert
#include <cmath>               // for sqrt
#include <fstream>             // IWYU pragma: keep for ifstream
#include <vector>              // for vector

#include "modle/common/common.hpp"         // for usize, i64, u64, bp_t, isize
#include "modle/common/numeric_utils.hpp"  // for parse_numeric_or_throw
#include "modle/common/random.hpp"         // for PRNG, uniform_int_distribution, unifo...
#include "modle/common/suppress_compiler_warnings.hpp"
#include "modle/common/utils.hpp"                      // for ndebug_defined, ndebug_not_defined
#include "modle/compressed_io/compressed_io.hpp"       // for CompressedReader
#include "modle/internal/contact_matrix_internal.hpp"  // for transpose_coords
#include "modle/stats/misc.hpp"                        // for compute_gauss_kernel

namespace modle {

template <class N>
N ContactMatrixDense<N>::get(const usize row, const usize col) const {
  const auto [i, j] = internal::transpose_coords(row, col);
  this->bound_check_coords(i, j);

  if (i >= this->nrows()) {
    return 0;
  }

  const auto lck = this->lock_pixel(i, j);
  return this->unsafe_at(i, j);
}

template <class N>
void ContactMatrixDense<N>::set(const usize row, const usize col, const N n) {
  const auto [i, j] = internal::transpose_coords(row, col);
  this->bound_check_coords(i, j);

  if (i >= this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, usize(1), std::memory_order_relaxed);
    return;
  }

  const auto lck = this->lock_pixel(i, j);
  this->unsafe_at(i, j) = n;
  this->_global_stats_outdated = true;
}

template <class N>
void ContactMatrixDense<N>::add(const usize row, const usize col, const N n) {
  assert(n > 0);
  const auto [i, j] = internal::transpose_coords(row, col);
  this->bound_check_coords(i, j);

  if (i >= this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, usize(1), std::memory_order_relaxed);
    return;
  }

  const auto lck = this->lock_pixel(i, j);
  this->unsafe_at(i, j) += n;
  this->_global_stats_outdated = true;
}

template <class N>
void ContactMatrixDense<N>::subtract(const usize row, const usize col, const N n) {
  assert(n >= 0);
  const auto [i, j] = internal::transpose_coords(row, col);
  this->bound_check_coords(i, j);

  if (i >= this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, usize(1), std::memory_order_relaxed);
    return;
  }

  const auto lck = this->lock_pixel(i, j);
  this->unsafe_at(i, j) -= n;
  this->_global_stats_outdated = true;
}

template <class N>
void ContactMatrixDense<N>::increment(usize row, usize col) {
  this->add(row, col, N(1));
}

template <class N>
void ContactMatrixDense<N>::decrement(usize row, usize col) {
  this->subtract(row, col, N(1));
}

template <class N>
double ContactMatrixDense<N>::get_fraction_of_missed_updates() const {
  if (this->empty() || this->get_n_of_missed_updates() == N(0)) {
    return 0.0;
  }
  const auto lck = this->lock();
  const auto missed_updates = static_cast<double>(this->get_n_of_missed_updates());
  return missed_updates / (static_cast<double>(this->unsafe_get_tot_contacts()) + missed_updates);
}

template <class N>
double ContactMatrixDense<N>::get_avg_contact_density() const {
  return static_cast<double>(this->get_tot_contacts()) / static_cast<double>(this->npixels());
}

template <class N>
auto ContactMatrixDense<N>::get_tot_contacts() const -> SumT {
  if (this->_global_stats_outdated) {
    const auto lck = this->lock();
    return this->unsafe_get_tot_contacts();
  }
  return this->_tot_contacts.load();
}

template <class N>
usize ContactMatrixDense<N>::get_nnz() const {
  if (this->_global_stats_outdated) {
    const auto lck = this->lock();
    return this->unsafe_get_nnz();
  }
  return this->_nnz.load();
}

template <class N>
N ContactMatrixDense<N>::get_min_count() const noexcept {
  if (this->get_tot_contacts() == 0) {
    return 0;
  }
  const auto lck = this->lock();
  return *std::min_element(this->_contacts.begin(), this->_contacts.end());
}

template <class N>
N ContactMatrixDense<N>::get_max_count() const noexcept {
  if (this->get_tot_contacts() == 0) {
    return 0;
  }
  const auto lck = this->lock();
  return *std::max_element(this->_contacts.begin(), this->_contacts.end());
}

template <class N>
void ContactMatrixDense<N>::reset() {
  const auto lck = this->exclusive_lock();
  this->unsafe_reset();
}

template <class N>
ContactMatrixDense<double> ContactMatrixDense<N>::blur(const double sigma, const double truncate,
                                                       BS::thread_pool* tpool) const {
  ContactMatrixDense<double> bmatrix(this->nrows(), this->ncols());
  if (this->empty()) {
    return bmatrix;
  }

  const auto gauss_kernel = stats::compute_gauss_kernel2d(sigma, truncate);

  const auto block_size = static_cast<usize>(std::sqrt(static_cast<double>(gauss_kernel.size())));
  assert(block_size * block_size == gauss_kernel.size());

  auto apply_kernel = [this, &bmatrix, &gauss_kernel, block_size](const usize i0, const usize i1) {
    std::vector<N> pixels(gauss_kernel.size());
    for (usize i = i0; i < i1; ++i) {
      for (usize j = i; j < std::min(i + this->nrows(), this->ncols()); ++j) {
        this->unsafe_get_block(i, j, block_size, pixels);
        bmatrix.unsafe_set(i, j, stats::cross_correlation(gauss_kernel, pixels));
      }
    }
  };

  const auto lck = this->lock();
  if (tpool) {
    auto fut = tpool->template parallelize_loop(usize(0), this->ncols(), apply_kernel);
    fut.wait();
  } else {
    apply_kernel(0, this->ncols());
  }
  return bmatrix;
}

template <class N>
ContactMatrixDense<double> ContactMatrixDense<N>::diff_of_gaussians(
    const double sigma1, const double sigma2, const double truncate, const double min_value,
    const double max_value, BS::thread_pool* tpool) const {
  assert(sigma1 <= sigma2);
  ContactMatrixDense<double> bmatrix(this->nrows(), this->ncols());
  if (this->empty()) {
    return bmatrix;
  }

  const auto gauss_kernel1 = stats::compute_gauss_kernel2d(sigma1, truncate);
  const auto gauss_kernel2 = stats::compute_gauss_kernel2d(sigma2, truncate);

  [[maybe_unused]] const auto block_size1 =
      static_cast<usize>(std::sqrt(static_cast<double>(gauss_kernel1.size())));
  [[maybe_unused]] const auto block_size2 =
      static_cast<usize>(std::sqrt(static_cast<double>(gauss_kernel2.size())));
  assert(block_size1 * block_size1 <= gauss_kernel1.size());
  assert(block_size2 * block_size2 <= gauss_kernel2.size());

  auto compute_gauss_diff = [this, &bmatrix, &gauss_kernel1, &gauss_kernel2, block_size1,
                             block_size2, min_value, max_value](const usize i0, const usize i1) {
    std::vector<N> pixels(std::max(gauss_kernel1.size(), gauss_kernel2.size()));
    for (usize i = i0; i < i1; ++i) {
      for (usize j = i; j < std::min(i + this->nrows(), this->ncols()); ++j) {
        this->unsafe_get_block(i, j, block_size1, pixels);
        const auto n1 = stats::cross_correlation(gauss_kernel1, pixels);
        this->unsafe_get_block(i, j, block_size2, pixels);
        const auto n2 = stats::cross_correlation(gauss_kernel2, pixels);

        bmatrix.set(i, j, std::clamp(n1 - n2, min_value, max_value));
      }
    }
  };

  const auto lck = this->lock();
  if (tpool) {
    auto fut = tpool->template parallelize_loop(usize(0), this->ncols(), compute_gauss_diff);
    fut.wait();
  } else {
    compute_gauss_diff(0, this->ncols());
  }
  return bmatrix;
}

template <class N>
template <class FP, class>
ContactMatrixDense<FP> ContactMatrixDense<N>::normalize(const double lb, const double ub) const {
  const auto lck = this->lock();
  return this->unsafe_normalize(lb, ub);
}

template <class N>
inline void ContactMatrixDense<N>::normalize_inplace(const N lb, const N ub) noexcept {
  const auto lck = this->lock();
  this->unsafe_normalize_inplace(lb, ub);
}

template <class N>
ContactMatrixDense<N> ContactMatrixDense<N>::clamp(const N lb, const N ub) const {
  const auto lck = this->lock();
  ContactMatrixDense<N> m(this->nrows(), this->ncols());
  ContactMatrixDense<N>::unsafe_clamp(*this, m, lb, ub);
  return m;
}

template <class N>
inline void ContactMatrixDense<N>::clamp_inplace(const N lb, const N ub) noexcept {
  const auto lck = this->lock();
  ContactMatrixDense<N>::unsafe_clamp(*this, *this, lb, ub);
}

template <class N>
template <class N1, class N2>
ContactMatrixDense<N1> ContactMatrixDense<N>::discretize(const IITree<N2, N1>& mappings) const {
  const auto lck = this->lock();
  return this->unsafe_discretize<N1>(mappings);
}

template <class N>
template <class M>
void ContactMatrixDense<N>::discretize_inplace(const IITree<M, N>& mappings) noexcept {
  const auto lck = this->lock();
  this->unsafe_discretize_inplace(mappings);
}

template <class N>
template <class M, class>
ContactMatrixDense<M> ContactMatrixDense<N>::as() const {
  const auto lck = this->lock();
  return this->unsafe_as<M>();
}

template <class N>
bool ContactMatrixDense<N>::empty() const {
  return this->get_tot_contacts() == 0;
}

template <class N>
ContactMatrixDense<N> ContactMatrixDense<N>::from_txt(const std::filesystem::path& path, char sep) {
  assert(std::filesystem::exists(path));
  compressed_io::Reader r(path);

  ContactMatrixDense<N> m{};

  std::string buff;
  std::vector<std::string_view> toks;
  usize i = 0;
  for (; r.getline(buff); ++i) {
    toks = absl::StrSplit(buff, sep);
    if (i == 0) {
      m.unsafe_resize(toks.size(), toks.size());
    }
    for (usize j = i; j < m.ncols(); ++j) {
      m.unsafe_set(i, j, utils::parse_numeric_or_throw<N>(toks[j]));
    }
  }
  assert(i != 0);

  return m;
}

}  // namespace modle

// IWYU pragma: private, include "modle/contact_matrix_dense.hpp"
// IWYU pragma: no_include <boost/core/checked_delete.hpp>
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <boost/exception/exception.hpp>
// IWYU pragma: no_include <boost/move/utility_core.hpp>
// IWYU pragma: no_include <H5Public.h>
