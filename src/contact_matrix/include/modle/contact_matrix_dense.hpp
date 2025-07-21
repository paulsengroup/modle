// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>

#include <BS_thread_pool.hpp>
#include <atomic>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <filesystem>
#include <iostream>
#include <limits>
#include <mutex>
#include <type_traits>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/common/utils.hpp"

namespace modle {

template <class N, class T>
class IITree;

template <class N = contacts_t>
class ContactMatrixDense {
  static_assert(std::is_arithmetic_v<N>,
                "ContactMatrixDense requires a numeric type as template argument.");

 public:
  using value_type = N;

  // This allows methods from different instantiations of the ContactMatrixDense template class:
  // Example: allowing ContactMatrixDense<int> to access private member variables and functions from
  // class ContactMatrixDense<double>
  template <class M>
  friend class ContactMatrixDense;

 private:
  using mutex_t = std::mutex;
  using SumT = typename std::conditional<std::is_floating_point_v<N>, double, i64>::type;
  u64 _nrows{0};
  u64 _ncols{0};
  std::vector<N> _contacts{};
  mutable std::vector<mutex_t> _mtxes{};
  mutable std::atomic<SumT> _tot_contacts{0};
  mutable std::atomic<usize> _nnz{0};
  mutable std::atomic<bool> _global_stats_outdated{false};
  std::atomic<usize> _updates_missed{0};

 public:
  // Constructors
  ContactMatrixDense() = default;
#if defined(__clang__) && __clang_major__ < 9
  ContactMatrixDense(ContactMatrixDense<N>&& other) = default;
#else
  ContactMatrixDense(ContactMatrixDense<N>&& other) noexcept = default;
#endif
  inline ContactMatrixDense(const ContactMatrixDense<N>& other);
  inline ContactMatrixDense(usize nrows, usize ncols);
  inline ContactMatrixDense(bp_t length, bp_t diagonal_width, bp_t bin_size);
  inline ContactMatrixDense(absl::Span<const N> contacts, usize nrows, usize ncols,
                            usize tot_contacts = 0, usize updates_missed = 0);
  ~ContactMatrixDense() = default;

  // Operators
  inline ContactMatrixDense<N>& operator=(const ContactMatrixDense<N>& other);
#if defined(__clang__) && __clang_major__ < 9
  ContactMatrixDense<N>& operator=(ContactMatrixDense<N>&& other) = default;
#else
  ContactMatrixDense<N>& operator=(ContactMatrixDense<N>&& other) noexcept = default;
#endif

  // Thread-safe count getters and setters
  [[nodiscard]] inline N get(usize row, usize col) const;
  inline void set(usize row, usize col, N n);
  inline void add(usize row, usize col, N n);
  inline void subtract(usize row, usize col, N n);
  inline void increment(usize row, usize col);
  inline void decrement(usize row, usize col);

  // Thread-UNsafe count getters and setters
  [[nodiscard]] inline N unsafe_get(usize row, usize col) const;
  inline void unsafe_set(usize row, usize col, N n);
  inline void unsafe_add(usize row, usize col, N n);
  inline void unsafe_subtract(usize row, usize col, N n);
  inline void unsafe_increment(usize row, usize col);
  inline void unsafe_decrement(usize row, usize col);

  inline void unsafe_get_column(usize col, std::vector<N>& buff, usize row_offset = 0) const;
  inline void unsafe_get_row(usize row, std::vector<N>& buff, usize col_offset = 0) const;
  // block_size is required to be an odd number at the moment
  inline void unsafe_get_block(usize row, usize col, usize block_size, std::vector<N>& buff) const;
  [[nodiscard]] inline N unsafe_get_block(usize row, usize col, usize block_size) const;

  // Shape/statistics getters
  [[nodiscard]] constexpr usize ncols() const;
  [[nodiscard]] constexpr usize nrows() const;
  [[nodiscard]] constexpr usize npixels() const;
  [[nodiscard]] constexpr usize get_n_of_missed_updates() const noexcept;
  [[nodiscard]] inline double get_fraction_of_missed_updates() const;
  [[nodiscard]] inline SumT get_tot_contacts() const;
  [[nodiscard]] inline usize get_nnz() const;
  [[nodiscard]] inline double get_avg_contact_density() const;
  [[nodiscard]] constexpr usize get_matrix_size_in_bytes() const;
  [[nodiscard]] inline N get_min_count() const noexcept;
  [[nodiscard]] inline N get_max_count() const noexcept;

  [[nodiscard]] constexpr double unsafe_get_fraction_of_missed_updates() const noexcept;
  [[nodiscard]] inline SumT unsafe_get_tot_contacts() const noexcept;
  [[nodiscard]] inline usize unsafe_get_nnz() const noexcept;
  [[nodiscard]] inline double unsafe_get_avg_contact_density() const;
  [[nodiscard]] inline N unsafe_get_min_count() const noexcept;
  [[nodiscard]] inline N unsafe_get_max_count() const noexcept;

  // Debug
  inline void unsafe_print(std::ostream& out_stream = std::cout, bool full = false) const;
  inline void unsafe_print(bool full) const;
  [[nodiscard]] inline std::vector<std::vector<N>> unsafe_generate_symmetric_matrix() const;
  [[nodiscard]] static inline ContactMatrixDense<N> from_txt(const std::filesystem::path& path,
                                                             char sep = '\t');

  // Misc
  inline void clear_missed_updates_counter();
  inline void unsafe_reset();
  // Note: upon resizing the content of a ContactMatrixDense is invalidated
  inline void unsafe_resize(usize nrows, usize ncols);
  inline void unsafe_resize(bp_t length, bp_t diagonal_width, bp_t bin_size);
  [[nodiscard]] inline bool empty() const;
  [[nodiscard]] inline bool unsafe_empty() const;
  [[nodiscard]] inline absl::Span<const N> get_raw_count_vector() const;
  [[nodiscard]] inline absl::Span<N> get_raw_count_vector();

  [[nodiscard]] inline ContactMatrixDense<double> blur(
      double sigma, double truncate = 3.5, BS::light_thread_pool* tpool = nullptr) const;
  [[nodiscard]] inline ContactMatrixDense<double> diff_of_gaussians(
      const double sigma1, const double sigma2, const double truncate = 3.5,
      const double min_value = std::numeric_limits<double>::lowest(),
      const double max_value = (std::numeric_limits<double>::max)(),
      BS::light_thread_pool* tpool = nullptr) const;

  template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
  [[nodiscard]] inline ContactMatrixDense<FP> unsafe_normalize(double lb = 0.0,
                                                               double ub = 1.0) const;
  template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
  [[nodiscard]] inline ContactMatrixDense<FP> normalize(double lb = 0.0, double ub = 1.0) const;
  [[nodiscard]] inline ContactMatrixDense<N> unsafe_clamp(N lb, N ub) const;
  [[nodiscard]] inline ContactMatrixDense<N> clamp(N lb, N ub) const;
  template <class N1, class N2>
  [[nodiscard]] inline ContactMatrixDense<N1> discretize(const IITree<N2, N1>& mappings) const;
  template <class N1, class N2>
  [[nodiscard]] inline ContactMatrixDense<N1> unsafe_discretize(
      const IITree<N2, N1>& mappings) const;

  // Convert a matrix of type N to a matrix of type M
  // When N is a floating point type and M isn't, contacts are round before casting them to M
  template <class M, class = std::enable_if_t<!std::is_same_v<N, M>>>
  [[nodiscard]] inline ContactMatrixDense<M> as() const;
  template <class M, class = std::enable_if_t<!std::is_same_v<N, M>>>
  [[nodiscard]] inline ContactMatrixDense<M> unsafe_as() const;

  inline void unsafe_normalize_inplace(N lb = 0, N ub = 1) noexcept;
  inline void normalize_inplace(N lb = 0, N ub = 1) noexcept;
  inline void unsafe_clamp_inplace(N lb, N ub) noexcept;
  inline void clamp_inplace(N lb, N ub) noexcept;
  template <class M>
  inline void discretize_inplace(const IITree<M, N>& mappings) noexcept;
  template <class M>
  inline void unsafe_discretize_inplace(const IITree<M, N>& mappings) noexcept;

 private:
  [[nodiscard]] inline N& unsafe_at(usize i, usize j);
  [[nodiscard]] inline const N& unsafe_at(usize i, usize j) const;

  inline void bound_check_coords(usize row, usize col) const;

  [[nodiscard]] inline utils::LockRangeExclusive<mutex_t> lock() const;
  [[nodiscard]] inline std::unique_lock<mutex_t> lock_pixel(usize row, usize col) const;
  [[nodiscard]] static inline usize hash_coordinates(usize i, usize j) noexcept;
  [[nodiscard]] static constexpr usize compute_number_of_mutexes(usize rows, usize cols) noexcept;
  [[nodiscard]] inline usize get_pixel_mutex_idx(usize row, usize col) const noexcept;

  template <class M, class = std::enable_if_t<std::is_floating_point_v<M>>>
  static inline void unsafe_normalize(const ContactMatrixDense<N>& input_matrix,
                                      ContactMatrixDense<M>& output_matrix, M lb = 0,
                                      M ub = 1) noexcept;

  static inline void unsafe_clamp(const ContactMatrixDense<N>& input_matrix,
                                  ContactMatrixDense<N>& output_matrix, N lb, N ub) noexcept;
  template <class N1, class N2>
  static inline void unsafe_discretize(const ContactMatrixDense<N>& input_matrix,
                                       ContactMatrixDense<N1>& output_matrix,
                                       const IITree<N2, N1>& mappings) noexcept;

  inline void unsafe_update_global_stats() const noexcept;
};
}  // namespace modle

#include "../../contact_matrix_dense_impl.hpp"         // IWYU pragma: export
#include "../../contact_matrix_dense_safe_impl.hpp"    // IWYU pragma: export
#include "../../contact_matrix_dense_unsafe_impl.hpp"  // IWYU pragma: export
// IWYU pragma: "../../contact_matrix_dense_impl.hpp"
// IWYU pragma: "../../contact_matrix_dense_safe_impl.hpp"
// IWYU pragma: "../../contact_matrix_dense_unsafe_impl.hpp"
