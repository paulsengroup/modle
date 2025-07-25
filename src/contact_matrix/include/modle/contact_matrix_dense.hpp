// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <BS_thread_pool.hpp>
#include <atomic>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <filesystem>
#include <iostream>
#include <limits>
#include <mutex>
#include <span>
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
  using SumT = typename std::conditional<std::is_floating_point_v<N>, double, std::int64_t>::type;
  std::uint64_t _nrows{0};
  std::uint64_t _ncols{0};
  std::vector<N> _contacts{};
  mutable std::vector<mutex_t> _mtxes{};
  mutable std::atomic<SumT> _tot_contacts{0};
  mutable std::atomic<std::size_t> _nnz{0};
  mutable std::atomic<bool> _global_stats_outdated{false};
  std::atomic<std::size_t> _updates_missed{0};

 public:
  // Constructors
  ContactMatrixDense() = default;
  ContactMatrixDense(ContactMatrixDense<N>&& other) noexcept = default;
  inline ContactMatrixDense(const ContactMatrixDense<N>& other);
  inline ContactMatrixDense(std::size_t nrows_, std::size_t ncols_);
  inline ContactMatrixDense(bp_t length, bp_t diagonal_width, bp_t bin_size);
  inline ContactMatrixDense(std::span<const N> contacts, std::size_t nrows_, std::size_t ncols_,
                            std::size_t tot_contacts = 0, std::size_t updates_missed = 0);
  ~ContactMatrixDense() = default;

  // Operators
  inline ContactMatrixDense<N>& operator=(const ContactMatrixDense<N>& other);
  ContactMatrixDense<N>& operator=(ContactMatrixDense<N>&& other) noexcept = default;

  // Thread-safe count getters and setters
  [[nodiscard]] inline N get(std::size_t row, std::size_t col) const;
  inline void set(std::size_t row, std::size_t col, N n);
  inline void add(std::size_t row, std::size_t col, N n);
  inline void subtract(std::size_t row, std::size_t col, N n);
  inline void increment(std::size_t row, std::size_t col);
  inline void decrement(std::size_t row, std::size_t col);

  // Thread-UNsafe count getters and setters
  [[nodiscard]] inline N unsafe_get(std::size_t row, std::size_t col) const;
  inline void unsafe_set(std::size_t row, std::size_t col, N n);
  inline void unsafe_add(std::size_t row, std::size_t col, N n);
  inline void unsafe_subtract(std::size_t row, std::size_t col, N n);
  inline void unsafe_increment(std::size_t row, std::size_t col);
  inline void unsafe_decrement(std::size_t row, std::size_t col);

  inline void unsafe_get_column(std::size_t col, std::vector<N>& buff,
                                std::size_t row_offset = 0) const;
  inline void unsafe_get_row(std::size_t row, std::vector<N>& buff,
                             std::size_t col_offset = 0) const;
  // block_size is required to be an odd number at the moment
  inline void unsafe_get_block(std::size_t row, std::size_t col, std::size_t block_size,
                               std::vector<N>& buff) const;
  [[nodiscard]] inline N unsafe_get_block(std::size_t row, std::size_t col,
                                          std::size_t block_size) const;

  // Shape/statistics getters
  [[nodiscard]] constexpr std::size_t ncols() const;
  [[nodiscard]] constexpr std::size_t nrows() const;
  [[nodiscard]] constexpr std::size_t npixels() const;
  [[nodiscard]] constexpr std::size_t get_n_of_missed_updates() const noexcept;
  [[nodiscard]] inline double get_fraction_of_missed_updates() const;
  [[nodiscard]] inline SumT get_tot_contacts() const;
  [[nodiscard]] inline std::size_t get_nnz() const;
  [[nodiscard]] inline double get_avg_contact_density() const;
  [[nodiscard]] constexpr std::size_t get_matrix_size_in_bytes() const;
  [[nodiscard]] inline N get_min_count() const noexcept;
  [[nodiscard]] inline N get_max_count() const noexcept;

  [[nodiscard]] constexpr double unsafe_get_fraction_of_missed_updates() const noexcept;
  [[nodiscard]] inline SumT unsafe_get_tot_contacts() const noexcept;
  [[nodiscard]] inline std::size_t unsafe_get_nnz() const noexcept;
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
  inline void unsafe_resize(std::size_t nrows, std::size_t ncols);
  inline void unsafe_resize(bp_t length, bp_t diagonal_width, bp_t bin_size);
  [[nodiscard]] inline bool empty() const;
  [[nodiscard]] inline bool unsafe_empty() const;
  [[nodiscard]] inline std::span<const N> get_raw_count_vector() const;
  [[nodiscard]] inline std::span<N> get_raw_count_vector();

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
  [[nodiscard]] inline N& unsafe_at(std::size_t i, std::size_t j);
  [[nodiscard]] inline const N& unsafe_at(std::size_t i, std::size_t j) const;

  inline void bound_check_coords(std::size_t row, std::size_t col) const;

  [[nodiscard]] inline utils::LockRangeExclusive<mutex_t> lock() const;
  [[nodiscard]] inline std::unique_lock<mutex_t> lock_pixel(std::size_t row, std::size_t col) const;
  [[nodiscard]] static inline std::size_t hash_coordinates(std::size_t i, std::size_t j) noexcept;
  [[nodiscard]] static constexpr std::size_t compute_number_of_mutexes(std::size_t rows,
                                                                       std::size_t cols) noexcept;
  [[nodiscard]] inline std::size_t get_pixel_mutex_idx(std::size_t row,
                                                       std::size_t col) const noexcept;

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
