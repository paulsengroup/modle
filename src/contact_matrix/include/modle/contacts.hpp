// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span

#include <atomic>                                   // for atomic
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/filesystem/path.hpp>                // for path
#include <iostream>                                 // for cout, ostream
#include <limits>                                   // for numeric_limits
#include <mutex>                                    // for unique_lock
#include <shared_mutex>                             // for shared_mutex, shared_lock
#include <thread_pool/thread_pool.hpp>              // for thread_pool
#include <type_traits>                              // for enable_if_t
#include <utility>                                  // for pair
#include <vector>                                   // for vector

#include "modle/common/common.hpp"  // for usize, bp_t, u64, contacts_t, i64

namespace modle {

template <class N, class T>
class IITree;

template <class N = contacts_t>
class ContactMatrix {
  static_assert(std::is_arithmetic_v<N>,
                "ContactMatrix requires a numeric type as template argument.");

 public:
  using value_type = N;
  using pixel_write_lock =
      std::pair<std::shared_lock<std::shared_mutex>, std::unique_lock<std::shared_mutex>>;
  using pixel_read_lock =
      std::pair<std::shared_lock<std::shared_mutex>, std::shared_lock<std::shared_mutex>>;

  // This allows methods from different instantiations of the ContactMatrix template class:
  // Example: allowing ContactMatrix<int> to access private member variables and functions from
  // class ContactMatrix<double>
  template <class M>
  friend class ContactMatrix;

 private:
  u64 _nrows{0};
  u64 _ncols{0};
  std::vector<N> _contacts{};
  mutable std::vector<std::shared_mutex> _mtxes{};
  mutable std::shared_mutex _global_mtx{};
  mutable std::atomic<N> _tot_contacts{0};
  mutable std::atomic<bool> _tot_contacts_outdated{false};
  std::atomic<i64> _updates_missed{0};

 public:
  // Constructors
  inline ContactMatrix() = default;
#if defined(__clang__) && __clang_major__ < 7
  inline ContactMatrix(ContactMatrix<N>&& other) = default;
#else
  inline ContactMatrix(ContactMatrix<N>&& other) noexcept = default;
#endif
  inline ContactMatrix(const ContactMatrix<N>& other);
  inline ContactMatrix(usize nrows, usize ncols, bool fill_with_random_numbers = false,
                       u64 seed = 8336046165695760686);
  inline ContactMatrix(bp_t length, bp_t diagonal_width, bp_t bin_size,
                       bool fill_with_random_numbers = false);
  inline ContactMatrix(absl::Span<const N> contacts, usize nrows, usize ncols,
                       usize tot_contacts = 0, usize updates_missed = 0);
  inline ~ContactMatrix() = default;

  // Operators
  inline ContactMatrix<N>& operator=(const ContactMatrix<N>& other);
#if defined(__clang__) && __clang_major__ < 7
  inline ContactMatrix<N>& operator=(ContactMatrix<N>&& other) = default;
#else
  inline ContactMatrix<N>& operator=(ContactMatrix<N>&& other) noexcept = default;
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
  [[nodiscard]] inline constexpr usize ncols() const;
  [[nodiscard]] inline constexpr usize nrows() const;
  [[nodiscard]] inline constexpr usize npixels() const;
  [[nodiscard]] inline usize unsafe_npixels_after_masking() const;
  [[nodiscard]] inline constexpr usize get_n_of_missed_updates() const noexcept;
  [[nodiscard]] inline constexpr double unsafe_get_fraction_of_missed_updates() const noexcept;
  [[nodiscard]] inline double get_fraction_of_missed_updates() const;
  [[nodiscard]] inline usize get_tot_contacts() const;
  [[nodiscard]] inline usize unsafe_get_tot_contacts() const noexcept;
  [[nodiscard]] inline double get_avg_contact_density() const;
  [[nodiscard]] inline double unsafe_get_avg_contact_density() const;
  [[nodiscard]] inline constexpr usize get_matrix_size_in_bytes() const;
  [[nodiscard]] inline constexpr double get_matrix_size_in_mb() const;
  [[nodiscard]] inline N get_min_count() const noexcept;
  [[nodiscard]] inline N unsafe_get_min_count() const noexcept;
  [[nodiscard]] inline N get_max_count() const noexcept;
  [[nodiscard]] inline N unsafe_get_max_count() const noexcept;

  // Debug
  inline void unsafe_print(std::ostream& out_stream = std::cout, bool full = false) const;
  inline void unsafe_print(bool full) const;
  [[nodiscard]] inline std::vector<std::vector<N>> unsafe_generate_symmetric_matrix() const;
  inline void unsafe_import_from_txt(const boost::filesystem::path& path, char sep = '\t');

  // Mask operations
  inline void unsafe_generate_mask_for_bins_without_contacts(boost::dynamic_bitset<>& mask) const;
  [[nodiscard]] inline boost::dynamic_bitset<> unsafe_generate_mask_for_bins_without_contacts()
      const;

  // Misc
  inline void clear_missed_updates_counter();
  inline void reset();
  inline void unsafe_reset();
  // Note: upon resizing the content of a ContactMatrix is invalidated
  inline void unsafe_resize(usize nrows, usize ncols);
  inline void unsafe_resize(bp_t length, bp_t diagonal_width, bp_t bin_size);
  [[nodiscard]] inline bool empty() const;
  [[nodiscard]] inline bool unsafe_empty() const;
  [[nodiscard]] inline absl::Span<const N> get_raw_count_vector() const;
  [[nodiscard]] inline absl::Span<N> get_raw_count_vector();
  inline void unsafe_compute_row_wise_contact_histogram(std::vector<u64>& buff) const;
  [[nodiscard]] inline std::vector<u64> unsafe_compute_row_wise_contact_histogram() const;
  inline void unsafe_deplete_contacts(double depletion_multiplier = 1.0);

  [[nodiscard]] inline ContactMatrix<double> blur(double sigma, double cutoff = 0.005,
                                                  thread_pool* tpool = nullptr) const;
  [[nodiscard]] inline ContactMatrix<double> gaussian_diff(
      double sigma1, double sigma2, double min_value = std::numeric_limits<double>::lowest(),
      double max_value = (std::numeric_limits<double>::max)(), thread_pool* tpool = nullptr) const;
  template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
  [[nodiscard]] inline ContactMatrix<FP> unsafe_normalize(double lb = 0.0, double ub = 1.0) const;
  template <class FP = double, class = std::enable_if_t<std::is_floating_point_v<FP>>>
  [[nodiscard]] inline ContactMatrix<FP> normalize(double lb = 0.0, double ub = 1.0) const;
  [[nodiscard]] inline ContactMatrix<N> unsafe_clamp(N lb, N ub) const;
  [[nodiscard]] inline ContactMatrix<N> clamp(N lb, N ub) const;
  template <class N1, class N2>
  [[nodiscard]] inline ContactMatrix<N1> discretize(const IITree<N2, N1>& mappings) const;
  template <class N1, class N2>
  [[nodiscard]] inline ContactMatrix<N1> unsafe_discretize(const IITree<N2, N1>& mappings) const;

  // Convert a matrix of type N to a matrix of type M
  // When N is a floating point type and M isn't, contacts are round before casting them to M
  template <class M, class = std::enable_if_t<!std::is_same_v<N, M>>>
  [[nodiscard]] inline ContactMatrix<M> as() const;
  template <class M, class = std::enable_if_t<!std::is_same_v<N, M>>>
  [[nodiscard]] inline ContactMatrix<M> unsafe_as() const;

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

  [[nodiscard]] constexpr static std::pair<usize, usize> transpose_coords(usize row,
                                                                          usize col) noexcept;
  inline void bound_check_coords(usize row, usize col) const;
  inline void check_for_overflow_on_add(usize row, usize col, N n) const;
  inline void check_for_overflow_on_subtract(usize row, usize col, N n) const;

  [[nodiscard]] inline std::unique_lock<std::shared_mutex> lock() const;
  [[nodiscard]] inline pixel_read_lock lock_pixel_read(usize row, usize col) const;
  [[nodiscard]] inline pixel_write_lock lock_pixel_write(usize row, usize col) const;
  [[nodiscard]] static inline usize hash_coordinates(usize i, usize j) noexcept;
  [[nodiscard]] static constexpr usize compute_number_of_mutexes(usize rows, usize cols) noexcept;
  [[nodiscard]] inline usize get_pixel_mutex_idx(usize row, usize col) const noexcept;

  template <class M, class = std::enable_if_t<std::is_floating_point_v<M>>>
  static inline void unsafe_normalize(const ContactMatrix<N>& input_matrix,
                                      ContactMatrix<M>& output_matrix, M lb = 0, M ub = 1) noexcept;

  static inline void unsafe_clamp(const ContactMatrix<N>& input_matrix,
                                  ContactMatrix<N>& output_matrix, N lb, N ub) noexcept;
  template <class N1, class N2>
  static inline void unsafe_discretize(const ContactMatrix<N>& input_matrix,
                                       ContactMatrix<N1>& output_matrix,
                                       const IITree<N2, N1>& mappings) noexcept;
};
}  // namespace modle

#include "../../contacts_impl.hpp"         // IWYU pragma: export
#include "../../contacts_safe_impl.hpp"    // IWYU pragma: export
#include "../../contacts_unsafe_impl.hpp"  // IWYU pragma: export
// IWYU pragma: "../../contacts_impl.hpp"
// IWYU pragma: "../../contacts_safe_impl.hpp"
// IWYU pragma: "../../contacts_unsafe_impl.hpp"
