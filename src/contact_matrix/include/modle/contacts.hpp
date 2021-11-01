// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span

#include <atomic>                                   // for atomic
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/filesystem/path.hpp>                // for path
#include <iostream>                                 // for cout, ostream
#include <mutex>                                    // for mutex
#include <utility>                                  // for pair
#include <vector>                                   // for vector

#include "modle/common/common.hpp"  // for usize, bp_t, u64, i64, contacts_t
#include "modle/common/utils.hpp"   // for ndebug_defined

namespace modle {

template <class N = contacts_t>
class ContactMatrix {
  static_assert(std::is_arithmetic_v<N>,
                "ContactMatrix requires a numeric type as template argument.");

 public:
  // Constructors
  inline ContactMatrix() = default;
#if defined(__clang__) && __clang_major__ < 7
  inline ContactMatrix(ContactMatrix<N>&& other) noexcept;
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
  [[nodiscard]] inline ContactMatrix<N>& operator=(const ContactMatrix<N>& other);
#if defined(__clang__) && __clang_major__ < 7
  [[nodiscard]] inline ContactMatrix<N>& operator=(ContactMatrix<N>&& other) noexcept;
#else
  [[nodiscard]] inline ContactMatrix<N>& operator=(ContactMatrix<N>&& other) noexcept = default;
#endif

  // Counts getters and setters
  [[nodiscard]] inline N unsafe_get(usize row, usize col) const noexcept(utils::ndebug_defined());
  inline void unsafe_get_column(usize col, std::vector<N>& buff, usize row_offset = 0) const
      noexcept(utils::ndebug_defined());
  inline void unsafe_get_row(usize row, std::vector<N>& buff, usize col_offset = 0) const
      noexcept(utils::ndebug_defined());

  // block_size is required to be an odd number at the moment
  [[nodiscard]] inline N unsafe_get(usize row, usize col, usize block_size) const
      noexcept(utils::ndebug_defined());
  inline void unsafe_get(usize row, usize col, usize block_size, std::vector<N>& buff) const
      noexcept(utils::ndebug_defined());

  inline void unsafe_set(usize row, usize col, N n) noexcept(utils::ndebug_defined());
  inline void set(usize row, usize col, N n) noexcept(utils::ndebug_defined());

  // The following four methods are thread-safe but operations are not atomic, as without locking it
  // is not possible to ensure that updates to pixels and total counters are atomic.
  // In order to avoid overflows, special care should be taken when N is unsigned and add/subtract
  // or increment/decrement are called from different threads
  inline void add(usize row, usize col, N n) noexcept(utils::ndebug_defined());
  inline void subtract(usize row, usize col, N) noexcept(utils::ndebug_defined());
  inline void increment(usize row, usize col) noexcept(utils::ndebug_defined());
  inline void decrement(usize row, usize col) noexcept(utils::ndebug_defined());

  // Shape/statistics getters
  [[nodiscard]] inline constexpr usize ncols() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr usize nrows() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr usize npixels() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline usize unsafe_npixels_after_masking() const;
  [[nodiscard]] inline constexpr usize get_n_of_missed_updates() const noexcept;
  [[nodiscard]] inline constexpr double unsafe_get_fraction_of_missed_updates() const noexcept;
  [[nodiscard]] inline constexpr usize get_tot_contacts() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline double get_avg_contact_density() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr usize get_matrix_size_in_bytes() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr double get_matrix_size_in_mb() const
      noexcept(utils::ndebug_defined());

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
  inline void unsafe_reset();
  // Note: upon resizing the content of a ContactMatrix is invalidated
  inline void unsafe_resize(usize nrows, usize ncols);
  inline void unsafe_resize(bp_t length, bp_t diagonal_width, bp_t bin_size);
  [[nodiscard]] inline bool empty() const;
  [[nodiscard]] inline absl::Span<const N> get_raw_count_vector() const;
  [[nodiscard]] inline absl::Span<N> get_raw_count_vector();
  inline void unsafe_compute_row_wise_contact_histogram(std::vector<u64>& buff) const;
  [[nodiscard]] inline std::vector<u64> unsafe_compute_row_wise_contact_histogram() const;
  inline void unsafe_deplete_contacts(double depletion_multiplier = 1.0);
  [[nodiscard]] inline ContactMatrix<double> blur(double sigma = 1, double cutoff = 0.005) const
      noexcept(utils::ndebug_defined());

  using value_type = N;

 private:
  u64 _nrows{0};
  u64 _ncols{0};
  std::vector<N> _contacts{};
  std::vector<std::mutex> _mtxes{};  // Each mutex protects one row of bins
  std::atomic<N> _tot_contacts{0};
  std::atomic<i64> _updates_missed{0};

  [[nodiscard]] inline N& at(usize i, usize j) noexcept(utils::ndebug_defined());
  [[nodiscard]] inline const N& at(usize i, usize j) const noexcept(utils::ndebug_defined());

  [[nodiscard]] inline static std::pair<usize, usize> transpose_coords(
      usize row, usize col) noexcept(utils::ndebug_defined());
  inline void bound_check_column(usize col) const;
  inline void check_for_overflow_on_add(usize row, usize col, N n) const;
  inline void check_for_overflow_on_subtract(usize row, usize col, N n) const;

  [[nodiscard]] inline N internal_get(usize row, usize col, usize block_size, std::mutex* mtx) const
      noexcept(utils::ndebug_defined());
  inline void internal_set(usize row, usize col, N n,
                           std::mutex* mtx) noexcept(utils::ndebug_defined());

  [[nodiscard]] constexpr double internal_get_fraction_of_missed_updates(
      std::mutex* mtx) const noexcept;
};
}  // namespace modle

#include "../../contacts_impl.hpp"  // IWYU pragma: export
// IWYU pragma: no_include "../../contacts_impl.hpp"
