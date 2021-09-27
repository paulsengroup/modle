// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span

#include <atomic>                                   // for atomic
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cstddef>                                  // for size_t
#include <cstdint>                                  // for uint64_t, uint32_t
#include <iostream>                                 // for cout
#include <mutex>                                    // for mutex
#include <type_traits>                              // for enable_if_t
#include <utility>                                  // for pair

#include "modle/common/common.hpp"
#include "modle/common/utils.hpp"  // for ndebug_defined

namespace modle {

template <class N = contacts_t>
class ContactMatrix {
  static_assert(std::is_arithmetic_v<N>,
                "ContactMatrix requires an integral type as template argument.");

 public:
  // Constructors
  inline ContactMatrix() = default;
#if defined(__clang__) && __clang_major__ < 7
  inline ContactMatrix(ContactMatrix<N>&& other) noexcept;
#else
  inline ContactMatrix(ContactMatrix<N>&& other) noexcept = default;
#endif
  inline ContactMatrix(const ContactMatrix<N>& other);
  inline ContactMatrix(size_t nrows, size_t ncols, bool fill_with_random_numbers = false,
                       uint64_t seed = 8336046165695760686);
  inline ContactMatrix(bp_t length, bp_t diagonal_width, bp_t bin_size,
                       bool fill_with_random_numbers = false);
  inline ContactMatrix(absl::Span<const N> contacts, size_t nrows, size_t ncols,
                       size_t tot_contacts = 0, size_t updates_missed = 0);
  inline ~ContactMatrix() = default;

  // Operators
  [[nodiscard]] inline ContactMatrix<N>& operator=(const ContactMatrix<N>& other);
#if defined(__clang__) && __clang_major__ < 7
  [[nodiscard]] inline ContactMatrix<N>& operator=(ContactMatrix<N>&& other) noexcept;
#else
  [[nodiscard]] inline ContactMatrix<N>& operator=(ContactMatrix<N>&& other) noexcept = default;
#endif

  // Counts getters and setters
  [[nodiscard]] inline N unsafe_get(size_t row, size_t col) const noexcept(utils::ndebug_defined());

  // block_size is required to be an odd number at the moment
  [[nodiscard]] inline N unsafe_get(size_t row, size_t col, size_t block_size) const
      noexcept(utils::ndebug_defined());
  inline void unsafe_set(size_t row, size_t col, N n) noexcept(utils::ndebug_defined());
  inline void set(size_t row, size_t col, N n) noexcept(utils::ndebug_defined());

  // The following four methods are thread-safe but operations are not atomic, as without locking it
  // is not possible to ensure that updates to pixels and total counters are atomic.
  // In order to avoid overflows, special care should be taken when N is unsigned and add/subtract
  // or increment/decrement are called from different threads
  inline void add(size_t row, size_t col, N n) noexcept(utils::ndebug_defined());
  inline void subtract(size_t row, size_t col, N) noexcept(utils::ndebug_defined());
  inline void increment(size_t row, size_t col) noexcept(utils::ndebug_defined());
  inline void decrement(size_t row, size_t col) noexcept(utils::ndebug_defined());

  // Shape/statistics getters
  [[nodiscard]] inline constexpr size_t ncols() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr size_t nrows() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr size_t npixels() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline size_t unsafe_npixels_after_masking() const;
  [[nodiscard]] inline constexpr size_t get_n_of_missed_updates() const noexcept;
  [[nodiscard]] inline constexpr double unsafe_get_fraction_of_missed_updates() const noexcept;
  [[nodiscard]] inline constexpr size_t get_tot_contacts() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline double get_avg_contact_density() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr size_t get_matrix_size_in_bytes() const
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
  inline void unsafe_resize(size_t nrows, size_t ncols);
  inline void unsafe_resize(bp_t length, bp_t diagonal_width, bp_t bin_size);
  [[nodiscard]] inline bool empty() const;
  [[nodiscard]] inline absl::Span<const N> get_raw_count_vector() const;
  [[nodiscard]] inline absl::Span<N> get_raw_count_vector();
  inline void unsafe_compute_row_wise_contact_histogram(std::vector<uint64_t>& buff) const;
  [[nodiscard]] inline std::vector<uint64_t> unsafe_compute_row_wise_contact_histogram() const;
  inline void unsafe_deplete_contacts(double depletion_multiplier = 1.0);

  using value_type = N;

 private:
  uint64_t _nrows{0};
  uint64_t _ncols{0};
  std::vector<N> _contacts{};
  std::vector<std::mutex> _mtxes{};  // Each mutex protects one row of bins
  std::atomic<int64_t> _tot_contacts{0};
  std::atomic<int64_t> _updates_missed{0};

  [[nodiscard]] inline N& at(size_t i, size_t j) noexcept(utils::ndebug_defined());
  [[nodiscard]] inline const N& at(size_t i, size_t j) const noexcept(utils::ndebug_defined());

  [[nodiscard]] inline static std::pair<size_t, size_t> transpose_coords(
      size_t row, size_t col) noexcept(utils::ndebug_defined());
  inline void bound_check_column(size_t col) const;
  inline void check_for_overflow_on_add(size_t row, size_t col, N n) const;
  inline void check_for_overflow_on_subtract(size_t row, size_t col, N n) const;

  [[nodiscard]] inline N internal_get(size_t row, size_t col, size_t block_size,
                                      std::mutex* mtx) const noexcept(utils::ndebug_defined());
  inline void internal_set(size_t row, size_t col, N n,
                           std::mutex* mtx) noexcept(utils::ndebug_defined());

  [[nodiscard]] constexpr double internal_get_fraction_of_missed_updates(
      std::mutex* mtx) const noexcept;
};
}  // namespace modle

#include "../../contacts_impl.hpp"  // IWYU pragma: export
// IWYU pragma: no_include "../../contacts_impl.hpp"
