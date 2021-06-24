#pragma once

// IWYU pragma: no_include "modle/src/contact_matrix/contacts_impl.hpp"

#include <absl/types/span.h>  // for Span

#include <atomic>                                   // for atomic
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cstddef>                                  // IWYU pragma: keep for size_t
#include <cstdint>                                  // for uint64_t
#include <mutex>                                    // for mutex
#include <utility>                                  // for pair
#include <vector>                                   // for vector

#include "modle/common/common.hpp"
#include "modle/common/utils.hpp"  // for ndebug_defined

namespace modle {

template <typename I = uint32_t>
class ContactMatrix {
  static_assert(std::is_integral_v<I>,
                "ContactMatrix requires an integral type as template argument.");

 public:
  // Constructors
  inline ContactMatrix() = default;
  inline ContactMatrix(ContactMatrix<I>&& other) noexcept = default;
  inline ContactMatrix(const ContactMatrix<I>& other);
  inline ContactMatrix(size_t nrows, size_t ncols, bool fill_with_random_numbers = false);
  template <typename I2, typename = std::enable_if_t<std::is_integral_v<I2>>>
  inline ContactMatrix(I2 length, I2 diagonal_width, I2 bin_size,
                       bool fill_with_random_numbers = false);
  inline ContactMatrix(absl::Span<const I> contacts, size_t nrows, size_t ncols,
                       size_t tot_contacts = 0, size_t updates_missed = 0);
  inline ~ContactMatrix() = default;

  // Operators
  [[nodiscard]] inline ContactMatrix& operator=(const ContactMatrix& other);
  [[nodiscard]] inline ContactMatrix& operator=(ContactMatrix&& other) noexcept = default;

  // Counts getters and setters
  [[nodiscard]] inline I get(size_t row, size_t col) const noexcept(utils::ndebug_defined());
  template <typename I2>
  inline void set(size_t row, size_t col, I2 n) noexcept(utils::ndebug_defined());
  template <typename I2>
  inline void add(size_t row, size_t col, I2 n) noexcept(utils::ndebug_defined());
  template <typename I2>
  inline void add(absl::Span<std::pair<size_t /* rows */, size_t /* cols */>> pixels, I2 n,
                  size_t size_thresh = 256) noexcept(utils::ndebug_defined());  // NOLINT

  template <typename I2>
  inline void subtract(size_t row, size_t col, I2 n) noexcept(utils::ndebug_defined());
  inline void increment(size_t row, size_t col) noexcept(utils::ndebug_defined());
  inline void increment(absl::Span<std::pair<size_t /* rows */, size_t /* cols */>> pixels,
                        size_t size_thresh = 256) noexcept(utils::ndebug_defined());  // NOLINT
  inline void decrement(size_t row, size_t col) noexcept(utils::ndebug_defined());

  // Shape/statistics getters
  [[nodiscard]] inline constexpr size_t ncols() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr size_t nrows() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr size_t npixels() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline size_t npixels_after_masking() const;
  [[nodiscard]] inline constexpr size_t get_n_of_missed_updates() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr size_t get_tot_contacts() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline double get_avg_contact_density() const noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr size_t get_matrix_size_in_bytes() const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline constexpr double get_matrix_size_in_mb() const
      noexcept(utils::ndebug_defined());

  // Debug
  inline void print(bool full = false) const;
  [[nodiscard]] inline std::vector<std::vector<I>> generate_symmetric_matrix() const;

  // Mask operations
  inline void generate_mask_for_bins_without_contacts(boost::dynamic_bitset<>& mask) const;
  [[nodiscard]] inline boost::dynamic_bitset<> generate_mask_for_bins_without_contacts() const;

  // Misc
  inline void clear_missed_updates_counter();
  inline void reset();
  inline void resize(size_t nrows, size_t ncols);
  template <typename I2, typename = std::enable_if_t<std::is_integral_v<I2>>>
  inline void resize(I2 length, I2 diagonal_width, I2 bin_size);
  [[nodiscard]] inline bool empty() const;
  [[nodiscard]] inline absl::Span<const I> get_raw_count_vector() const;
  [[nodiscard]] inline std::vector<I>& get_raw_count_vector();
  inline void compute_row_wise_contact_histogram(std::vector<uint64_t>& buff) const;
  [[nodiscard]] inline std::vector<uint64_t> compute_row_wise_contact_histogram() const;
  inline void deplete_contacts(double depletion_multiplier = 1.0);

 private:
  uint64_t _nrows{0};
  uint64_t _ncols{0};
  std::vector<I> _contacts{};
  std::atomic<size_t> _tot_contacts{0};
  std::atomic<size_t> _updates_missed{0};
  std::vector<std::mutex> _locks{};

  [[nodiscard]] inline I& at(size_t i, size_t j) noexcept(utils::ndebug_defined());
  [[nodiscard]] inline const I& at(size_t i, size_t j) const noexcept(utils::ndebug_defined());

  template <typename I2>
  inline void add_small_buff(absl::Span<std::pair<size_t /* rows */, size_t /* cols */>> pixels,
                             I2 n) noexcept(utils::ndebug_defined());
  template <typename I2>
  inline void add_large_buff(absl::Span<std::pair<size_t /* rows */, size_t /* cols */>> pixels,
                             I2 n) noexcept(utils::ndebug_defined());

  [[nodiscard]] inline static std::pair<size_t, size_t> transpose_coords(
      size_t row, size_t col) noexcept(utils::ndebug_defined());
  inline void bound_check_column(size_t col) const;
  template <typename I2>
  inline void check_for_overflow_on_add(size_t row, size_t col, I2 n) const;
  template <typename I2>
  inline void check_overflow_on_subtract(size_t row, size_t col, I2 n) const;
};

}  // namespace modle

#include "../../contacts_impl.hpp"  // IWYU pragma: keep
