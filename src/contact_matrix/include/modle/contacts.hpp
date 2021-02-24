#pragma once

#include <absl/types/span.h>
#ifdef USE_XOSHIRO
#include <XoshiroCpp.hpp>
#endif

#include <boost/dynamic_bitset.hpp>
#include <cstdint>  // for uint_*t
#include <fstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace modle {

#ifdef USE_XOSHIRO
using PRNG = XoshiroCpp::Xoshiro256PlusPlus;
using seeder = XoshiroCpp::SplitMix64;
#else
using PRNG = std::mt19937_64;
using seeder = std::seed_seq;
#endif

template <typename I>
class ContactMatrix {
  static_assert(std::is_integral_v<I>,
                "ContactMatrix requires an integral type as template argument.");

 public:
  // Constructors
  inline ContactMatrix() = default;
  inline ContactMatrix(std::size_t nrows, std::size_t ncols, bool fill_with_random_numbers = false);

  // Counts getters and setters
  [[nodiscard]] inline I get(std::size_t row, std::size_t col) const;
  template <typename I2>
  inline void set(std::size_t row, std::size_t col, I2 n);
  template <typename I2>
  inline void add(std::size_t row, std::size_t col, I2 n);
  template <typename I2>
  inline void subtract(std::size_t row, std::size_t col, I2 n);
  inline void increment(std::size_t row, std::size_t col);
  inline void decrement(std::size_t row, std::size_t col);

  // Shape/statistics getters
  [[nodiscard]] inline std::size_t ncols() const;
  [[nodiscard]] inline std::size_t nrows() const;
  [[nodiscard]] inline std::size_t npixels() const;
  [[nodiscard]] inline std::size_t npixels_after_masking() const;
  [[nodiscard]] inline uint64_t get_n_of_missed_updates() const;
  [[nodiscard]] inline uint64_t get_tot_contacts() const;
  [[nodiscard]] inline double get_avg_contact_density() const;
  [[nodiscard]] inline uint64_t get_matrix_size_in_bytes() const;
  [[nodiscard]] inline double get_matrix_size_in_mb() const;

  // Debug
  inline void print(bool full = false) const;
  [[nodiscard]] inline std::vector<std::vector<I>> generate_symmetric_matrix() const;

  // Mask operations
  inline void generate_mask_for_bins_without_contacts(boost::dynamic_bitset<>& mask) const;
  [[nodiscard]] inline boost::dynamic_bitset<> generate_mask_for_bins_without_contacts() const;

  // Misc
  inline void clear_missed_updates_counter();
  inline void reset();
  [[nodiscard]] inline bool empty() const;
  [[nodiscard]] inline absl::Span<const I> get_raw_count_vector() const;
  inline void compute_row_wise_contact_histogram(std::vector<uint64_t>& buff) const;
  [[nodiscard]] inline std::vector<uint64_t> compute_row_wise_contact_histogram() const;
  inline void deplete_contacts(double depletion_multiplier = 1.0);
  inline void add_noise(double mean, double std, PRNG& rand_eng);
  inline void add_noise(std::size_t bin_size, double mean, double stddev, PRNG& rand_eng);

 private:
  uint64_t _nrows{0};
  uint64_t _ncols{0};
  std::vector<I> _contacts{};
  uint64_t _tot_contacts{0};
  uint64_t _updates_missed{0};

  [[nodiscard]] inline I& at(std::size_t i, std::size_t j);
  [[nodiscard]] inline const I& at(std::size_t i, std::size_t j) const;
  [[nodiscard]] inline static std::pair<std::size_t, std::size_t> transpose_coords(std::size_t row,
                                                                                   std::size_t col);
};

}  // namespace modle

#include "../../contacts_impl.hpp"
