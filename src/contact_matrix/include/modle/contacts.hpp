#pragma once

#include <boost/iostreams/filtering_stream.hpp>
#include <cstdint>  // for uint_*t
#include <fstream>
#include <random>  // for normal_distribution
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace modle {

template <typename I>
class ContactMatrix {
  static_assert(std::is_integral<I>::value,
                "ContactMatrix requires an integral type as template argument.");

 public:
  struct Header {
    std::string chr_name{};
    uint64_t bin_size{};
    std::size_t start{};
    std::size_t end{};
    std::size_t diagonal_width{};
    std::size_t ncols{};
    std::size_t nrows{};
    [[nodiscard]] std::string to_string() const;
    // #chr_name   bin_size   start   end   diagonal_width
    static constexpr uint8_t N_OF_EXPECTED_TOKENS = 5;
  };

  ContactMatrix(uint64_t nrows, uint64_t ncols, bool fill_with_random_numbers = false);
  ContactMatrix(const std::string& path_to_file, uint64_t nrows, uint64_t ncols, char sep = '\t');
  explicit ContactMatrix(const std::string& path_to_file,
                         std::normal_distribution<double>* noise_generator = nullptr,
                         uint64_t seed = 0, char sep = '\t');
  [[nodiscard]] I get(std::size_t row, std::size_t col) const;
  template <typename I2>
  void set(std::size_t row, std::size_t col, I2 n);
  template <typename I2>
  void increment(std::size_t row, std::size_t col, I2 n);
  template <typename I2>
  void decrement(std::size_t row, std::size_t col, I2 n);
  void increment(std::size_t row, std::size_t col);
  void decrement(std::size_t row, std::size_t col);
  void print(bool full = false) const;
  [[nodiscard]] std::vector<std::vector<I>> generate_symmetric_matrix() const;
  [[nodiscard]] std::size_t n_cols() const;
  [[nodiscard]] std::size_t n_rows() const;
  std::pair<uint32_t /* bytes_in */, uint32_t /* bytes_out */> write_full_matrix_to_tsv(
      const std::string& path_to_file, int bzip2_block_size = 9) const;
  std::pair<uint32_t /* bytes_in */, uint32_t /* bytes_out */> write_to_tsv(
      const std::string& path_to_file, std::string_view header = "",
      int bzip2_block_size = 9) const;
  void clear_missed_updates_counter();
  [[nodiscard]] const std::vector<I>& get_raw_count_vector() const;
  [[nodiscard]] uint64_t get_n_of_missed_updates() const;
  [[nodiscard]] uint64_t get_tot_contacts() const;
  [[nodiscard]] static Header parse_header(std::string_view path_to_file);
  [[nodiscard]] static Header parse_header(std::string_view path_to_file,
                                           boost::iostreams::filtering_istream& in,
                                           bool rewind_file = false);

 private:
  uint64_t _nrows;
  uint64_t _ncols;
  std::vector<I> _contacts;
  uint64_t _tot_contacts{0};
  uint64_t _updates_missed{0};

  [[nodiscard]] I& at(std::size_t i, std::size_t j);
  [[nodiscard]] const I& at(std::size_t i, std::size_t j) const;
};

}  // namespace modle

#include "../../contacts_impl.hpp"
