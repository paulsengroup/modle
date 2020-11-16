#pragma once

#include <boost/iostreams/filtering_stream.hpp>
#include <cstdint>
#include <fstream>
#include <random>
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
  [[nodiscard]] I get(uint64_t row, uint64_t col) const;
  void set(uint64_t row, uint64_t col, I n = 1);
  void increment(uint64_t row, uint64_t col, I n = 1);
  void decrement(uint64_t row, uint64_t col, I n = 1);
  void print(bool full = false) const;
  [[nodiscard]] std::vector<std::vector<I>> generate_symmetric_matrix() const;
  [[nodiscard]] uint64_t n_cols() const;
  [[nodiscard]] uint64_t n_rows() const;
  std::pair<uint32_t /* bytes_in */, uint32_t /* bytes_out */> write_full_matrix_to_tsv(
      const std::string& path_to_file, int bzip2_block_size = 9) const;
  std::pair<uint32_t /* bytes_in */, uint32_t /* bytes_out */> write_to_tsv(
      const std::string& path_to_file, std::string_view header = "",
      int bzip2_block_size = 9) const;
  void clear_missed_updates_counter();
  const std::vector<I>& get_raw_count_vector() const;
  [[nodiscard]] static Header parse_header(std::string_view path_to_file);
  [[nodiscard]] static Header parse_header(std::string_view path_to_file,
                                           boost::iostreams::filtering_istream& in,
                                           bool rewind_file = false);

 private:
  uint64_t _nrows;
  uint64_t _ncols;
  std::vector<I> _contacts;
  uint64_t _updates_missed{0};

  [[nodiscard]] I& at(uint64_t i, uint64_t j);
  [[nodiscard]] const I& at(uint64_t i, uint64_t j) const;
};

}  // namespace modle

#include "../../contacts_impl.hpp"
