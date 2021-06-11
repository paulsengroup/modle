#pragma once

// IWYU pragma: no_include <H5DaccProp.h>
// IWYU pragma: no_include <H5DataSet.h>
// IWYU pragma: no_include <H5DataSpace.h>
// IWYU pragma: no_include <H5DcreatProp.h>
// IWYU pragma: no_include <H5File.h>
// IWYU pragma: no_include <H5Group.h>
// IWYU pragma: no_include <H5PredType.h>
// IWYU pragma: no_include <H5StrType.h>
// IWYU pragma: no_include <H5public.h>
// IWYU pragma: no_include "modle/src/libio/cooler_impl.hpp"

#include <H5Cpp.h>                                // IWYU pragma: keep
#include <absl/types/span.h>                      // for Span
#include <readerwriterqueue/readerwriterqueue.h>  // for BlockingReadWriterQueue

#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for int64_t, uint_fast8_t, uint32_t, int32_t
#include <filesystem>   // for path
#include <memory>       // for unique_ptr, allocator
#include <string>       // for string
#include <string_view>  // for string_view
#include <utility>      // for pair
#include <vector>       // for vector

namespace modle {
template <typename I>
class ContactMatrix;
}

namespace modle::cooler {
class Cooler {
 public:
  enum IO_MODE : uint_fast8_t { READ_ONLY, WRITE_ONLY };
  enum Flavor {
    UNK = 0,
    AUTO = 1,
    COOL = 2,
    MCOOL = 3,
    SCOOL = 4  // For the time being, we do not support SCOOL
  };

  struct Pixel {
    size_t row;
    size_t col;
    size_t count;
    [[nodiscard]] bool operator==(const Pixel &other) const {
      return this->row == other.row && this->col == other.col && this->count == other.count;
    }
  };

 private:
  struct InternalBuffers {
    explicit InternalBuffers(size_t buff_size = 1024 * 1024 / sizeof(int64_t));  // 1 MB
    std::vector<int64_t> bin_pos_buff;
    std::vector<int32_t> bin_chrom_buff;
    std::vector<int64_t> pixel_b1_idx_buff;
    std::vector<int64_t> pixel_b2_idx_buff;
    std::vector<int64_t> pixel_count_buff;
    std::vector<int64_t> idx_bin1_offset_buff;
    std::vector<int64_t> idx_chrom_offset_buff;

    [[nodiscard]] size_t capacity() const;
  };

 public:
  static constexpr uint_fast8_t DEFAULT_COMPRESSION_LEVEL = 6;
  static constexpr size_t DEFAULT_HDF5_BUFFER_SIZE = 1024 * 1024ULL;      // 1MB
  static constexpr size_t DEFAULT_HDF5_CHUNK_SIZE = 1024 * 1024ULL;       // 1MB
  static constexpr size_t DEFAULT_HDF5_CACHE_SIZE = 16 * 1024 * 1024ULL;  // 16MB

  Cooler() = delete;
  explicit Cooler(std::filesystem::path path_to_file, IO_MODE mode = READ_ONLY, size_t bin_size = 0,
                  size_t max_str_length = 0, std::string_view assembly_name = "",
                  Flavor flavor = AUTO, bool validate = true,
                  uint_fast8_t compression_lvl = DEFAULT_COMPRESSION_LEVEL,
                  size_t chunk_size = DEFAULT_HDF5_CHUNK_SIZE,
                  size_t cache_size = DEFAULT_HDF5_CACHE_SIZE);
  explicit Cooler(const std::string &path_to_file, IO_MODE mode = READ_ONLY, size_t bin_size = 0,
                  size_t max_str_length = 0, std::string_view assembly_name = "",
                  Flavor flavor = AUTO, bool validate = true,
                  uint_fast8_t compression_lvl = DEFAULT_COMPRESSION_LEVEL,
                  size_t chunk_size = DEFAULT_HDF5_CHUNK_SIZE,
                  size_t cache_size = DEFAULT_HDF5_CACHE_SIZE);
  explicit Cooler(std::string_view path_to_file, IO_MODE mode = READ_ONLY, size_t bin_size = 0,
                  size_t max_str_length = 0, std::string_view assembly_name = "",
                  Flavor flavor = AUTO, bool validate = true,
                  uint_fast8_t compression_lvl = DEFAULT_COMPRESSION_LEVEL,
                  size_t chunk_size = DEFAULT_HDF5_CHUNK_SIZE,
                  size_t cache_size = DEFAULT_HDF5_CACHE_SIZE);

  ~Cooler();

  Cooler(const Cooler &) = delete;
  Cooler &operator=(const Cooler &) = delete;
  Cooler(Cooler &&) = default;
  Cooler &operator=(Cooler &&) = default;

  // Getters
  [[nodiscard]] bool is_read_only() const;
  [[nodiscard]] const std::filesystem::path &get_path() const;

  [[nodiscard]] size_t get_nchroms();

  void get_chrom_names(std::vector<std::string> &buff);
  [[nodiscard]] std::vector<std::string> get_chrom_names();

  template <typename I>
  inline void get_chrom_sizes(std::vector<I> &buff);
  [[nodiscard]] std::vector<int64_t> get_chrom_sizes();

  [[nodiscard]] std::vector<std::pair<std::string, size_t>> get_chroms();

  [[nodiscard]] bool is_cool() const;
  [[nodiscard]] bool is_mcool() const;
  [[nodiscard]] bool is_scool() const;

  [[nodiscard]] size_t get_bin_size() const;

  [[nodiscard]] bool has_contacts_for_chrom(std::string_view chrom_name,
                                            bool try_common_chrom_prefixes = false);
  [[nodiscard]] bool has_contacts_for_chrom(size_t chrom_idx) const;

  // Write to file
  void write_metadata();

  template <typename I1, typename I2>
  inline void write_or_append_cmatrix_to_file(const ContactMatrix<I1> &cmatrix,
                                              std::string_view chrom_name, I2 chrom_start,
                                              I2 chrom_end, I2 chrom_length, bool quiet = false);

  template <typename I1, typename I2>
  inline void write_or_append_cmatrix_to_file(const ContactMatrix<I1> *cmatrix,
                                              std::string_view chrom_name, I2 chrom_start,
                                              I2 chrom_end, I2 chrom_length, bool quiet = false);
  template <typename I1, typename I2>
  inline void write_or_append_cmatrices_to_file(const std::vector<ContactMatrix<I1>> &cmatrices,
                                                const std::vector<std::string> &chrom_names,
                                                const std::vector<I2> &chrom_starts,
                                                const std::vector<I2> &chrom_ends,
                                                const std::vector<I2> &chrom_sizes,
                                                bool quiet = false);
  template <typename I1, typename I2>
  inline void write_or_append_cmatrices_to_file(const std::vector<ContactMatrix<I1> *> &cmatrices,
                                                const std::vector<std::string> &chrom_names,
                                                const std::vector<I2> &chrom_starts,
                                                const std::vector<I2> &chrom_ends,
                                                const std::vector<I2> &chrom_sizes,
                                                bool quiet = false);
  template <typename I1, typename I2>
  inline void write_or_append_cmatrices_to_file(absl::Span<ContactMatrix<I1> *const> cmatrices,
                                                absl::Span<const std::string> chrom_names,
                                                absl::Span<const I2> chrom_starts,
                                                absl::Span<const I2> chrom_ends,
                                                absl::Span<const I2> chrom_sizes,
                                                bool quiet = false);
  // Read from file
  [[nodiscard]] ContactMatrix<uint32_t> cooler_to_cmatrix(
      std::string_view chrom_name, size_t nrows,
      std::pair<size_t, size_t> chrom_boundaries = {0, -1}, bool try_common_chrom_prefixes = true,
      bool prefer_using_balanced_counts = true);
  [[nodiscard]] ContactMatrix<uint32_t> cooler_to_cmatrix(
      std::string_view chrom_name, size_t diagonal_width, size_t bin_size,
      std::pair<size_t, size_t> chrom_boundaries = {0, -1}, bool try_common_chrom_prefixes = true,
      bool prefer_using_balanced_counts = true);

  size_t stream_contacts_for_chrom(moodycamel::BlockingReaderWriterQueue<Pixel> &queue,
                                   std::string_view chrom_name, size_t nrows,
                                   std::pair<size_t, size_t> chrom_boundaries = {0, -1},
                                   bool try_common_chrom_prefixes = true,
                                   bool prefer_using_balanced_counts = true);
  size_t stream_contacts_for_chrom(moodycamel::BlockingReaderWriterQueue<Pixel> &queue,
                                   std::string_view chrom_name, size_t diagonal_width,
                                   size_t bin_size,
                                   std::pair<size_t, size_t> chrom_boundaries = {0, -1},
                                   bool try_common_chrom_prefixes = true,
                                   bool prefer_using_balanced_counts = true);
  // Misc
  [[nodiscard]] static bool validate_file_format(H5::H5File &f, Flavor expected_flavor,
                                                 IO_MODE mode = READ_ONLY, size_t bin_size = 0,
                                                 bool throw_on_failure = true);

 private:
  H5::StrType STR_TYPE;
  H5::PredType INT64_TYPE{H5::PredType::NATIVE_INT64};
  H5::PredType INT32_TYPE{H5::PredType::NATIVE_INT32};
  std::filesystem::path _path_to_file;
  std::string _root_path{"/"};
  IO_MODE _mode;
  size_t _bin_size;
  std::string _assembly_name;
  Flavor _flavor;
  std::unique_ptr<H5::H5File> _fp{nullptr};
  std::vector<H5::Group> _groups{};
  std::vector<H5::DataSet> _datasets{};
  std::vector<hsize_t> _dataset_file_offsets{};
  std::unique_ptr<H5::DataSpace> _mem_space{nullptr};
  std::vector<std::pair<H5::DataSpace, hsize_t>> _fspaces{};  // fspaces + file offset

  uint8_t _compression_lvl;
  hsize_t _chunk_size;
  hsize_t _cache_size;

  std::unique_ptr<H5::DSetCreatPropList> _cprop_str{nullptr};
  std::unique_ptr<H5::DSetCreatPropList> _cprop_int32{nullptr};
  std::unique_ptr<H5::DSetCreatPropList> _cprop_int64{nullptr};

  std::unique_ptr<H5::DSetAccPropList> _aprop_str{nullptr};
  std::unique_ptr<H5::DSetAccPropList> _aprop_int32{nullptr};
  std::unique_ptr<H5::DSetAccPropList> _aprop_int64{nullptr};
  std::unique_ptr<H5::DSetAccPropList> _aprop_float64{nullptr};

  std::unique_ptr<InternalBuffers> _buff{nullptr};

  std::vector<int64_t> _idx_chrom_offset{};
  std::vector<int64_t> _idx_bin1_offset{};
  hsize_t _nbins{0};
  int64_t _nchroms{0};
  int64_t _nnz{0};

  enum Groups : uint_fast8_t { chrom = 0, BIN = 1, PXL = 2, IDX = 3 };
  enum Datasets : uint_fast8_t {
    chrom_LEN = 0,
    chrom_NAME = 1,
    BIN_CHROM = 2,
    BIN_START = 3,
    BIN_END = 4,
    BIN_WEIGHT = 5,
    PXL_B1 = 6,
    PXL_B2 = 7,
    PXL_COUNT = 8,
    IDX_BIN1 = 9,
    IDX_CHR = 10,
    DEFAULT_DATASETS_NR = 11
  };

  [[nodiscard]] static H5::StrType generate_default_str_type(size_t max_str_length);
  template <typename T1, typename T2>
  [[nodiscard]] inline static std::unique_ptr<H5::DSetCreatPropList> generate_default_cprop(
      hsize_t chunk_size, uint8_t compression_lvl, T1 type, T2 fill_value);
  template <typename T>
  [[nodiscard]] inline static std::unique_ptr<H5::DSetAccPropList> generate_default_aprop(
      T type, hsize_t chunk_size, hsize_t cache_size);

  [[nodiscard]] static std::unique_ptr<H5::H5File> open_file(const std::filesystem::path &path,
                                                             IO_MODE mode, size_t bin_size = 0,
                                                             size_t max_str_length = 0,
                                                             Flavor flavor = AUTO,
                                                             bool validate = true);
  [[nodiscard]] static std::vector<H5::Group> open_groups(H5::H5File &f,
                                                          bool create_if_not_exist = false,
                                                          size_t bin_size = 0);
  void init_default_datasets();
  void open_default_datasets();

  template <typename I1, typename I2, typename I3>
  [[nodiscard]] inline hsize_t write_bins(I1 chrom, I2 length, I3 bin_size,
                                          std::vector<int32_t> &buff32,
                                          std::vector<int64_t> &buff64, hsize_t file_offset,
                                          hsize_t buff_size = DEFAULT_HDF5_BUFFER_SIZE /
                                                              sizeof(int64_t));

  size_t read_chrom_offset_idx();
  size_t read_bin1_offset_idx();

  [[nodiscard]] absl::Span<const int64_t> get_bin1_offset_idx_for_chrom(
      size_t chrom_idx, std::pair<size_t, size_t> chrom_subrange = {0, -1});
  [[nodiscard]] absl::Span<const int64_t> get_bin1_offset_idx_for_chrom(
      std::string_view chrom_name, std::pair<size_t, size_t> chrom_subrange = {0, -1});

  [[nodiscard]] std::pair<int64_t, int64_t> read_chrom_pixels_boundaries(
      std::string_view chrom_name);
  [[nodiscard]] std::pair<int64_t, int64_t> read_chrom_pixels_boundaries(size_t chrom_idx);

  [[nodiscard]] ContactMatrix<uint32_t> cooler_to_cmatrix(std::pair<hsize_t, hsize_t> bin_range,
                                                          absl::Span<const int64_t> bin1_offset_idx,
                                                          size_t nrows,
                                                          double bias_scaling_factor = 1.0,
                                                          bool prefer_using_balanced_counts = true);
  [[nodiscard]] ContactMatrix<uint32_t> cooler_to_cmatrix(
      std::pair<hsize_t, hsize_t> bin_range, const std::vector<int64_t> &bin1_offset_idx,
      size_t nrows, double scaling_factor = 1.0, bool prefer_using_balanced_counts = true);

  [[nodiscard]] size_t stream_contacts_for_chrom(
      moodycamel::BlockingReaderWriterQueue<Pixel> &queue, std::pair<hsize_t, hsize_t> bin_range,
      absl::Span<const int64_t> bin1_offset_idx, size_t nrows, double bias_scaling_factor = 1.0,
      bool prefer_using_balanced_counts = true);
  [[nodiscard]] size_t stream_contacts_for_chrom(
      moodycamel::BlockingReaderWriterQueue<Pixel> &queue, std::pair<hsize_t, hsize_t> bin_range,
      const std::vector<int64_t> &bin1_offset_idx, size_t nrows, double scaling_factor = 1.0,
      bool prefer_using_balanced_counts = true);

  [[nodiscard]] size_t get_chrom_idx(std::string_view query_chrom_name,
                                     bool try_common_chrom_prefixes = false);

  [[nodiscard]] static std::string flavor_to_string(Flavor f);
  [[nodiscard]] static Flavor detect_file_flavor(H5::H5File &f);
  [[nodiscard]] static bool validate_cool_flavor(H5::H5File &f, size_t bin_size,
                                                 std::string_view root_path = "/",
                                                 bool throw_on_failure = true,
                                                 bool check_version = true);
  [[nodiscard]] static bool validate_multires_cool_flavor(H5::H5File &f, size_t bin_size,
                                                          std::string_view root_path = "/",
                                                          bool throw_on_failure = true);
};

}  // namespace modle::cooler

#include "../../cooler_impl.hpp"  // IWYU pragma: keep
