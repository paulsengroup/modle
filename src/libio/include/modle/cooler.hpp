#pragma once

#include <H5Cpp.h>
#include <absl/types/span.h>

#include <filesystem>
#include <memory>
#include <mutex>
#include <string_view>
#include <utility>
#include <vector>

#include "modle/contacts.hpp"

namespace modle::cooler {
class Cooler {
 private:
  [[nodiscard]] inline static H5::StrType generate_default_str_type(std::size_t max_str_length);

 public:
  const H5::StrType STR_TYPE;
  const H5::PredType INT64_TYPE{H5::PredType::NATIVE_INT64};  // NOLINT
  const H5::PredType INT32_TYPE{H5::PredType::NATIVE_INT32};  // NOLINT
  enum IO_MODE : uint8_t { READ_ONLY, WRITE_ONLY };
  enum Flavor {
    UNK = 0,
    AUTO = 1,
    COOL = 2,
    MCOOL = 3,
    SCOOL = 4  // For the time being, we do not support SCOOL
  };
  enum Groups : uint8_t { CHR = 0, BIN = 1, PXL = 2, IDX = 3 };
  enum Datasets : uint8_t {
    CHR_LEN = 0,
    CHR_NAME = 1,
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

  Cooler() = delete;
  inline explicit Cooler(std::filesystem::path path_to_file, IO_MODE mode = READ_ONLY,
                         std::size_t bin_size = 0, std::size_t max_str_length = 64,  // NOLINT
                         std::string_view assembly_name = "", Flavor flavor = AUTO,
                         bool validate = true,
                         uint8_t compression_lvl = 6,                    // NOLINT
                         std::size_t chunk_size = 1024 * 1024ULL,        // 1 MB NOLINT
                         std::size_t cache_size = 16 * 1024 * 1024ULL);  // 16 MB NOLINT
  inline explicit Cooler(const std::string &path_to_file, IO_MODE mode = READ_ONLY,
                         std::size_t bin_size = 0, std::size_t max_str_length = 64,  // NOLINT
                         std::string_view assembly_name = "", Flavor flavor = AUTO,
                         bool validate = true,
                         uint8_t compression_lvl = 6,                    // NOLINT
                         std::size_t chunk_size = 1024 * 1024ULL,        // 1 MB NOLINT
                         std::size_t cache_size = 16 * 1024 * 1024ULL);  // 16 MB NOLINT
  inline explicit Cooler(std::string_view path_to_file, IO_MODE mode = READ_ONLY,
                         std::size_t bin_size = 0, std::size_t max_str_length = 64,  // NOLINT
                         std::string_view assembly_name = "", Flavor flavor = AUTO,
                         bool validate = true,
                         uint8_t compression_lvl = 6,                    // NOLINT
                         std::size_t chunk_size = 1024 * 1024ULL,        // 1 MB NOLINT
                         std::size_t cache_size = 16 * 1024 * 1024ULL);  // 16 MB NOLINT

  inline ~Cooler();

  [[nodiscard]] inline static Flavor detect_file_flavor(H5::H5File &f);

  [[nodiscard]] inline static bool validate_file_format(H5::H5File &f, Flavor expected_flavor,
                                                        IO_MODE mode = READ_ONLY,
                                                        std::size_t bin_size = 0,
                                                        bool throw_on_failure = true);

  [[nodiscard]] inline static bool validate_cool_flavor(H5::H5File &f, std::size_t bin_size,
                                                        std::string_view root_path = "/",
                                                        bool throw_on_failure = true,
                                                        bool check_version = true);

  [[nodiscard]] inline static bool validate_multires_cool_flavor(H5::H5File &f,
                                                                 std::size_t bin_size,
                                                                 std::string_view root_path = "/",
                                                                 bool throw_on_failure = true);
  template <typename I1, typename I2>
  inline void write_or_append_cmatrix_to_file(const ContactMatrix<I1> &cmatrix,
                                              std::string_view chr_name, I2 chr_start, I2 chr_end,
                                              I2 chr_length, bool quiet = false);

  template <typename I1, typename I2>
  inline void write_or_append_cmatrices_to_file(const std::vector<ContactMatrix<I1>> &cmatrices,
                                                const std::vector<std::string> &chr_names,
                                                const std::vector<I2> &chr_starts,
                                                const std::vector<I2> &chr_ends,
                                                const std::vector<I2> &chr_sizes,
                                                bool quiet = false);
  template <typename I1, typename I2>
  inline void write_or_append_cmatrices_to_file(const std::vector<ContactMatrix<I1> *> &cmatrices,
                                                const std::vector<std::string> &chr_names,
                                                const std::vector<I2> &chr_starts,
                                                const std::vector<I2> &chr_ends,
                                                const std::vector<I2> &chr_sizes,
                                                bool quiet = false);

  template <typename I1, typename I2>
  inline void write_or_append_cmatrices_to_file(absl::Span<ContactMatrix<I1> *const> cmatrices,
                                                absl::Span<const std::string> chr_names,
                                                absl::Span<const I2> chr_starts,
                                                absl::Span<const I2> chr_ends,
                                                absl::Span<const I2> chr_sizes, bool quiet = false);

  inline void write_metadata();
  [[nodiscard]] inline bool is_read_only() const;

  [[nodiscard]] inline ContactMatrix<uint32_t> cooler_to_cmatrix(
      std::string_view chr_name, std::size_t nrows, bool try_common_chr_prefixes = true,
      bool prefer_using_balanced_counts = true);

  [[nodiscard]] inline ContactMatrix<uint32_t> cooler_to_cmatrix(
      std::string_view chr_name, std::size_t diagonal_width, std::size_t bin_size,
      bool try_common_chr_prefixes = true, bool prefer_using_balanced_counts = true);
  [[nodiscard]] inline std::size_t get_nchroms();

  inline void get_chr_names(std::vector<std::string> &buff);
  [[nodiscard]] inline std::vector<std::string> get_chr_names();

  template <typename I>
  inline void get_chr_sizes(std::vector<I> &buff);
  [[nodiscard]] inline std::vector<int64_t> get_chr_sizes();

  [[nodiscard]] inline bool is_cool() const;
  [[nodiscard]] inline bool is_mcool() const;
  [[nodiscard]] inline bool is_scool() const;
  [[nodiscard]] inline std::size_t get_bin_size() const;

 private:
  std::filesystem::path _path_to_file;
  std::string _root_path{"/"};
  IO_MODE _mode;
  std::size_t _bin_size;
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

  std::mutex _mutex;  // This mutex protects _datasets and _idx_*_offset during initialization
  std::vector<int64_t> _idx_chrom_offset{};
  std::vector<int64_t> _idx_bin1_offset{};
  int64_t _nbins{0};
  int64_t _nchroms{0};
  int64_t _nnz{0};

  [[nodiscard]] inline static std::unique_ptr<H5::H5File> open_file(
      const std::filesystem::path &path, IO_MODE mode, std::size_t bin_size = 0,
      Flavor flavor = AUTO, bool validate = true);
  [[nodiscard]] inline static std::vector<H5::Group> open_groups(H5::H5File &f,
                                                                 bool create_if_not_exist = false,
                                                                 std::size_t bin_size = 0);
  [[nodiscard]] inline static std::string flavor_to_string(Flavor f);

  template <typename T1, typename T2>
  [[nodiscard]] inline static std::unique_ptr<H5::DSetCreatPropList> generate_default_cprop(
      hsize_t chunk_size, uint8_t compression_lvl, T1 type, T2 fill_value);
  template <typename T>
  [[nodiscard]] inline static std::unique_ptr<H5::DSetAccPropList> generate_default_aprop(
      T type, hsize_t chunk_size, hsize_t cache_size);
  inline void init_default_datasets();
  inline void open_cooler_datasets();
  [[nodiscard]] inline hsize_t write_bins(int32_t chrom, int64_t length, int64_t bin_size,
                                          std::vector<int32_t> &buff32,
                                          std::vector<int64_t> &buff64, hsize_t file_offset,
                                          hsize_t buff_size = 1024 * 1024ULL /  // NOLINT
                                                              sizeof(int64_t));
  [[nodiscard]] inline ContactMatrix<uint32_t> cooler_to_cmatrix(
      int64_t bin_offset, const std::vector<int64_t> &bin1_offset_idx, std::size_t nrows,
      double scaling_factor = 1.0, bool prefer_using_balanced_counts = true);

  [[nodiscard]] inline ContactMatrix<uint32_t> cooler_to_cmatrix(
      int64_t bin_offset, absl::Span<const int64_t> bin1_offset_idx, std::size_t nrows,
      double scaling_factor = 1.0, bool prefer_using_balanced_counts = true);

  [[nodiscard]] inline std::size_t get_chr_idx(std::string_view chr_name);

  inline std::size_t read_chr_offset_idx();

  inline std::size_t read_bin1_offset_idx();
  [[nodiscard]] inline absl::Span<const int64_t> get_bin1_offset_idx_for_chr(std::size_t chr_idx);
  [[nodiscard]] inline absl::Span<const int64_t> get_bin1_offset_idx_for_chr(
      std::string_view chr_name);
  [[nodiscard]] inline std::pair<int64_t, int64_t> read_chrom_pixels_boundaries(
      std::string_view chr_name);
  [[nodiscard]] inline std::pair<int64_t, int64_t> read_chrom_pixels_boundaries(
      std::size_t chr_idx);
};

}  // namespace modle::cooler

#include "../../cooler_impl.hpp"
