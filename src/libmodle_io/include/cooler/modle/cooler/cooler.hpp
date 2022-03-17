// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

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
#include <absl/types/variant.h>                   // for variant, monostate
#include <readerwriterqueue/readerwriterqueue.h>  // for BlockingReaderWriterQueue

#include <boost/filesystem/path.hpp>  // for path
#include <memory>                     // for unique_ptr, allocator
#include <string>                     // for string
#include <string_view>                // for string_view
#include <utility>                    // for pair
#include <vector>                     // for vector

#include "modle/common/common.hpp"  // for i64, i32, u8f, u32
#include "modle/contacts.hpp"       // for ContactMatrix

namespace modle {
template <class N>
class ContactMatrix;
}

namespace modle::cooler {

template <class N = contacts_t>
class Cooler {
  static_assert(std::is_arithmetic_v<N>);

 public:
  static constexpr bool IS_FP = std::is_floating_point_v<N>;
  using value_type = std::conditional_t<IS_FP, double, i32>;
  enum class IO_MODE : u8f {
    READ_ONLY,
    WRITE_ONLY
  };  // TODO Replace WRITE_ONLY with TRUNCATE and APPEND?
  enum class FLAVOR : u8f {
    UNK = 0,
    AUTO = 1,
    COOL = 2,
    MCOOL = 3,
    SCOOL = 4  // For the time being, we do not support SCOOL
  };

 private:
  struct InternalBuffers;

  H5::StrType STR_TYPE;
  boost::filesystem::path _path_to_file;
  std::string _root_path{"/"};
  IO_MODE _mode;
  usize _bin_size;
  std::string _assembly_name;
  FLAVOR _flavor;
  std::unique_ptr<H5::H5File> _fp{nullptr};
  std::vector<H5::Group> _groups{};
  std::vector<H5::DataSet> _datasets{};
  std::vector<hsize_t> _dataset_file_offsets{};
  std::unique_ptr<H5::DataSpace> _mem_space{nullptr};
  std::vector<std::pair<H5::DataSpace, hsize_t>> _fspaces{};  // fspaces + file offset

  u8 _compression_lvl;
  hsize_t _chunk_size;
  hsize_t _cache_size;

  std::unique_ptr<H5::DSetCreatPropList> _cprop_str{nullptr};
  std::unique_ptr<H5::DSetCreatPropList> _cprop_int32{nullptr};
  std::unique_ptr<H5::DSetCreatPropList> _cprop_int64{nullptr};
  std::unique_ptr<H5::DSetCreatPropList> _cprop_float64{nullptr};

  std::unique_ptr<H5::DSetAccPropList> _aprop_str{nullptr};
  std::unique_ptr<H5::DSetAccPropList> _aprop_int32{nullptr};
  std::unique_ptr<H5::DSetAccPropList> _aprop_int64{nullptr};
  std::unique_ptr<H5::DSetAccPropList> _aprop_float64{nullptr};

  std::unique_ptr<InternalBuffers> _buff{nullptr};

  std::vector<i64> _idx_chrom_offset{};
  std::vector<i64> _idx_bin1_offset{};
  hsize_t _nbins{0};
  i64 _nchroms{0};
  i64 _nnz{0};
  using sum_t = typename std::conditional<IS_FP, double, i64>::type;
  sum_t _sum{0};

  enum Groups : u8f { chrom = 0, BIN = 1, PXL = 2, IDX = 3 };
  enum Datasets : u8f {
    CHROM_LEN = 0,
    CHROM_NAME = 1,
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

 public:
  struct Pixel {
    usize row;
    usize col;
    N count;
    [[nodiscard]] bool operator==(const Pixel &other) const {
      return this->row == other.row && this->col == other.col && this->count == other.count;
    }
  };

 private:
  struct InternalBuffers {
    using M = Cooler::value_type;
    inline explicit InternalBuffers(usize buff_size = 1024 * 1024 / sizeof(M));  // 1 MB
    std::vector<i64> bin_pos_buff;
    std::vector<i32> bin_chrom_buff;
    std::vector<i64> pixel_b1_idx_buff;
    std::vector<i64> pixel_b2_idx_buff;
    std::vector<M> pixel_count_buff;
    std::vector<i64> idx_bin1_offset_buff;
    std::vector<i64> idx_chrom_offset_buff;

    [[nodiscard]] inline usize capacity() const;
  };

 public:
  static constexpr u8f DEFAULT_COMPRESSION_LEVEL = 6;
  static constexpr usize DEFAULT_HDF5_BUFFER_SIZE = 1024 * 1024ULL;      // 1MB
  static constexpr usize DEFAULT_HDF5_CHUNK_SIZE = 1024 * 1024ULL;       // 1MB
  static constexpr usize DEFAULT_HDF5_CACHE_SIZE = 16 * 1024 * 1024ULL;  // 16MB

  Cooler() = delete;
  inline explicit Cooler(boost::filesystem::path path_to_file, IO_MODE mode = IO_MODE::READ_ONLY,
                         usize bin_size = 0, usize max_str_length = 0,
                         std::string_view assembly_name = "", FLAVOR flavor = FLAVOR::AUTO,
                         bool validate = true, u8f compression_lvl = DEFAULT_COMPRESSION_LEVEL,
                         usize chunk_size = DEFAULT_HDF5_CHUNK_SIZE,
                         usize cache_size = DEFAULT_HDF5_CACHE_SIZE);

  inline ~Cooler();

  Cooler(const Cooler &) = delete;
  Cooler &operator=(const Cooler &) = delete;
#if defined(__clang__) && __clang_major__ < 8
  Cooler(Cooler &&) = default;
  Cooler &operator=(Cooler &&) = default;
#else
  Cooler(Cooler &&) noexcept = default;
  Cooler &operator=(Cooler &&) noexcept = default;
#endif

  // Getters
  [[nodiscard]] constexpr bool is_read_only() const noexcept;
  [[nodiscard]] inline const boost::filesystem::path &get_path() const;

  [[nodiscard]] inline usize get_nchroms();

  inline void get_chrom_names(std::vector<std::string> &buff);
  [[nodiscard]] inline std::vector<std::string> get_chrom_names();

  template <class I>
  inline void get_chrom_sizes(std::vector<I> &buff);
  [[nodiscard]] inline std::vector<i64> get_chrom_sizes();

  [[nodiscard]] inline std::vector<std::pair<std::string, usize>> get_chroms();

  [[nodiscard]] constexpr bool is_cool() const noexcept;
  [[nodiscard]] constexpr bool is_mcool() const noexcept;
  [[nodiscard]] constexpr bool is_scool() const noexcept;

  [[nodiscard]] constexpr usize get_bin_size() const noexcept;

  [[nodiscard]] inline bool has_contacts_for_chrom(std::string_view chrom_name,
                                                   bool try_common_chrom_prefixes = false);
  [[nodiscard]] inline bool has_contacts_for_chrom(usize chrom_idx) const;

  // Write to file
  inline void write_metadata();
  inline void write_metadata_attribute(std::string_view metadata_str);

  template <class I>
  inline void write_or_append_empty_cmatrix_to_file(std::string_view chrom_name, I chrom_start,
                                                    I chrom_end, I chrom_length);

  template <class M, class I,
            class = std::enable_if_t<std::is_arithmetic_v<N> && std::is_integral_v<I>>>
  inline void write_or_append_cmatrix_to_file(const ContactMatrix<M> &cmatrix,
                                              std::string_view chrom_name, I chrom_start,
                                              I chrom_end, I chrom_length);

  template <class M, class I,
            class = std::enable_if_t<std::is_arithmetic_v<N> && std::is_integral_v<I>>>
  inline void write_or_append_cmatrix_to_file(const ContactMatrix<M> *cmatrix,
                                              std::string_view chrom_name, I chrom_start,
                                              I chrom_end, I chrom_length);
  // Read from file
  [[nodiscard]] inline ContactMatrix<N> cooler_to_cmatrix(
      std::string_view chrom_name, usize nrows, std::pair<usize, usize> chrom_boundaries = {0, -1},
      bool try_common_chrom_prefixes = true, bool prefer_using_balanced_counts = true);
  [[nodiscard]] inline ContactMatrix<N> cooler_to_cmatrix(
      std::string_view chrom_name, usize diagonal_width, usize bin_size,
      std::pair<usize, usize> chrom_boundaries = {0, -1}, bool try_common_chrom_prefixes = true,
      bool prefer_using_balanced_counts = true);

  inline usize stream_contacts_for_chrom(moodycamel::BlockingReaderWriterQueue<Pixel> &queue,
                                         std::string_view chrom_name, usize nrows,
                                         std::pair<usize, usize> chrom_boundaries = {0, -1},
                                         bool try_common_chrom_prefixes = true,
                                         bool prefer_using_balanced_counts = true);
  inline usize stream_contacts_for_chrom(moodycamel::BlockingReaderWriterQueue<Pixel> &queue,
                                         std::string_view chrom_name, usize diagonal_width,
                                         usize bin_size,
                                         std::pair<usize, usize> chrom_boundaries = {0, -1},
                                         bool try_common_chrom_prefixes = true,
                                         bool prefer_using_balanced_counts = true);
  // Misc
  [[nodiscard]] static inline bool validate_file_format(H5::H5File &f, FLAVOR expected_flavor,
                                                        IO_MODE mode = IO_MODE::READ_ONLY,
                                                        usize bin_size = 0,
                                                        bool throw_on_failure = true);

 private:
  [[nodiscard]] static inline H5::StrType generate_default_str_type(usize max_str_length);
  template <class T1, class T2>
  [[nodiscard]] static inline std::unique_ptr<H5::DSetCreatPropList> generate_default_cprop(
      hsize_t chunk_size, u8 compression_lvl, T1 type, T2 fill_value);
  template <class T>
  [[nodiscard]] static inline std::unique_ptr<H5::DSetAccPropList> generate_default_aprop(
      T type, hsize_t chunk_size, hsize_t cache_size);

  [[nodiscard]] static inline std::unique_ptr<H5::H5File> open_file(
      const boost::filesystem::path &path, IO_MODE mode, usize bin_size = 0,
      usize max_str_length = 0, FLAVOR flavor = FLAVOR::AUTO, bool validate = true);
  [[nodiscard]] static inline std::vector<H5::Group> open_groups(H5::H5File &f,
                                                                 bool create_if_not_exist = false,
                                                                 usize bin_size = 0);
  inline void init_default_datasets();
  inline void open_default_datasets();

  template <class I1, class I2, class I3>
  [[nodiscard]] inline hsize_t write_bins(I1 chrom, I2 length, I3 bin_size,
                                          std::vector<i32> &buff32, std::vector<i64> &buff64,
                                          hsize_t file_offset,
                                          hsize_t buff_size = DEFAULT_HDF5_BUFFER_SIZE /
                                                              sizeof(i64));

  inline usize read_chrom_offset_idx();
  inline usize read_bin1_offset_idx();

  [[nodiscard]] inline absl::Span<const i64> get_bin1_offset_idx_for_chrom(
      usize chrom_idx, std::pair<usize, usize> chrom_subrange = {0, -1});
  [[nodiscard]] inline absl::Span<const i64> get_bin1_offset_idx_for_chrom(
      std::string_view chrom_name, std::pair<usize, usize> chrom_subrange = {0, -1});

  [[nodiscard]] inline std::pair<i64, i64> read_chrom_pixels_boundaries(
      std::string_view chrom_name);
  [[nodiscard]] inline std::pair<i64, i64> read_chrom_pixels_boundaries(usize chrom_idx);

  [[nodiscard]] inline ContactMatrix<N> cooler_to_cmatrix(std::pair<hsize_t, hsize_t> bin_range,
                                                          absl::Span<const i64> bin1_offset_idx,
                                                          usize nrows,
                                                          double bias_scaling_factor = 1.0,
                                                          bool prefer_using_balanced_counts = true);
  [[nodiscard]] inline ContactMatrix<N> cooler_to_cmatrix(std::pair<hsize_t, hsize_t> bin_range,
                                                          const std::vector<i64> &bin1_offset_idx,
                                                          usize nrows, double scaling_factor = 1.0,
                                                          bool prefer_using_balanced_counts = true);

  [[nodiscard]] inline usize stream_contacts_for_chrom(
      moodycamel::BlockingReaderWriterQueue<Pixel> &queue, std::pair<hsize_t, hsize_t> bin_range,
      absl::Span<const i64> bin1_offset_idx, usize nrows, double bias_scaling_factor = 1.0,
      bool prefer_using_balanced_counts = true);
  [[nodiscard]] inline usize stream_contacts_for_chrom(
      moodycamel::BlockingReaderWriterQueue<Pixel> &queue, std::pair<hsize_t, hsize_t> bin_range,
      const std::vector<i64> &bin1_offset_idx, usize nrows, double scaling_factor = 1.0,
      bool prefer_using_balanced_counts = true);

  [[nodiscard]] inline usize get_chrom_idx(std::string_view query_chrom_name,
                                           bool try_common_chrom_prefixes = false);

  [[nodiscard]] static inline std::string flavor_to_string(FLAVOR f);
  [[nodiscard]] static inline FLAVOR detect_file_flavor(H5::H5File &f);
  [[nodiscard]] static inline bool validate_cool_flavor(H5::H5File &f, usize bin_size,
                                                        std::string_view root_path = "/",
                                                        bool throw_on_failure = true,
                                                        bool check_version = true);
  [[nodiscard]] static inline bool validate_multires_cool_flavor(H5::H5File &f, usize bin_size,
                                                                 std::string_view root_path = "/",
                                                                 bool throw_on_failure = true);
};

}  // namespace modle::cooler

#include "../../../../cooler_impl.hpp"        // IWYU pragma: export
#include "../../../../cooler_read_impl.hpp"   // IWYU pragma: export
#include "../../../../cooler_write_impl.hpp"  // IWYU pragma: export
