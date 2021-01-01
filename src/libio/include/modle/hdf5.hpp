#pragma once

#include <H5Cpp.h>

#include <cstdint>
#include <memory>
#include <string_view>
#include <vector>

#include "modle/contacts.hpp"

namespace modle::hdf5 {

/*inline constexpr uint8_t DEFAULT_COOLER_COMPRESSION_LVL = 6;
enum CoolerDatasets {
  CHR_LEN = 0,
  CHR_NAME = 1,
  BIN_CHROM = 2,
  BIN_START = 3,
  BIN_END = 4,
  PXL_B1 = 5,
  PXL_B2 = 6,
  PXL_COUNT = 7,
  IDX_BIN1 = 8,
  IDX_CHR = 9
}; */

[[nodiscard]] inline std::string construct_error_stack(H5::Exception *e = nullptr);

template <typename S>
[[nodiscard]] inline hsize_t write_str(const S &str, H5::DataSet &dataset,
                                       const H5::StrType &str_type, hsize_t file_offset);
template <typename CS>
[[nodiscard]] inline hsize_t write_strings(const CS &strings, H5::DataSet &dataset,
                                           const H5::StrType &str_type, hsize_t file_offset);

template <typename N>
[[nodiscard]] inline hsize_t write_number(N &num, H5::DataSet &dataset, hsize_t file_offset);

template <typename CN>
[[nodiscard]] inline hsize_t write_numbers(CN &numbers, H5::DataSet &dataset, hsize_t file_offset);

template <typename N>
[[nodiscard]] inline hsize_t read_number(H5::DataSet &dataset, N &buff, hsize_t file_offset);

template <typename CN>
[[nodiscard]] inline hsize_t read_numbers(H5::DataSet &dataset, CN &buff, hsize_t file_offset);

[[nodiscard]] inline hsize_t read_str(H5::DataSet &dataset, std::string &buff, hsize_t file_offset);

[[nodiscard]] inline hsize_t read_strings(H5::DataSet &dataset, std::vector<std::string> &buff,
                                          hsize_t file_offset);

[[nodiscard]] inline std::string read_str(H5::DataSet &dataset, hsize_t file_offset);

[[nodiscard]] inline std::vector<std::string> read_strings(H5::DataSet &dataset,
                                                           hsize_t file_offset);

template <typename T>
inline void read_attribute(H5::H5File &f, std::string_view attr_name, T &buff,
                           std::string_view path = "");
template <typename T>
inline void read_attribute(std::string_view path_to_file, std::string_view attr_name, T &buff,
                           std::string_view path = "");

[[nodiscard]] inline std::string read_attribute_str(H5::H5File &f, std::string_view attr_name,
                                                    std::string_view path = "");
[[nodiscard]] inline std::string read_attribute_str(std::string_view path_to_file,
                                                    std::string_view attr_name,
                                                    std::string_view path = "");

[[nodiscard]] inline int64_t read_attribute_int(H5::H5File &f, std::string_view attr_name,
                                                std::string_view path = "");
[[nodiscard]] inline int64_t read_attribute_int(std::string_view path_to_file,
                                                std::string_view attr_name,
                                                std::string_view path = "");

inline bool group_exists(H5::H5File &f, std::string_view name, std::string_view root_path = "/");

template <typename T>
inline bool dataset_exists(H5::H5File &f, std::string_view name, T type,
                           std::string_view root_path = "/");
/*
// BEGINNING OF OLD HEADER

// The followings are high-level functions used to read/write Cooler files
template <typename I>
inline void write_modle_cmatrix_to_cooler(const std::vector<ContactMatrix<I>> &cmatrices,
                                          const std::vector<std::string> &chr_names,
                                          const std::vector<uint64_t> &chr_starts,
                                          const std::vector<uint64_t> &chr_ends,
                                          const std::vector<uint64_t> &chr_sizes, uint64_t bin_size,
                                          std::string_view output_file, bool force_overwrite);
template <typename I>
inline void write_modle_cmatrices_to_cooler(const std::vector<const ContactMatrix<I> *> &cmatrices,
                                            const std::vector<std::string> &chr_names,
                                            const std::vector<uint64_t> &chr_starts,
                                            const std::vector<uint64_t> &chr_ends,
                                            const std::vector<uint64_t> &chr_sizes,
                                            uint64_t bin_size, std::string_view output_file,
                                            bool force_overwrite);

// IMPRTANT: for the time being, all cooler_to_cmatrix overload do not support reading portions of a
// cooler file. This means that given a chromosome CHR of length 100Mbp, where we have simulated
// loop extrusion between 50 and 75Mb, if we write to disk the resultant contact matrix using one of
// the write_modle_cmatri* functions, reading the cooler file back in memory will produce a larger
// contact matrix. This happens because the new in-memory contact matrix will represent bins from
// 0-100Mbp. The counts will still be accurate, but the two cmatrices are not identical at the bit
// level.
[[nodiscard]] inline ContactMatrix<uint32_t> cooler_to_cmatrix(std::string_view path_to_file,
                                                               std::string_view chr_name,
                                                               std::size_t diagonal_width,
                                                               std::size_t bin_size,
                                                               bool try_common_chr_prefixes = true);

[[nodiscard]] inline ContactMatrix<uint32_t> cooler_to_cmatrix(std::string_view path_to_file,
                                                               std::string_view chr_name,
                                                               std::size_t nrows,
                                                               bool try_common_chr_prefixes = true);

[[nodiscard]] inline ContactMatrix<uint32_t> cooler_to_cmatrix(
    H5::H5File &f, int64_t bin_offset, const std::vector<int64_t> &bin1_offset_idx,
    std::size_t diagonal_width, std::size_t bin_size);

[[nodiscard]] inline ContactMatrix<uint32_t> cooler_to_cmatrix(
    H5::H5File &f, int64_t bin_offset, const std::vector<int64_t> &bin1_offset_idx,
    std::size_t nrows);

// These functions are used to WRITE Cooler files

inline hsize_t write_bins(H5::H5File &f, int32_t chrom, int64_t length, int64_t bin_size,
                          std::vector<int32_t> &buff32, std::vector<int64_t> &buff64,
                          hsize_t file_offset, hsize_t BUFF_SIZE = ONE_MB / sizeof(int64_t));

inline void write_metadata(H5::H5File &f, int32_t bin_size, std::string_view assembly_name = "");

// TODO: Change this to not use the header
H5::EnumType write_chroms(H5::H5File &f,
                          const typename std::vector<ContactMatrix<int64_t>::Header> &headers,
                          std::string_view path_to_chrom_sizes);

// These functions are used to READ Cooler file

// This function also validates Cooler files
[[nodiscard]] inline H5::H5File open_for_reading(std::string_view path_to_file,
                                                 bool validate = true);

// rdcc_w0 == 1 means that HDF5lib will evict from cache the least recently
// used chunk which has been fully read or written);
[[nodiscard]] inline std::vector<H5::DataSet> open_cooler_datasets(
    H5::H5File &f, hsize_t CACHE_INT64_NBYTES = 64, hsize_t CACHE_INT32_NBYTES = 32,  // NOLINT
    hsize_t CACHE_STR_NBYTES = 10, double CACHE_RDCC_W0 = 1.0);                       // NOLINT

template <typename I>
inline hsize_t read_int(H5::H5File &f, std::string_view dataset_name, I &BUFF, hsize_t file_offset);
inline hsize_t read_str(H5::H5File &f, std::string_view dataset_name, std::string &BUFF,
                        hsize_t file_offset);

template <typename I>
inline hsize_t read_vect_of_int(H5::H5File &f, std::string_view dataset_name, std::vector<I> &BUFF,
                                hsize_t file_offset);

// The first overload takes a buffer as parameter, and should be used whenever possible.
// The second overload will allocate and return a buffer.
inline std::size_t read_vect_of_str(H5::H5File &f, std::string_view dataset_name,
                                    std::vector<std::string> &BUFF, hsize_t file_offset);

[[nodiscard]] inline std::vector<std::string> read_vect_of_str(H5::H5File &f,
                                                               std::string_view dataset_name,
                                                               hsize_t file_offset);

// This function will read the offset for all the chromosomes
[[nodiscard]] inline std::vector<int64_t> read_chr_offset_idx(H5::H5File &f);

// This function will read the offset for a given
[[nodiscard]] inline std::size_t read_chr_offset_idx(H5::H5File &f, std::string_view chr_name);

// This function will read the index for all the chromosomes
[[nodiscard]] inline std::vector<int64_t> read_bin1_offset_idx(H5::H5File &f);

// These functions will read the bin index for a given chromosome
[[nodiscard]] inline std::vector<int64_t> read_bin1_offset_idx(H5::H5File &f, std::size_t chr_idx);
[[nodiscard]] inline std::vector<int64_t> read_bin1_offset_idx(H5::H5File &f,
                                                               std::string_view chr_name);

[[nodiscard]] inline std::pair<int64_t, int64_t> read_chrom_pixels_boundaries(H5::H5File &f,
                                                                              std::size_t chr_idx);

[[nodiscard]] inline std::pair<int64_t, int64_t> read_chrom_pixels_boundaries(
    H5::H5File &f, std::string_view chr_name);
/*
template <typename I>
inline void write_modle_cmatrices_to_cooler(const ContactMatrix<I> &cmatrix,
                                            std::string_view chr_name, uint64_t chr_start,
                                            uint64_t chr_end, uint64_t chr_length,
                                            std::string_view output_file, bool force_overwrite);

template <typename T>
inline void read_attribute(H5::H5File &f, std::string_view attr_name, T &buff,
                           std::string_view path = "");
template <typename T>
inline void read_attribute(std::string_view path_to_file, std::string_view attr_name, T &buff,
                           std::string_view path = "");

[[nodiscard]] inline std::string read_attribute_str(H5::H5File &f, std::string_view attr_name,
                                                    std::string_view path = "");
[[nodiscard]] inline std::string read_attribute_str(std::string_view path_to_file,
                                                    std::string_view attr_name,
                                                    std::string_view path = "");

[[nodiscard]] inline int64_t read_attribute_int(H5::H5File &f, std::string_view attr_name,
                                                std::string_view path = "");
[[nodiscard]] inline int64_t read_attribute_int(std::string_view path_to_file,
                                                std::string_view attr_name,
                                                std::string_view path = "");

// The following functions perform low-level actions on HDF5 files

template <typename DataType>
[[nodiscard]] inline H5::PredType getH5_type();

[[nodiscard]] inline H5::H5File init_cooler_file(std::string_view path_to_output_file,
                                                 bool force_overwrite = false);

[[nodiscard]] inline std::vector<H5::DataSet> init_cooler_datasets(
    H5::H5File &f, hsize_t MAX_STR_LENGTH, hsize_t CHUNK_DIMS = ONE_MB,
    uint8_t COMPRESSION_LEVEL = DEFAULT_COOLER_COMPRESSION_LVL);

template <typename S>
inline hsize_t write_str(const S &str, H5::H5File &f, std::string_view dataset_name,
                         hsize_t MAX_STR_LENGTH, hsize_t file_offset = 0,
                         hsize_t CHUNK_DIMS = ONE_MB / 2,
                         uint8_t COMPRESSION_LEVEL = DEFAULT_COOLER_COMPRESSION_LVL);
template <typename I>
inline hsize_t write_int(I num, H5::H5File &f, std::string_view dataset_name,
                         hsize_t file_offset = 0, hsize_t CHUNK_DIMS = ONE_MB / sizeof(I),
                         uint8_t COMPRESSION_LEVEL = DEFAULT_COOLER_COMPRESSION_LVL);

template <typename S>
inline hsize_t write_vect_of_str(std::vector<S> &data, H5::H5File &f, std::string_view dataset_name,
                                 hsize_t file_offset = 0, hsize_t CHUNK_DIMS = ONE_MB / 2,
                                 uint8_t COMPRESSION_LEVEL = DEFAULT_COOLER_COMPRESSION_LVL);

inline hsize_t write_vect_of_enums(std::vector<int32_t> &data, const H5::EnumType &ENUM,
                                   H5::H5File &f, std::string_view dataset_name,
                                   hsize_t file_offset = 0,
                                   hsize_t CHUNK_DIMS = ONE_MB / sizeof(int32_t),
                                   uint8_t COMPRESSION_LEVEL = DEFAULT_COOLER_COMPRESSION_LVL);

template <typename I>
inline hsize_t write_vect_of_int(std::vector<I> &data, H5::H5File &f, std::string_view dataset_name,
                                 hsize_t file_offset = 0, hsize_t CHUNK_DIMS = ONE_MB / sizeof(I),
                                 uint8_t COMPRESSION_LEVEL = DEFAULT_COOLER_COMPRESSION_LVL);
template <typename I>
inline hsize_t write_vect_of_int(std::vector<I> &data, H5::DataSet &dataset, hsize_t file_offset,
                                 std::shared_ptr<H5::DataSpace> mem_space);



[[nodiscard]] inline H5::EnumType init_enum_from_strs(const std::vector<std::string> &data,
                                                      int32_t offset = 0);
*/
}  // namespace modle::hdf5

#include "../../hdf5_impl.hpp"
