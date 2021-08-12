#pragma once

// IWYU pragma: private, include "modle/cooler.hpp"
// IWYU pragma: no_include <boost/core/checked_delete.hpp>
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <cpp-sort/sorter_facade.h>
// IWYU pragma: no_include <H5DaccProp.h>
// IWYU pragma: no_include <H5DataSet.h>
// IWYU pragma: no_include <H5DataSpace.h>
// IWYU pragma: no_include <H5DcreatProp.h>
// IWYU pragma: no_include <H5Exception.h>
// IWYU pragma: no_include <H5File.h>
// IWYU pragma: no_include <H5Group.h>
// IWYU pragma: no_include <H5PredType.h>
// IWYU pragma: no_include <H5Public.h>
// IWYU pragma: no_include "H5Public.h"
// IWYU pragma: no_include <H5Ppublic.h>
// IWYU pragma: no_include "H5Ppublic.h"
// IWYU pragma: no_include <H5public.h>
// IWYU pragma: no_include <H5StrType.h>
// IWYU pragma: no_include "hdf5_impl.hpp"
// IWYU pragma: no_include <ssrteam>

#include <H5Cpp.h>                 // IWYU pragma: keep
#include <absl/strings/str_cat.h>  // for StrCat
#include <absl/time/clock.h>       // for Now
#include <absl/time/time.h>        // for FormatDuration, operator-, Time
#include <absl/types/span.h>       // for MakeConstSpan, Span
#include <fmt/format.h>            // for FMT_STRING, print, format
#include <fmt/ostream.h>           // for formatbuf<>::int_type

#include <algorithm>    // for max, min, fill
#include <cassert>      // for assert
#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for int64_t, int32_t, uint64_t, uint8_t
#include <cstdio>       // for stderr
#include <memory>       // for unique_ptr, make_unique
#include <stdexcept>    // for runtime_error
#include <string>       // for string
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNI...
#include "modle/common/utils.hpp"                       // for throw_with_trace
#include "modle/contacts.hpp"                           // IWYU pragma: keep for ContactMatrix
#include "modle/hdf5.hpp"                               // for has_attribute, read_attribute_int

namespace modle::cooler {

template <typename I>
void Cooler::get_chrom_sizes(std::vector<I> &buff) {
  static_assert(std::is_integral_v<I>, "buff should be a vector of integers.");
  assert(this->_fp);
  const auto nchroms = this->get_nchroms();
  buff.resize(nchroms);
  if (!this->_datasets.empty()) {
    const auto &d = this->_datasets[chrom_LEN];
    (void)hdf5::read_numbers(d, buff, 0);
  } else {
    if (this->is_cool()) {
      auto d = this->_fp->openDataSet("/chroms/length", *this->_aprop_str);
      (void)hdf5::read_numbers(d, buff, 0);
    } else if (this->is_mcool()) {
      assert(this->_bin_size != 0);
      auto d = this->_fp->openDataSet(
          absl::StrCat("/resolutions/", this->_bin_size, "/chroms/length"), *this->_aprop_str);
      (void)hdf5::read_numbers(d, buff, 0);
    } else {
      assert(!this->is_scool());
      buff.clear();
    }
  }
}

template <typename I>
void Cooler::write_or_append_empty_cmatrix_to_file(std::string_view chrom_name, I chrom_start,
                                                   I chrom_end, I chrom_length, bool quiet) {
  ContactMatrix<int64_t> *null_matrix{nullptr};
  Cooler::write_or_append_cmatrix_to_file(null_matrix, chrom_name, chrom_start, chrom_end,
                                          chrom_length, quiet);
}

template <typename I1, typename I2>
void Cooler::write_or_append_cmatrix_to_file(const ContactMatrix<I1> &cmatrix,
                                             std::string_view chrom_name, I2 chrom_start,
                                             I2 chrom_end, I2 chrom_length, bool quiet) {
  Cooler::write_or_append_cmatrix_to_file(&cmatrix, chrom_name, chrom_start, chrom_end,
                                          chrom_length, quiet);
}

template <typename I1, typename I2>
void Cooler::write_or_append_cmatrix_to_file(const ContactMatrix<I1> *cmatrix,
                                             std::string_view chrom_name, I2 chrom_start_,
                                             I2 chrom_end_, I2 chrom_length_, bool quiet) {
  static_assert(std::is_integral_v<I1>, "I1 should be an integral type.");
  static_assert(std::is_integral_v<I2>, "I2 should be an integral type.");
  assert(chrom_start_ >= 0);   // NOLINT
  assert(chrom_end_ >= 0);     // NOLINT
  assert(chrom_length_ >= 0);  // NOLINT
  assert(this->_bin_size != 0);
  const auto chrom_start = static_cast<size_t>(chrom_start_);
  const auto chrom_end = static_cast<size_t>(chrom_end_);
  const auto chrom_length = static_cast<size_t>(chrom_length_);

  ++this->_nchroms;

  auto t0 = absl::Now();

  // Declare one buffer for each dataset
  if (!this->_buff) {
    this->_buff = std::make_unique<InternalBuffers>();
  }
  auto &b = this->_buff;

  // Declare aliases to file offsets for HDF5 datasets
  auto &chrom_name_h5_foffset = this->_dataset_file_offsets[chrom_NAME];
  auto &chrom_length_h5_foffset = this->_dataset_file_offsets[chrom_LEN];
  auto &bin_chrom_name_h5_foffset = this->_dataset_file_offsets[BIN_CHROM];
  auto &bin_start_h5_foffset = this->_dataset_file_offsets[BIN_START];
  auto &bin_end_h5_foffset = this->_dataset_file_offsets[BIN_END];
  auto &pixel_b1_id_h5_foffset = this->_dataset_file_offsets[PXL_B1];
  auto &pixel_b2_id_h5_foffset = this->_dataset_file_offsets[PXL_B2];
  auto &pixel_count_h5_foffset = this->_dataset_file_offsets[PXL_COUNT];
  auto &idx_bin1_offset_h5_foffset = this->_dataset_file_offsets[IDX_BIN1];
  auto &idx_chrom_offset_h5_foffset = this->_dataset_file_offsets[IDX_CHR];

  auto &d = this->_datasets;
  try {
    // Number of non-zero (values)
    this->_nnz = hdf5::has_attribute(*this->_fp, "nnz", this->_root_path)
                     ? hdf5::read_attribute_int(*this->_fp, "nnz", this->_root_path)
                     : 0;
    // Idx of the first non-zero pixel of the current chrom
    auto pxl_offset = static_cast<hsize_t>(
        hdf5::has_attribute(*this->_fp, "nbins", this->_root_path)
            ? hdf5::read_attribute_int(*this->_fp, "nbins", this->_root_path) - 1
            : 0);
    auto chrom_offset =
        static_cast<hsize_t>(hdf5::has_attribute(*this->_fp, "nchroms", this->_root_path)
                                 ? hdf5::read_attribute_int(*this->_fp, "nchroms", this->_root_path)
                                 : 0);

    auto write_pixels_to_file =
        [&]() {  // Lambda used to write pixel data to file. Mostly useful to reduce code bloat
          pixel_b1_id_h5_foffset =
              hdf5::write_numbers(b->pixel_b1_idx_buff, d[PXL_B1], pixel_b1_id_h5_foffset);
          pixel_b2_id_h5_foffset =
              hdf5::write_numbers(b->pixel_b2_idx_buff, d[PXL_B2], pixel_b2_id_h5_foffset);
          pixel_count_h5_foffset =
              hdf5::write_numbers(b->pixel_count_buff, d[PXL_COUNT], pixel_count_h5_foffset);

          b->pixel_b1_idx_buff.clear();
          b->pixel_b2_idx_buff.clear();
          b->pixel_count_buff.clear();
        };

    for (auto chrom_idx = 0UL; chrom_offset + chrom_idx < static_cast<size_t>(this->_nchroms);
         ++chrom_idx) {
      // Declare several aliases/variables to improve code readability in later sections
      if (!quiet) {
        fmt::print(stderr, FMT_STRING("Writing contacts for '{}' ({:.2f} Mbp)..."), chrom_name,
                   static_cast<double>(chrom_length) / 1.0e6);  // NOLINT
        if (cmatrix && cmatrix->get_n_of_missed_updates() != 0) {
          const auto &n = cmatrix->get_n_of_missed_updates();
          fmt::print(
              stderr, FMT_STRING(" WARNING: Detected {} missed updates ({:.4f}% of the total)."), n,
              static_cast<double>(100U * n) / static_cast<double>(cmatrix->get_tot_contacts()));
        }
      }
      const auto t1 = absl::Now();

      // Write chrom name and size
      chrom_name_h5_foffset =
          hdf5::write_str(chrom_name, d[chrom_NAME], STR_TYPE, chrom_name_h5_foffset);
      chrom_length_h5_foffset =
          hdf5::write_number(chrom_length, d[chrom_LEN], chrom_length_h5_foffset);
      // Add idx of the first bin belonging to the chromosome that is being processed
      b->idx_chrom_offset_buff.emplace_back(this->_nbins);

      // Write all fixed-size bins for the current chromosome. Maybe in the future we can switch to
      // use variable-size bins
      this->_nbins = this->write_bins(chrom_offset + chrom_idx, chrom_length, this->_bin_size,
                                      b->bin_chrom_buff, b->bin_pos_buff, this->_nbins);
      // Update file offsets of bin_* datasets
      bin_chrom_name_h5_foffset = this->_nbins;
      bin_start_h5_foffset = this->_nbins;
      bin_end_h5_foffset = this->_nbins;

      // In case we are simulating a subset of a chromosome (i.e. end - start != chrom size), write
      // the index for all the bins corresponding to genomic coordinates before the start position
      // Example: suppose we are writing contacts for a chromosome "C" of size 10 Mbp. Suppose we
      // also know that there are no contacts in the first and last 2 Mbps. In this case chrom_start
      // will be 2 Mbp and chrom_end will be 8 Mbp. This for loop write the index for all the bins
      // corresponding to the genomic region 0-2 Mbp. A more elegant solution would be to use
      // variable bin-size and write a single 2 Mbp bin.
      b->idx_bin1_offset_buff.resize(
          b->idx_bin1_offset_buff.size() + (chrom_start / this->_bin_size), this->_nnz);
      idx_bin1_offset_h5_foffset =
          hdf5::write_numbers(b->idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
      b->idx_bin1_offset_buff.clear();

      if (!cmatrix || cmatrix->ncols() == 0) {
        DISABLE_WARNING_PUSH
        DISABLE_WARNING_SIGN_CONVERSION
        b->idx_bin1_offset_buff.resize((chrom_length + this->_bin_size - 1) / this->_bin_size,
                                       this->_nnz);
        DISABLE_WARNING_POP
        idx_bin1_offset_h5_foffset =
            hdf5::write_numbers(b->idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
        b->idx_bin1_offset_buff.clear();
      }

      pxl_offset += chrom_start / this->_bin_size;
      if (cmatrix) {  // when cmatrix == nullptr we only write chrom/bins/indexes (no pixels)
        for (auto i = 0UL; i < cmatrix->ncols(); ++i) {  // Iterate over columns in the cmatrix
          // Write first pixel that refers to a given bin1 to the index
          b->idx_bin1_offset_buff.push_back(this->_nnz);

          // Iterate over rows of the cmatrix. The first condition serves the purpose to avoid
          // reading data from regions that are full of zeros by design (because cmatrix only stores
          // contacts for a certain width along the diagonal). The second condition makes sure we
          // are not reading past the array end
          for (auto j = i; j < i + cmatrix->nrows() && j < cmatrix->ncols(); ++j) {
            if (const auto m = cmatrix->get(i, j); m != 0) {  // Only write non-zero pixels
              if constexpr (utils::ndebug_not_defined()) {
                // Make sure we are always reading from the upper-triangle of the underlying square
                // contact matrix
                if (pxl_offset + i > pxl_offset + j) {
                  utils::throw_with_trace(std::runtime_error(fmt::format(
                      FMT_STRING("Cooler::write_or_append_cmatrix_to_file(): b1 > b2: b1={}; "
                                 "b2={}; offset={}; m={}\n"),
                      b->pixel_b1_idx_buff.back(), b->pixel_b2_idx_buff.back(), pxl_offset, m)));
                }
              }
              b->pixel_b1_idx_buff.push_back(static_cast<int64_t>(pxl_offset + i));
              b->pixel_b2_idx_buff.push_back(static_cast<int64_t>(pxl_offset + j));
              b->pixel_count_buff.push_back(m);
              ++this->_nnz;

              if (b->pixel_b1_idx_buff.size() ==
                  b->capacity()) {  // Write pixels to disk when buffer is full
                write_pixels_to_file();
              }
            }
          }

          if (b->idx_bin1_offset_buff.size() ==
              b->capacity()) {  // Write bin idx when buffer is full
            idx_bin1_offset_h5_foffset = hdf5::write_numbers(b->idx_bin1_offset_buff, d[IDX_BIN1],
                                                             idx_bin1_offset_h5_foffset);
            b->idx_bin1_offset_buff.clear();
          }
        }
      }

      // In case we are simulating a subset of a chromosome (i.e. end - start != chrom size), write
      // the index for all the bins corresponding to genomic coordinates after the end position. See
      // previous comment for an example
      b->idx_bin1_offset_buff.resize(
          b->idx_bin1_offset_buff.size() +
              ((static_cast<size_t>(chrom_length) - chrom_end) / this->_bin_size) +
              ((static_cast<size_t>(chrom_length) - chrom_end) % this->_bin_size != 0),
          this->_nnz);
      idx_bin1_offset_h5_foffset =
          hdf5::write_numbers(b->idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
      b->idx_bin1_offset_buff.clear();

      pxl_offset = this->_nbins;
      if (!quiet) {
        fmt::print(stderr, FMT_STRING(" DONE in {}!\n"), absl::FormatDuration(absl::Now() - t1));
      }
    }

    // Write all non-empty buffers to disk
    if (!b->pixel_b1_idx_buff.empty()) {
      write_pixels_to_file();
    }

    // Writing the chrom_index only at the end should be ok, given that most of the time we are
    // processing 20-100 chrom
    idx_chrom_offset_h5_foffset = hdf5::write_number(b->idx_chrom_offset_buff.back(), d[IDX_CHR],
                                                     idx_chrom_offset_h5_foffset);
    idx_bin1_offset_h5_foffset =
        hdf5::write_numbers(b->idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);

    auto &f = *this->_fp;
    hdf5::write_or_create_attribute(f, "nchroms", this->_nchroms);
    const auto nbins = static_cast<int64_t>(this->_nbins);
    hdf5::write_or_create_attribute(f, "nbins", nbins);
    hdf5::write_or_create_attribute(f, "nnz", this->_nnz);
    // this->_fp->flush(H5F_SCOPE_GLOBAL); // This is probably unnecessary
  } catch (const H5::Exception &err) {
    throw std::runtime_error(hdf5::construct_error_stack());
  }

  if (!quiet) {
    fmt::print(stderr, "All contacts have been written to file {}. Saved {:.2f}M pixels in {}.\n",
               this->_path_to_file, static_cast<double>(this->_nnz) / 1.0e6,
               absl::FormatDuration(absl::Now() - t0));
  }
}

template <typename T1, typename T2>
std::unique_ptr<H5::DSetCreatPropList> Cooler::generate_default_cprop(hsize_t chunk_size,
                                                                      uint8_t compression_lvl,
                                                                      T1 type, T2 fill_value) {
  static_assert(
      std::is_same_v<T2, int64_t> || std::is_same_v<T2, int32_t> ||
          std::is_constructible_v<H5std_string, T2>,
      "fill_value should have one of the following types: int64_t, int32_t or std::string.");
  static_assert(
      (std::is_constructible_v<H5std_string, T2> && std::is_same_v<T1, H5::StrType>) ||
          std::is_same_v<T1, H5::PredType>,
      "Incompatible data type for variables type and fill_value: if T2 is string "
      "constructible, then T1 must be H5::StrType, else T2 is integral and type is H5::PredType");
  (void)fill_value;
  (void)type;

  H5::DSetCreatPropList prop{};
  prop.setChunk(1, &chunk_size);
  prop.setDeflate(compression_lvl);
  if constexpr (!std::is_constructible_v<H5std_string, T2>) {
    prop.setFillValue(type, &fill_value);
  }

  return std::make_unique<H5::DSetCreatPropList>(prop);
}

template <typename T>
std::unique_ptr<H5::DSetAccPropList> Cooler::generate_default_aprop(T type, hsize_t chunk_size,
                                                                    hsize_t cache_size) {
  static_assert(std::is_same_v<T, H5::StrType> || std::is_same_v<T, H5::PredType>,
                "type should be of type H5::StrType or H5::PredType");
  H5::DSetAccPropList prop{};
  // https://support.hdfgroup.org/HDF5/doc/Advanced/Chunking/index.html
  constexpr size_t default_multiplier{100};
  constexpr double rdcc_w0{0.99};
  if constexpr (std::is_same_v<T, H5::StrType>) {
    prop.setChunkCache(default_multiplier * (cache_size / chunk_size), cache_size, 0.01);
    (void)type;
    (void)rdcc_w0;
  } else {
    if (type == H5::PredType::NATIVE_INT64) {  // int64_t
      prop.setChunkCache(default_multiplier * (cache_size / chunk_size), cache_size, rdcc_w0);
    } else if (type == H5::PredType::NATIVE_INT32) {  // int32_t
      prop.setChunkCache(default_multiplier * (cache_size / (chunk_size + chunk_size)), cache_size,
                         rdcc_w0);
    } else if (type == H5::PredType::NATIVE_DOUBLE) {
      prop.setChunkCache(default_multiplier * (cache_size / chunk_size), cache_size, rdcc_w0);
    } else {
      throw std::runtime_error(
          "Cooler::generate_default_aprop(), type should have type H5::StrType or be one of "
          "H5::PredType:NATIVE_INT, H5::PredType::NATIVE_INT64");
    }
  }
  return std::make_unique<H5::DSetAccPropList>(prop);
}

template <typename I1, typename I2, typename I3>
hsize_t Cooler::write_bins(I1 chrom_, I2 length_, I3 bin_size_, std::vector<int32_t> &buff32,
                           std::vector<int64_t> &buff64, hsize_t file_offset, hsize_t buff_size) {
  assert(bin_size_ != 0);
  buff32.resize(buff_size);
  buff64.resize(buff_size);
  const auto chromosome = static_cast<int32_t>(chrom_);
  const auto length = static_cast<int64_t>(length_);
  const auto bin_size = static_cast<int64_t>(bin_size_);

  std::fill(buff32.begin(), buff32.end(), chromosome);

  int64_t start = 0;
  int64_t end = std::min(length, bin_size);
  hsize_t bins_processed = 0;
  const int64_t nbins = (length / bin_size) + (length % bin_size != 0);
  const int64_t nchunks =
      (nbins / static_cast<int64_t>(buff_size)) + (nbins % static_cast<int64_t>(buff_size) != 0);
  auto &d = this->_datasets;

  for (auto i = 0; i < nchunks; ++i) {
    const auto chunk_size = std::min(buff_size, static_cast<hsize_t>(nbins) - bins_processed);
    if (chunk_size != buff_size) {
      buff32.resize(chunk_size);
      buff64.resize(chunk_size);
    }
    (void)hdf5::write_numbers(buff32, d[BIN_CHROM], file_offset);

    std::generate(buff64.begin(), buff64.end(), [&]() { return (start += bin_size) - bin_size; });
    (void)hdf5::write_numbers(buff64, d[BIN_START], file_offset);

    std::generate(buff64.begin(), buff64.end(), [&]() { return (end += bin_size) - bin_size; });
    if (chunk_size != buff_size) {
      buff64.back() = length;
    }
    file_offset = hdf5::write_numbers(buff64, d[BIN_END], file_offset);
    bins_processed += chunk_size;
  }
  assert(buff64.back() == length);

  return file_offset;
}

}  // namespace modle::cooler
