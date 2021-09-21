// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

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
#include <fmt/format.h>            // for format, FMT_STRING
#include <fmt/ostream.h>           // for formatbuf<>::int_type
#include <spdlog/spdlog.h>         // for warn

#include <algorithm>    // for max, min, generate, fill
#include <cassert>      // for assert
#include <cstddef>      // for size_t
#include <cstdint>      // for int64_t, int32_t, uint8_t
#include <cstdio>       // for stderr
#include <memory>       // for unique_ptr, make_unique
#include <stdexcept>    // for runtime_error
#include <string>       // for string
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...
#include "modle/common/utils.hpp"                       // for ndebug_not_defined
#include "modle/contacts.hpp"                           // for ContactMatrix
#include "modle/hdf5.hpp"                               // for write_numbers, has_attribute, rea...

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
                                                   I chrom_end, I chrom_length) {
  ContactMatrix<int64_t> *null_matrix{nullptr};
  Cooler::write_or_append_cmatrix_to_file(null_matrix, chrom_name, chrom_start, chrom_end,
                                          chrom_length);
}

template <typename I1, typename I2>
void Cooler::write_or_append_cmatrix_to_file(const ContactMatrix<I1> &cmatrix,
                                             std::string_view chrom_name, I2 chrom_start,
                                             I2 chrom_end, I2 chrom_length) {
  Cooler::write_or_append_cmatrix_to_file(&cmatrix, chrom_name, chrom_start, chrom_end,
                                          chrom_length);
}

template <typename I1, typename I2>
void Cooler::write_or_append_cmatrix_to_file(const ContactMatrix<I1> *cmatrix,
                                             std::string_view chrom_name, I2 chrom_start_,
                                             I2 chrom_end_, I2 chrom_length_) {
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
    auto pxl_offset =
        static_cast<hsize_t>(hdf5::has_attribute(*this->_fp, "nbins", this->_root_path)
                                 ? hdf5::read_attribute_int(*this->_fp, "nbins", this->_root_path)
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
      if (cmatrix && cmatrix->get_n_of_missed_updates() != 0) {
        const auto &n = cmatrix->get_n_of_missed_updates();
        spdlog::warn(
            FMT_STRING(
                "Detected {} missed updates ({:.4f}% of the total number of contacts) for '{}'."),
            n, 100.0 * cmatrix->get_fraction_of_missed_updates(), chrom_name);
      }
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
      // will be 2 Mbp and chrom_end will be 8 Mbp. This for loop writes the index for all the bins
      // corresponding to the genomic region 0-2 Mbp. A more elegant solution would be to use
      // variable bin-size and write a single 2 Mbp bin.
      {
        const auto new_size = b->idx_bin1_offset_buff.size() +
                              ((chrom_start + this->_bin_size - 1) / this->_bin_size);
        b->idx_bin1_offset_buff.resize(new_size, this->_nnz);
      }
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

      pxl_offset += (chrom_start + this->_bin_size - 1) / this->_bin_size;
      if (cmatrix) {  // when cmatrix == nullptr we only write chrom/bins/indexes (no pixels)
        for (auto i = 0UL; i < cmatrix->ncols(); ++i) {  // Iterate over columns in the cmatrix
          // Write the first pixel that refers to a given bin1 to the index
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
                  throw std::runtime_error(fmt::format(
                      FMT_STRING("Cooler::write_or_append_cmatrix_to_file(): b1 > b2: b1={}; "
                                 "b2={}; offset={}; m={}\n"),
                      b->pixel_b1_idx_buff.back(), b->pixel_b2_idx_buff.back(), pxl_offset, m));
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
      {
        const auto new_size =
            b->idx_bin1_offset_buff.size() +
            ((static_cast<size_t>(chrom_length) - chrom_end + this->_bin_size - 1) /
             this->_bin_size);
        b->idx_bin1_offset_buff.resize(new_size, this->_nnz);
      }
      idx_bin1_offset_h5_foffset =
          hdf5::write_numbers(b->idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
      b->idx_bin1_offset_buff.clear();

      pxl_offset = this->_nbins;
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
  } catch ([[maybe_unused]] const H5::Exception &e) {
    throw std::runtime_error(hdf5::construct_error_stack());
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
std::unique_ptr<H5::DSetAccPropList> Cooler::generate_default_aprop([[maybe_unused]] T type,
                                                                    hsize_t chunk_size,
                                                                    hsize_t cache_size) {
  static_assert(std::is_same_v<T, H5::StrType> || std::is_same_v<T, H5::PredType>,
                "type should be of type H5::StrType or H5::PredType");
  H5::DSetAccPropList prop{};
  // https://support.hdfgroup.org/HDF5/doc/Advanced/Chunking/index.html
  // The w0 parameter affects how the library decides which chunk to evict when it needs room in the
  // cache. If w0 is set to 0, then the library will always evict the least recently used chunk in
  // cache. If w0 is set to 1, the library will always evict the least recently used chunk which has
  // been fully read or written, and if none have been fully read or written, it will evict the
  // least recently used chunk. If w0 is between 0 and 1, the behavior will be a blend of the two.

  // Therefore, if the application will access the same data more than once, w0 should be set closer
  // to 0, and if the application does not, w0 should be set closer to 1.
  const size_t default_multiplier{100};
  [[maybe_unused]] const auto rdcc_w0_streaming{0.99};
  [[maybe_unused]] const auto rdcc_w0_random{0.01};
  if constexpr (std::is_same_v<T, H5::StrType>) {
    prop.setChunkCache(default_multiplier * (cache_size / chunk_size), cache_size, rdcc_w0_random);
  } else {
    if (type == H5::PredType::NATIVE_INT64 ||
        type == H5::PredType::NATIVE_DOUBLE) {  // int64_t and double
      prop.setChunkCache(default_multiplier * (cache_size / chunk_size), cache_size,
                         rdcc_w0_streaming);
    } else if (type == H5::PredType::NATIVE_INT32) {  // int32_t
      prop.setChunkCache(default_multiplier * (cache_size / (chunk_size + chunk_size)), cache_size,
                         rdcc_w0_streaming);
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
  assert(bin_size_ != 0);  // NOLINT
  buff32.resize(buff_size);
  buff64.resize(buff_size);
  const auto chromosome = static_cast<int32_t>(chrom_);
  const auto length = static_cast<int64_t>(length_);
  const auto bin_size = static_cast<int64_t>(bin_size_);

  std::fill(buff32.begin(), buff32.end(), chromosome);

  int64_t start = 0;
  int64_t end = std::min(length, bin_size);
  hsize_t bins_processed = 0;
  const int64_t nbins = (length + bin_size - 1) / bin_size;
  const auto nchunks = static_cast<int64_t>((nbins + buff_size - 1) / buff_size);
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
  assert(buff64.back() == length);  // NOLINT

  return file_offset;
}

template <class N, class>
ContactMatrix<N> Cooler::cooler_to_cmatrix(std::string_view chrom_name, size_t nrows,
                                           std::pair<size_t, size_t> chrom_boundaries,
                                           bool try_common_chrom_prefixes,
                                           bool prefer_using_balanced_counts) {
  assert(this->_fp);                         // NOLINT
  assert(!this->_datasets.empty());          // NOLINT
  assert(!this->_idx_bin1_offset.empty());   // NOLINT
  assert(!this->_idx_chrom_offset.empty());  // NOLINT

  const auto chrom_idx = this->get_chrom_idx(chrom_name, try_common_chrom_prefixes);
  const auto chrom_size = static_cast<size_t>(this->get_chrom_sizes()[chrom_idx]);
  if (chrom_boundaries.second > chrom_size) {
    chrom_boundaries.second = chrom_size;
  }

  assert(chrom_boundaries.first < chrom_boundaries.second);  // NOLINT
  const auto bin_range = [&]() {
    const auto first_bin_chrom = static_cast<size_t>(this->_idx_chrom_offset[chrom_idx]);
    auto first_bin_range = first_bin_chrom + (chrom_boundaries.first / this->_bin_size);
    auto last_bin_range =
        first_bin_chrom + ((chrom_boundaries.second + this->_bin_size - 1) / this->_bin_size);

    return std::make_pair(first_bin_range, last_bin_range);
  }();

  const auto bin1_offset_idx = this->get_bin1_offset_idx_for_chrom(chrom_idx, chrom_boundaries);
  double pxl_count_scaling_factor{1.0};

  if (prefer_using_balanced_counts &&
      hdf5::has_dataset(*this->_fp, "bins/weight", this->_root_path)) {
    const auto &d = this->_datasets[BIN_WEIGHT];
    uint8_t cis_only;  // NOLINT
    try {
      hdf5::read_attribute(d, "cis_only", cis_only);
    } catch (const std::runtime_error &e) {
      if (absl::StrContains(e.what(), "cis_only")) {
        throw std::runtime_error(
            "File has a \"bins/weight\" dataset, but does not have an attribute named "
            "\"cis_only\". This most likely means that the file was generated by a very old "
            "version of Cooler. In order to proceed, you should rebalance the file using a "
            "recent version of Cooler, or in alternative remove the weight dataset (not "
            "recommended).");
      }
    }
    if (cis_only) {  // --cis-only balancing produces an array of sale factors
      std::vector<double> buff;
      hdf5::read_attribute(d, "scale", buff);
      pxl_count_scaling_factor = buff[chrom_idx];
    } else {  // standard or --trans-only balancing produces a single scale factor
      hdf5::read_attribute(d, "scale", pxl_count_scaling_factor);
    }
  }

  return cooler_to_cmatrix<N>(bin_range, bin1_offset_idx, nrows, pxl_count_scaling_factor,
                              prefer_using_balanced_counts);
}

template <class N, class>
ContactMatrix<N> Cooler::cooler_to_cmatrix(std::string_view chrom_name, size_t diagonal_width,
                                           size_t bin_size,
                                           std::pair<size_t, size_t> chrom_boundaries,
                                           bool try_common_chrom_prefixes,
                                           bool prefer_using_balanced_counts) {
  assert(this->_bin_size != 0);  // NOLINT
  if (bin_size != 0 && this->_bin_size != bin_size) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "Unable to read a Cooler file with bin size {} in a contact matrix with bin size {}"),
        this->_bin_size, bin_size));
  }
  const auto nrows = (diagonal_width + this->_bin_size - 1) / this->_bin_size;
  return cooler_to_cmatrix<N>(chrom_name, nrows, chrom_boundaries, try_common_chrom_prefixes,
                              prefer_using_balanced_counts);
}

template <class N, class>
ContactMatrix<N> Cooler::cooler_to_cmatrix(std::pair<hsize_t, hsize_t> bin_range,
                                           const std::vector<int64_t> &bin1_offset_idx,
                                           size_t nrows, double scaling_factor,
                                           bool prefer_using_balanced_counts) {
  return this->cooler_to_cmatrix<N>(bin_range, absl::MakeConstSpan(bin1_offset_idx), nrows,
                                    scaling_factor, prefer_using_balanced_counts);
}

template <class N, class>
ContactMatrix<N> Cooler::cooler_to_cmatrix(std::pair<hsize_t, hsize_t> bin_range,
                                           absl::Span<const int64_t> bin1_offset_idx, size_t nrows,
                                           double bias_scaling_factor,
                                           bool prefer_using_balanced_counts) {
  if (this->_datasets.empty()) {
    this->open_default_datasets();
  }

  const auto [first_bin, last_bin] = bin_range;
  assert(first_bin < last_bin);                            // NOLINT
  assert(last_bin <= first_bin + bin1_offset_idx.size());  // NOLINT
  ContactMatrix<N> cmatrix(nrows, last_bin - first_bin);
  std::vector<int64_t> bin1_BUFF(nrows);
  std::vector<int64_t> bin2_BUFF(nrows);
  std::vector<int64_t> count_BUFF(nrows);
  std::vector<double> bin_weights;

  const auto &d = this->_datasets;
  if (prefer_using_balanced_counts &&
      hdf5::has_dataset(*this->_fp, "bins/weight", this->_root_path)) {
    bin_weights.resize(last_bin - first_bin + 1);
    (void)hdf5::read_numbers(d[BIN_WEIGHT], bin_weights, static_cast<hsize_t>(first_bin));
  }

  for (size_t i = 1; i < bin1_offset_idx.size(); ++i) {
    const auto file_offset = static_cast<hsize_t>(bin1_offset_idx[i - 1]);
    const auto buff_size =
        std::min(static_cast<size_t>(bin1_offset_idx[i] - bin1_offset_idx[i - 1]), nrows);
    if (buff_size == 0) {
      continue;
    }

    bin1_BUFF.resize(buff_size);
    bin2_BUFF.resize(buff_size);
    count_BUFF.resize(buff_size);

    assert(file_offset + buff_size <= this->_idx_bin1_offset.back());  // NOLINT
    (void)hdf5::read_numbers(d[PXL_B1], bin1_BUFF, file_offset);
    (void)hdf5::read_numbers(d[PXL_B2], bin2_BUFF, file_offset);
    (void)hdf5::read_numbers(d[PXL_COUNT], count_BUFF, file_offset);

    assert(bin1_BUFF.size() == buff_size);   // NOLINT
    assert(bin2_BUFF.size() == buff_size);   // NOLINT
    assert(count_BUFF.size() == buff_size);  // NOLINT

    for (auto j = 0UL; j < buff_size; ++j) {
      assert(count_BUFF[j] != 0);  // NOLINT
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      DISABLE_WARNING_SIGN_COMPARE
      DISABLE_WARNING_CONVERSION
      const auto bin1 = bin1_BUFF[j] - first_bin;
      const auto bin2 = bin2_BUFF[j] - first_bin;
      if (bin2 >= i + nrows - 1 || bin2 >= bin1_offset_idx.size() - 1) {
        break;
      }
      if (bin_weights.empty()) {
        cmatrix.set(bin2, bin1, static_cast<N>(count_BUFF[j]));
      } else {
        // According to Cooler documentations, NaN means that a bin has been excluded by the
        // matrix balancing procedure. In this case we set the count to 0
        if (std::isnan(bin_weights[bin1]) || std::isnan(bin_weights[bin2])) {
          continue;  // Same as setting count to 0;
        }
        const auto bin1_bias = bin_weights[bin1];
        const auto bin2_bias = bin_weights[bin2];
        // See https://github.com/robomics/modle/issues/36 and
        // https://github.com/open2c/cooler/issues/35
        const auto count =
            static_cast<double>(count_BUFF[j]) / (bin1_bias * bin2_bias) / bias_scaling_factor;
        cmatrix.set(bin2, bin1, static_cast<N>(std::round(count)));
        DISABLE_WARNING_POP
      }
    }
  }

  return cmatrix;
}

}  // namespace modle::cooler
