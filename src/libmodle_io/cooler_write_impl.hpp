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

#include <H5Cpp.h>            // IWYU pragma: keep
#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatTime, UTCTimeZone
#include <fmt/format.h>       // for format

#include <string>       // for string
#include <string_view>  // for string_view

#include "modle/hdf5/hdf5.hpp"  // for read_attribute, read_numbers, wri...

namespace modle::cooler {

template <class N>
void Cooler<N>::write_metadata() {
  if (this->is_read_only()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Caught attempt to write metadata to an HDF5 file that is open in "
                               "read-only mode. File name: {}"),
                    this->_path_to_file.string()));
  }
  assert(this->_bin_size != 0);  // NOLINT
  H5::DataSpace attr_space(H5S_SCALAR);
  i64 int_buff{};
  std::string str_buff{};
  std::string name{};

  try {
    name = "format";
    str_buff = "HDF5::Cooler";
    hdf5::write_or_create_attribute(*this->_fp, name, str_buff);

    name = "format-version";
    int_buff = 3;
    hdf5::write_or_create_attribute(*this->_fp, name, int_buff);

    name = "bin-type";
    str_buff = "fixed";
    hdf5::write_or_create_attribute(*this->_fp, name, str_buff);

    name = "bin-size";
    int_buff = static_cast<i64>(this->_bin_size);
    hdf5::write_or_create_attribute(*this->_fp, name, int_buff);

    name = "storage-mode";
    str_buff = "symmetric-upper";
    hdf5::write_or_create_attribute(*this->_fp, name, str_buff);

    if (!this->_assembly_name.empty()) {
      name = "assembly-name";
      str_buff = this->_assembly_name;
      hdf5::write_or_create_attribute(*this->_fp, name, str_buff);
    }

    name = "generated-by";
    str_buff = modle_version_long;
    hdf5::write_or_create_attribute(*this->_fp, name, str_buff);

    name = "creation-date";
    str_buff = absl::FormatTime(absl::Now(), absl::UTCTimeZone());
    hdf5::write_or_create_attribute(*this->_fp, name, str_buff);
  } catch ([[maybe_unused]] const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while writing metadata to file {}: "
                               "error while writing attribute '{}':\n{}"),
                    this->_path_to_file, name, hdf5::construct_error_stack()));
  }
}

template <class N>
void Cooler<N>::write_metadata_attribute(std::string_view metadata_str) {
  if (this->is_read_only()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Caught attempt to write metadata to an HDF5 file that is open in "
                               "read-only mode. File name: {}"),
                    this->_path_to_file.string()));
  }

  assert(this->_bin_size != 0);   // NOLINT
  assert(!metadata_str.empty());  // NOLINT
  H5::DataSpace attr_space(H5S_SCALAR);
  const std::string name = "metadata";
  const std::string buff{metadata_str};

  try {
    hdf5::write_or_create_attribute(*this->_fp, name, buff);
  } catch ([[maybe_unused]] const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while writing metadata to file {}: "
                               "error while writing attribute '{}':\n{}"),
                    this->_path_to_file, name, hdf5::construct_error_stack()));
  }
}

template <class N>
template <class I>
void Cooler<N>::write_or_append_empty_cmatrix_to_file(std::string_view chrom_name, I chrom_start,
                                                      I chrom_end, I chrom_length) {
  ContactMatrix<N> *null_matrix{nullptr};
  Cooler::write_or_append_cmatrix_to_file(null_matrix, chrom_name, chrom_start, chrom_end,
                                          chrom_length);
}

template <class N>
template <class M, class I, class>
void Cooler<N>::write_or_append_cmatrix_to_file(const ContactMatrix<M> &cmatrix,
                                                std::string_view chrom_name, I chrom_start,
                                                I chrom_end, I chrom_length) {
  Cooler::write_or_append_cmatrix_to_file(&cmatrix, chrom_name, chrom_start, chrom_end,
                                          chrom_length);
}

template <class N>
template <class M, class I, class>
void Cooler<N>::write_or_append_cmatrix_to_file(const ContactMatrix<M> *cmatrix,
                                                std::string_view chrom_name, I chrom_start_,
                                                I chrom_end_, I chrom_length_) {
  static_assert(std::is_floating_point_v<N> == std::is_floating_point_v<M>,
                "Cooler<> and ContactMatrix<> template arguments should both be integral or "
                "floating point types.");
  assert(chrom_start_ >= 0);   // NOLINT
  assert(chrom_end_ >= 0);     // NOLINT
  assert(chrom_length_ >= 0);  // NOLINT
  assert(this->_bin_size != 0);
  const auto chrom_start = static_cast<usize>(chrom_start_);
  const auto chrom_end = static_cast<usize>(chrom_end_);
  const auto chrom_length = static_cast<usize>(chrom_length_);

  ++this->_nchroms;

  assert(this->_buff);  // NOLINT
  auto &b = *this->_buff;

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
              hdf5::write_numbers(b.pixel_b1_idx_buff, d[PXL_B1], pixel_b1_id_h5_foffset);
          pixel_b2_id_h5_foffset =
              hdf5::write_numbers(b.pixel_b2_idx_buff, d[PXL_B2], pixel_b2_id_h5_foffset);
          pixel_count_h5_foffset =
              hdf5::write_numbers(b.pixel_count_buff, d[PXL_COUNT], pixel_count_h5_foffset);

          b.pixel_b1_idx_buff.clear();
          b.pixel_b2_idx_buff.clear();
          b.pixel_count_buff.clear();
        };

    for (usize chrom_idx = 0; chrom_offset + chrom_idx < static_cast<usize>(this->_nchroms);
         ++chrom_idx) {
      if (cmatrix && cmatrix->get_n_of_missed_updates() != 0) {
        const auto &n = cmatrix->get_n_of_missed_updates();
        spdlog::warn(
            FMT_STRING(
                "Detected {} missed updates ({:.4f}% of the total number of contacts) for '{}'."),
            n, 100.0 * cmatrix->unsafe_get_fraction_of_missed_updates(), chrom_name);
      }
      // Write chrom name and size
      chrom_name_h5_foffset =
          hdf5::write_str(chrom_name, d[chrom_NAME], STR_TYPE, chrom_name_h5_foffset);
      chrom_length_h5_foffset =
          hdf5::write_number(chrom_length, d[chrom_LEN], chrom_length_h5_foffset);
      // Add idx of the first bin belonging to the chromosome that is being processed
      b.idx_chrom_offset_buff.emplace_back(this->_nbins);

      // Write all fixed-size bins for the current chromosome. Maybe in the future we can switch to
      // use variable-size bins
      this->_nbins = this->write_bins(chrom_offset + chrom_idx, chrom_length, this->_bin_size,
                                      b.bin_chrom_buff, b.bin_pos_buff, this->_nbins);
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
        const auto new_size =
            b.idx_bin1_offset_buff.size() + ((chrom_start + this->_bin_size - 1) / this->_bin_size);
        b.idx_bin1_offset_buff.resize(new_size, this->_nnz);
      }
      idx_bin1_offset_h5_foffset =
          hdf5::write_numbers(b.idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
      b.idx_bin1_offset_buff.clear();

      if (!cmatrix || cmatrix->ncols() == 0) {
        DISABLE_WARNING_PUSH
        DISABLE_WARNING_SIGN_CONVERSION
        b.idx_bin1_offset_buff.resize((chrom_length + this->_bin_size - 1) / this->_bin_size,
                                      this->_nnz);
        DISABLE_WARNING_POP
        idx_bin1_offset_h5_foffset =
            hdf5::write_numbers(b.idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
        b.idx_bin1_offset_buff.clear();
      }

      pxl_offset += (chrom_start + this->_bin_size - 1) / this->_bin_size;
      if (cmatrix) {  // when cmatrix == nullptr we only write chrom/bins/indexes (no pixels)
        for (usize i = 0; i < cmatrix->ncols(); ++i) {  // Iterate over columns in the cmatrix
          // Write the first pixel that refers to a given bin1 to the index
          b.idx_bin1_offset_buff.push_back(this->_nnz);

          // Iterate over rows of the cmatrix. The first condition serves the purpose to avoid
          // reading data from regions that are full of zeros by design (because cmatrix only stores
          // contacts for a certain width along the diagonal). The second condition makes sure we
          // are not reading past the array end
          for (auto j = i; j < i + cmatrix->nrows() && j < cmatrix->ncols(); ++j) {
            const auto m = [&]() {
              DISABLE_WARNING_PUSH
              DISABLE_WARNING_USELESS_CAST
              if constexpr (std::is_floating_point_v<M> && !std::is_floating_point_v<value_type>) {
                return static_cast<value_type>(std::round(cmatrix->unsafe_get(i, j)));
              } else {
                return static_cast<value_type>(cmatrix->unsafe_get(i, j));
              }
              DISABLE_WARNING_POP
            }();
            if (m != value_type(0)) {  // Only write non-zero pixels
              if constexpr (utils::ndebug_not_defined()) {
                // Make sure we are always reading from the upper-triangle of the underlying square
                // contact matrix
                if (pxl_offset + i > pxl_offset + j) {
                  throw std::runtime_error(fmt::format(
                      FMT_STRING("Cooler::write_or_append_cmatrix_to_file(): b1 > b2: b1={}; "
                                 "b2={}; offset={}; m={}\n"),
                      b.pixel_b1_idx_buff.back(), b.pixel_b2_idx_buff.back(), pxl_offset, m));
                }
              }
              b.pixel_b1_idx_buff.push_back(static_cast<i64>(pxl_offset + i));
              b.pixel_b2_idx_buff.push_back(static_cast<i64>(pxl_offset + j));

              b.pixel_count_buff.push_back(m);

              ++this->_nnz;

              if (b.pixel_b1_idx_buff.size() ==
                  b.capacity()) {  // Write pixels to disk when buffer is full
                write_pixels_to_file();
              }
            }
          }

          if (b.idx_bin1_offset_buff.size() == b.capacity()) {  // Write bin idx when buffer is full
            idx_bin1_offset_h5_foffset = hdf5::write_numbers(b.idx_bin1_offset_buff, d[IDX_BIN1],
                                                             idx_bin1_offset_h5_foffset);
            b.idx_bin1_offset_buff.clear();
          }
        }
      }

      // In case we are simulating a subset of a chromosome (i.e. end - start != chrom size), write
      // the index for all the bins corresponding to genomic coordinates after the end position. See
      // previous comment for an example
      {
        const auto new_size =
            b.idx_bin1_offset_buff.size() +
            ((static_cast<usize>(chrom_length) - chrom_end + this->_bin_size - 1) /
             this->_bin_size);
        b.idx_bin1_offset_buff.resize(new_size, this->_nnz);
      }
      idx_bin1_offset_h5_foffset =
          hdf5::write_numbers(b.idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
      b.idx_bin1_offset_buff.clear();

      pxl_offset = this->_nbins;
    }

    // Write all non-empty buffers to disk
    if (!b.pixel_b1_idx_buff.empty()) {
      write_pixels_to_file();
    }

    // Writing the chrom_index only at the end should be ok, given that most of the time we are
    // processing 20-100 chrom
    idx_chrom_offset_h5_foffset =
        hdf5::write_number(b.idx_chrom_offset_buff.back(), d[IDX_CHR], idx_chrom_offset_h5_foffset);
    idx_bin1_offset_h5_foffset =
        hdf5::write_numbers(b.idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);

    auto &f = *this->_fp;
    hdf5::write_or_create_attribute(f, "nchroms", this->_nchroms);
    const auto nbins = static_cast<i64>(this->_nbins);
    hdf5::write_or_create_attribute(f, "nbins", nbins);
    hdf5::write_or_create_attribute(f, "nnz", this->_nnz);
    // this->_fp->flush(H5F_SCOPE_GLOBAL); // This is probably unnecessary
  } catch ([[maybe_unused]] const H5::Exception &e) {
    throw std::runtime_error(hdf5::construct_error_stack());
  }
}

template <class N>
template <class I1, class I2, class I3>
hsize_t Cooler<N>::write_bins(I1 chrom_, I2 length_, I3 bin_size_, std::vector<i32> &buff32,
                              std::vector<i64> &buff64, hsize_t file_offset, hsize_t buff_size) {
  assert(bin_size_ != 0);  // NOLINT
  buff32.resize(buff_size);
  buff64.resize(buff_size);
  const auto chromosome = static_cast<i32>(chrom_);
  const auto length = static_cast<i64>(length_);
  const auto bin_size = static_cast<i64>(bin_size_);

  std::fill(buff32.begin(), buff32.end(), chromosome);

  i64 start = 0;
  i64 end = std::min(length, bin_size);
  hsize_t bins_processed = 0;
  const i64 nbins = (length + bin_size - 1) / bin_size;
  const auto nchunks = (static_cast<hsize_t>(nbins) + buff_size - 1) / buff_size;
  auto &d = this->_datasets;

  for (usize i = 0; i < nchunks; ++i) {
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

}  // namespace modle::cooler
