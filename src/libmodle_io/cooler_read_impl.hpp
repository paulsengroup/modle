// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
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

#include <H5Cpp.h>                                // IWYU pragma: keep
#include <readerwriterqueue/readerwriterqueue.h>  // for BlockingReaderWriterQueue

#include <cassert>      // for assert
#include <string_view>  // for string_view
#include <tuple>        // for ignore
#include <utility>      // for make_pair

#include "modle/common/common.hpp"  // for i64, i32, u8, u32
#include "modle/hdf5/hdf5.hpp"      // for write_numbers, read_numbers, has_...

namespace modle::cooler {

template <class N>
ContactMatrixDense<N> Cooler<N>::cooler_to_cmatrix(std::string_view chrom_name, usize nrows,
                                                   std::pair<usize, usize> chrom_boundaries,
                                                   bool try_common_chrom_prefixes,
                                                   bool prefer_using_balanced_counts) {
  assert(this->_fp);
  assert(!this->_datasets.empty());
  assert(!this->_idx_bin1_offset.empty());
  assert(!this->_idx_chrom_offset.empty());

  const auto chrom_idx = this->get_chrom_idx(chrom_name, try_common_chrom_prefixes);
  const auto chrom_size = static_cast<usize>(this->get_chrom_sizes()[chrom_idx]);
  if (chrom_boundaries.second > chrom_size) {
    chrom_boundaries.second = chrom_size;
  }

  assert(chrom_boundaries.first < chrom_boundaries.second);
  const auto bin_range = [&]() {
    const auto first_bin_chrom = static_cast<usize>(this->_idx_chrom_offset[chrom_idx]);
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
    u8 cis_only;
    try {
      hdf5::read_attribute(d, "cis_only", cis_only);
    } catch (const std::runtime_error &e) {
      if (absl::StrContains(e.what(), "cis_only")) {
        throw std::runtime_error(
            "File has a \"bins/weight\" dataset, but dataset does not have an attribute named "
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

  return cooler_to_cmatrix(bin_range, bin1_offset_idx, nrows, pxl_count_scaling_factor,
                           prefer_using_balanced_counts);
}

template <class N>
ContactMatrixDense<N> Cooler<N>::cooler_to_cmatrix(std::string_view chrom_name,
                                                   usize diagonal_width, usize bin_size,
                                                   std::pair<usize, usize> chrom_boundaries,
                                                   bool try_common_chrom_prefixes,
                                                   bool prefer_using_balanced_counts) {
  assert(this->_bin_size != 0);
  if (bin_size != 0 && this->_bin_size != bin_size) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "Unable to read a Cooler file with bin size {} in a contact matrix with bin size {}"),
        this->_bin_size, bin_size));
  }
  const auto nrows = (diagonal_width + this->_bin_size - 1) / this->_bin_size;
  return cooler_to_cmatrix(chrom_name, nrows, chrom_boundaries, try_common_chrom_prefixes,
                           prefer_using_balanced_counts);
}

template <class N>
ContactMatrixDense<N> Cooler<N>::cooler_to_cmatrix(std::pair<hsize_t, hsize_t> bin_range,
                                                   const std::vector<i64> &bin1_offset_idx,
                                                   usize nrows, double scaling_factor,
                                                   bool prefer_using_balanced_counts) {
  return this->cooler_to_cmatrix(bin_range, absl::MakeConstSpan(bin1_offset_idx), nrows,
                                 scaling_factor, prefer_using_balanced_counts);
}

template <class N>
ContactMatrixDense<N> Cooler<N>::cooler_to_cmatrix(std::pair<hsize_t, hsize_t> bin_range,
                                                   absl::Span<const i64> bin1_offset_idx,
                                                   usize nrows, double bias_scaling_factor,
                                                   bool prefer_using_balanced_counts) {
  if (this->_datasets.empty()) {
    this->open_default_datasets();
  }

  const auto [first_bin, last_bin] = bin_range;
  assert(first_bin < last_bin);
  assert(last_bin <= first_bin + bin1_offset_idx.size());
  ContactMatrixDense<N> cmatrix(nrows, last_bin - first_bin);
  std::vector<i64> bin1_BUFF(nrows);
  std::vector<i64> bin2_BUFF(nrows);
  std::vector<N> count_BUFF(nrows);
  std::vector<double> bin_weights;

  const auto &d = this->_datasets;
  if (prefer_using_balanced_counts &&
      hdf5::has_dataset(*this->_fp, "bins/weight", this->_root_path)) {
    bin_weights.resize(last_bin - first_bin);
    std::ignore = hdf5::read_numbers(d[BIN_WEIGHT], bin_weights, static_cast<hsize_t>(first_bin));
  }

  for (usize i = 1; i < bin1_offset_idx.size(); ++i) {
    const auto file_offset = static_cast<hsize_t>(bin1_offset_idx[i - 1]);
    const auto buff_size =
        std::min(static_cast<usize>(bin1_offset_idx[i] - bin1_offset_idx[i - 1]), nrows);
    if (buff_size == 0) {
      continue;
    }

    bin1_BUFF.resize(buff_size);
    bin2_BUFF.resize(buff_size);
    count_BUFF.resize(buff_size);

    assert(static_cast<i64>(file_offset + buff_size) <= this->_idx_bin1_offset.back());
    std::ignore = hdf5::read_numbers(d[PXL_B1], bin1_BUFF, file_offset);
    std::ignore = hdf5::read_numbers(d[PXL_B2], bin2_BUFF, file_offset);
    std::ignore = hdf5::read_numbers(d[PXL_COUNT], count_BUFF, file_offset);

    assert(bin1_BUFF.size() == buff_size);
    assert(bin2_BUFF.size() == buff_size);
    assert(count_BUFF.size() == buff_size);

    for (usize j = 0; j < buff_size; ++j) {
      assert(count_BUFF[j] != 0);
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
        cmatrix.set(bin2, bin1, count_BUFF[j]);
      } else {
        // According to Cooler documentations, NaN means that a bin has been excluded by the
        // matrix balancing procedure. In this case we set the count to 0
        if (std::isnan(bin_weights[bin1]) || std::isnan(bin_weights[bin2])) {
          continue;  // Same as setting count to 0;
        }
        const auto bin1_bias = bin_weights[bin1];
        const auto bin2_bias = bin_weights[bin2];
        // See https://github.com/open2c/cooler/issues/35
        const auto count =
            static_cast<double>(count_BUFF[j]) / (bin1_bias * bin2_bias) / bias_scaling_factor;
        if constexpr (std::is_integral_v<N>) {
          cmatrix.set(bin2, bin1, utils::conditional_static_cast<N>(std::round(count)));
        } else {
          cmatrix.set(bin2, bin1, count);
        }
        DISABLE_WARNING_POP
      }
    }
  }

  return cmatrix;
}

template <class N>
usize Cooler<N>::stream_contacts_for_chrom(
    moodycamel::BlockingReaderWriterQueue<Cooler::PixelT> &queue, std::string_view chrom_name,
    usize nrows, std::pair<usize, usize> chrom_boundaries, bool try_common_chrom_prefixes,
    bool prefer_using_balanced_counts) {
  assert(this->_fp);
  assert(!this->_datasets.empty());
  assert(!this->_idx_bin1_offset.empty());
  assert(!this->_idx_chrom_offset.empty());

  const auto chrom_idx = this->get_chrom_idx(chrom_name, try_common_chrom_prefixes);
  const auto chrom_size = static_cast<usize>(this->get_chrom_sizes()[chrom_idx]);
  if (chrom_boundaries.second > chrom_size) {
    chrom_boundaries.second = chrom_size;
  }
  const auto bin_range = std::make_pair(static_cast<usize>(this->_idx_chrom_offset[chrom_idx]) +
                                            (chrom_boundaries.first / this->_bin_size),
                                        static_cast<usize>(this->_idx_chrom_offset[chrom_idx]) +
                                            (chrom_boundaries.second / this->_bin_size));
  const auto bin1_offset_idx = this->get_bin1_offset_idx_for_chrom(chrom_idx, chrom_boundaries);
  double pxl_count_scaling_factor = 1.0;

  if (prefer_using_balanced_counts &&
      hdf5::has_dataset(*this->_fp, "bins/weight", this->_root_path)) {
    const auto &d = this->_datasets[BIN_WEIGHT];
    u8 cis_only;
    try {
      hdf5::read_attribute(d, "cis_only", cis_only);
    } catch (const std::runtime_error &e) {
      if (absl::StrContains(e.what(), "cis_only")) {
        throw std::runtime_error(
            "File has a \"bins/weight\" dataset, but does not have an attribute named "
            "\"cis_only\". This most likely means that the file was generated by a very old "
            "version of Cooler. In order to proceed, you should rebalance the file using a "
            "recent "
            "version of Cooler, or in alternative remove the weight dataset (not recommended).");
      }
    }                // NOLINTNEXTLINE(readability-implicit-bool-conversion)
    if (cis_only) {  // --cis-only balancing produces an array of sale factors
      std::vector<double> buff;
      hdf5::read_attribute(d, "scale", buff);
      pxl_count_scaling_factor = buff[chrom_idx];
    } else {  // standard or --trans-only balancing produces a single scale factor
      hdf5::read_attribute(d, "scale", pxl_count_scaling_factor);
    }
  }

  return this->stream_contacts_for_chrom(queue, bin_range, bin1_offset_idx, nrows,
                                         pxl_count_scaling_factor, prefer_using_balanced_counts);
}

template <class N>
usize Cooler<N>::stream_contacts_for_chrom(
    moodycamel::BlockingReaderWriterQueue<Cooler::PixelT> &queue, std::string_view chrom_name,
    usize diagonal_width, usize bin_size, std::pair<usize, usize> chrom_boundaries,
    bool try_common_chrom_prefixes, bool prefer_using_balanced_counts) {
  const auto nrows =
      (diagonal_width / bin_size) + static_cast<usize>(diagonal_width % bin_size == 0);
  return this->stream_contacts_for_chrom(queue, chrom_name, nrows, chrom_boundaries,
                                         try_common_chrom_prefixes, prefer_using_balanced_counts);
}

template <class N>
// NOLINTNEXTLINE(readability-function-cognitive-complexity) TODO: reduce complexity
usize Cooler<N>::stream_contacts_for_chrom(moodycamel::BlockingReaderWriterQueue<PixelT> &queue,
                                           std::pair<hsize_t, hsize_t> bin_range,
                                           absl::Span<const i64> bin1_offset_idx, usize nrows,
                                           double bias_scaling_factor,
                                           bool prefer_using_balanced_counts) {
  if (this->_datasets.empty()) {
    this->open_default_datasets();
  }

  const auto &[first_bin, last_bin] = bin_range;
  assert(first_bin < last_bin);
  assert(last_bin <= first_bin + bin1_offset_idx.size());
  std::vector<i64> bin1_BUFF(nrows);
  std::vector<i64> bin2_BUFF(nrows);
  std::vector<i64> count_BUFF(nrows);
  std::vector<double> bin_weights;
  usize pixel_count = 0;

  const auto &d = this->_datasets;
  if (prefer_using_balanced_counts &&
      hdf5::has_dataset(*this->_fp, "bins/weight", this->_root_path)) {
    bin_weights.resize(last_bin - first_bin);
    std::ignore = hdf5::read_numbers(d[BIN_WEIGHT], bin_weights, static_cast<hsize_t>(first_bin));
  }

  for (usize i = 1; i < bin1_offset_idx.size(); ++i) {
    const auto file_offset = static_cast<hsize_t>(bin1_offset_idx[i - 1]);
    const auto buff_size =
        std::min(static_cast<usize>(bin1_offset_idx[i] - bin1_offset_idx[i - 1]), nrows);
    if (buff_size == 0) {
      continue;
    }

    bin1_BUFF.resize(buff_size);
    bin2_BUFF.resize(buff_size);
    count_BUFF.resize(buff_size);

    std::ignore = hdf5::read_numbers(d[PXL_B1], bin1_BUFF, file_offset);
    std::ignore = hdf5::read_numbers(d[PXL_B2], bin2_BUFF, file_offset);
    std::ignore = hdf5::read_numbers(d[PXL_COUNT], count_BUFF, file_offset);

    assert(bin1_BUFF.size() == buff_size);
    assert(bin2_BUFF.size() == buff_size);
    assert(count_BUFF.size() == buff_size);

    for (usize j = 0; j < buff_size; ++j) {
      assert(count_BUFF[j] != 0);
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
        while (!queue.try_emplace(PixelT{std::min(bin1, bin2), std::max(bin1, bin2),
                                         static_cast<contacts_t>(count_BUFF[j])})) {
          std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        ++pixel_count;
      } else {
        // According to Cooler documentations, NaN means that a bin has been excluded by the
        // matrix balancing procedure. In this case we set the count to 0
        if (std::isnan(bin_weights[bin1]) || std::isnan(bin_weights[bin2])) {
          continue;  // Same as setting count to 0;
        }
        const auto bin1_bias = bin_weights[bin1];
        const auto bin2_bias = bin_weights[bin2];
        // See https://github.com/open2c/cooler/issues/35
        const auto count =
            static_cast<double>(count_BUFF[j]) / (bin1_bias * bin2_bias) / bias_scaling_factor;
        while (!queue.try_emplace(PixelT{std::min(bin1, bin2), std::max(bin1, bin2),
                                         static_cast<contacts_t>(std::round(count))})) {
          std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        ++pixel_count;
        DISABLE_WARNING_POP
      }
    }
  }
  return pixel_count;
}

template <class N>
usize Cooler<N>::stream_contacts_for_chrom(moodycamel::BlockingReaderWriterQueue<PixelT> &queue,
                                           std::pair<hsize_t, hsize_t> bin_range,
                                           const std::vector<i64> &bin1_offset_idx, usize nrows,
                                           double scaling_factor,
                                           bool prefer_using_balanced_counts) {
  return this->stream_contacts_for_chrom(queue, bin_range, absl::MakeConstSpan(bin1_offset_idx),
                                         nrows, scaling_factor, prefer_using_balanced_counts);
}
}  // namespace modle::cooler
