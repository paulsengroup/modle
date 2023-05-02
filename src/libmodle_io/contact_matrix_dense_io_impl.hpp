// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <coolerpp/coolerpp.hpp>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/contact_matrix_dense.hpp"

namespace modle::io {

namespace internal {

template <class PixelT, class N>
inline bool emplace_pixel(const coolerpp::File& f, std::vector<PixelT>& buff, usize chrom_id,
                          bp_t offset, usize bin1_id, usize bin2_id, N n) {
  using T = decltype(std::declval<PixelT>().count);

  const auto bin1_start = utils::conditional_static_cast<u32>(offset) +
                          (utils::conditional_static_cast<u32>(bin1_id) * f.bin_size());
  const auto bin2_start = utils::conditional_static_cast<u32>(offset) +
                          (utils::conditional_static_cast<u32>(bin2_id) * f.bin_size());
  buff.emplace_back(PixelT{{f.bins_ptr(), static_cast<u32>(chrom_id), bin1_start, bin2_start},
                           utils::conditional_static_cast<T>(n)});

  return buff.size() == buff.capacity();
}

template <class PixelT>
inline void write_pixels(coolerpp::File& f, std::vector<PixelT>& buff) {
  f.append_pixels(buff.begin(), buff.end(), utils::ndebug_not_defined(), buff.size());
  buff.clear();
}

template <class N, class PixelT, usize chunk_size>
inline void append_contact_matrix_to_cooler(coolerpp::File& f, usize chrom_id,
                                            const ContactMatrixDense<N>& matrix,
                                            std::vector<PixelT>& buff, bp_t offset) {
  assert(f.has_pixel_of_type<i32>() || f.has_pixel_of_type<double>());
  buff.reserve(std::min(chunk_size, matrix.npixels()));

  for (usize i = 0; i < matrix.ncols(); ++i) {
    for (usize j = i; j < matrix.ncols() && j - i < matrix.nrows(); ++j) {
      if (const auto n = matrix.unsafe_get(i, j); n != 0) {
        const auto buffer_is_full = emplace_pixel(f, buff, chrom_id, offset, i, j, n);
        if (buffer_is_full) {
          write_pixels(f, buff);
        }
      }
    }
  }

  if (!buff.empty()) {
    write_pixels(f, buff);
  }
}
}  // namespace internal

template <class N, class ChromIt>
inline coolerpp::File init_cooler_file(const std::filesystem::path& path, bool force_overwrite,
                                       ChromIt first_chrom, ChromIt last_chrom, usize bin_size,
                                       std::string_view assembly, std::string_view generated_by,
                                       std::string_view metadata_str) {
  const auto num_chroms = static_cast<usize>(std::distance(first_chrom, last_chrom));
  std::vector<std::string> chrom_names(num_chroms);
  std::vector<u32> chrom_sizes(num_chroms);

  std::transform(first_chrom, last_chrom, chrom_names.begin(),
                 [](const auto& chrom) { return std::string{chrom.name()}; });
  std::transform(first_chrom, last_chrom, chrom_sizes.begin(), [](const auto& chrom) {
    return utils::conditional_static_cast<u32>(chrom.size());
  });

  return init_cooler_file<N>(path, force_overwrite, chrom_names.begin(), chrom_names.end(),
                             chrom_sizes.begin(), bin_size, assembly, generated_by, metadata_str);
}

template <class N, class ChromIt>
[[nodiscard]] inline coolerpp::File init_cooler_file(const std::filesystem::path& path,
                                                     bool force_overwrite, ChromIt first_chrom,
                                                     ChromIt last_chrom,
                                                     coolerpp::StandardAttributes attrs) {
  const auto num_chroms = static_cast<usize>(std::distance(first_chrom, last_chrom));
  std::vector<std::string> chrom_names(num_chroms);
  std::vector<u32> chrom_sizes(num_chroms);

  std::transform(first_chrom, last_chrom, chrom_names.begin(),
                 [](const auto& chrom) { return std::string{chrom.name()}; });
  std::transform(first_chrom, last_chrom, chrom_sizes.begin(), [](const auto& chrom) {
    return utils::conditional_static_cast<u32>(chrom.size());
  });

  return init_cooler_file<N>(path, force_overwrite, chrom_names.begin(), chrom_names.end(),
                             chrom_sizes.begin(), std::move(attrs));
}

template <class N, class ChromNameIt, class ChromSizeIt>
inline coolerpp::File init_cooler_file(const std::filesystem::path& path, bool force_overwrite,
                                       ChromNameIt first_name, ChromNameIt last_name,
                                       ChromSizeIt first_size, usize bin_size,
                                       std::string_view assembly, std::string_view generated_by,
                                       std::string_view metadata_str) {
  assert(!assembly.empty());
  assert(!generated_by.empty());
  auto attrs = coolerpp::StandardAttributes::init(utils::conditional_static_cast<u32>(bin_size));
  attrs.assembly = assembly;
  attrs.generated_by = generated_by;
  if (!metadata_str.empty()) {
    attrs.metadata = std::string{metadata_str};
  }

  return init_cooler_file<N>(path, force_overwrite, first_name, last_name, first_size,
                             std::move(attrs));
}

template <class N, class ChromNameIt, class ChromSizeIt>
inline coolerpp::File init_cooler_file(const std::filesystem::path& path, bool force_overwrite,
                                       ChromNameIt first_name, ChromNameIt last_name,
                                       ChromSizeIt first_size, coolerpp::StandardAttributes attrs) {
  return coolerpp::File::create_new_cooler<N>(
      path.string(), coolerpp::ChromosomeSet{first_name, last_name, first_size}, attrs.bin_size,
      force_overwrite, std::move(attrs));
}

template <class N>
inline coolerpp::File init_cooler_file(const std::filesystem::path& path, bool force_overwrite,
                                       coolerpp::ChromosomeSet chroms,
                                       coolerpp::StandardAttributes attrs) {
  return coolerpp::File::create_new_cooler<N>(path.string(), std::move(chroms), attrs.bin_size,
                                              force_overwrite, std::move(attrs));
}

template <class N, usize chunk_size>
inline void append_contact_matrix_to_cooler(coolerpp::File& f, std::string_view chrom_name,
                                            const ContactMatrixDense<N>& matrix, bp_t offset) {
  append_contact_matrix_to_cooler<N, chunk_size>(
      f, static_cast<usize>(f.chromosomes().get_id(chrom_name)), matrix, offset);
}

template <class N, usize chunk_size>
inline void append_contact_matrix_to_cooler(coolerpp::File& f, usize chrom_id,
                                            const ContactMatrixDense<N>& matrix, bp_t offset) {
  using T = std::conditional_t<std::is_floating_point_v<N>, double, i32>;
  using PixelT = coolerpp::Pixel<T>;
  std::vector<PixelT> buff{};

  internal::append_contact_matrix_to_cooler<N, PixelT, chunk_size>(f, chrom_id, matrix, buff,
                                                                   offset);
}

template <class N>
inline ContactMatrixDense<N> read_contact_matrix_from_cooler(
    const std::filesystem::path& path_to_cooler, std::string_view chrom_name, bp_t start_pos,
    bp_t end_pos, bp_t diagonal_width) {
  return read_contact_matrix_from_cooler<N>(coolerpp::File::open_read_only(path_to_cooler.string()),
                                            chrom_name, start_pos, end_pos, diagonal_width);
}

template <class N>
inline ContactMatrixDense<N> read_contact_matrix_from_cooler(
    const std::filesystem::path& path_to_cooler, usize chrom_id, bp_t start_pos, bp_t end_pos,
    bp_t diagonal_width) {
  return read_contact_matrix_from_cooler<N>(coolerpp::File::open_read_only(path_to_cooler.string()),
                                            chrom_id, start_pos, end_pos, diagonal_width);
}

template <class N>
inline ContactMatrixDense<N> read_contact_matrix_from_cooler(const coolerpp::File& f,
                                                             std::string_view chrom_name,
                                                             bp_t start_pos, bp_t end_pos,
                                                             bp_t diagonal_width) {
  return read_contact_matrix_from_cooler<N>(f,
                                            static_cast<usize>(f.chromosomes().get_id(chrom_name)),
                                            start_pos, end_pos, diagonal_width);
}

template <class N>
inline ContactMatrixDense<N> read_contact_matrix_from_cooler(const coolerpp::File& f,
                                                             usize chrom_id, bp_t start_pos,
                                                             bp_t end_pos, bp_t diagonal_width) {
  try {
    assert(start_pos < end_pos);
    assert(f.bin_size() != 0);
    const auto& chrom = f.chromosomes().at(static_cast<u32>(chrom_id));

    if (end_pos > chrom.size) {
      if (end_pos != std::numeric_limits<bp_t>::max()) {
        throw std::runtime_error(fmt::format(FMT_STRING("Interval {}:{}-{} lies outside of {}"),
                                             chrom.name, start_pos, end_pos, chrom));
      }
      end_pos = chrom.size;
    }

    const auto interval_span = end_pos - start_pos;
    if (diagonal_width == 0) {
      diagonal_width = interval_span;
    }

    ContactMatrixDense<N> m{interval_span, diagonal_width, f.bin_size()};
    for (auto&& pixel :
         f.fetch<N>(chrom.name, static_cast<u32>(start_pos), static_cast<u32>(end_pos))) {
      if (pixel.coords.bin2().start - pixel.coords.bin1().start < diagonal_width) {
        const auto i = (pixel.coords.bin1().start - static_cast<u32>(start_pos)) / f.bin_size();
        const auto j = (pixel.coords.bin2().start - static_cast<u32>(start_pos)) / f.bin_size();
        m.unsafe_set(i, j, pixel.count);
      }
    }
    return m;
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to read contact matrix for {}:{}-{} from Cooler at URI \"{}\": {}"),
        f.chromosomes().at(static_cast<u32>(chrom_id)), start_pos, end_pos, f.uri(), e.what()));
  }
}

inline bool query_returns_no_pixels(const coolerpp::File& f, std::string_view chrom_name,
                                    bp_t start_pos, bp_t end_pos) {
  const auto match =
      f.fetch<double>(chrom_name, static_cast<u32>(start_pos), static_cast<u32>(end_pos));

  return match.begin() == match.end();
}

}  // namespace modle::io
