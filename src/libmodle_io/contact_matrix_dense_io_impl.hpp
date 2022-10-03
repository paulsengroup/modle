#pragma once

#include <cassert>
#include <coolerpp/coolerpp.hpp>
#include <filesystem>
#include <string_view>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/contact_matrix_dense.hpp"

namespace modle::io {

template <class N>
inline void append_contact_matrix_to_cooler(coolerpp::File& f, usize chrom_id,
                                            const ContactMatrixDense<N>& matrix) {
  constexpr usize chunk_size{64U << 10U};

  using PixelT = coolerpp::Pixel<N>;
  std::vector<PixelT> pixels;
  pixels.reserve(std::min(chunk_size, matrix.npixels()));

  const auto bin_table = f.bins();

  for (usize i = 0; i < matrix.ncols(); ++i) {
    for (usize j = i; j < matrix.ncols() && j - i < matrix.nrows(); ++j) {
      if (const auto n = matrix.unsafe_get(i, j); n != 0) {
        const auto bin1_start = utils::conditional_static_cast<u32>(i) * f.bin_size();
        const auto bin2_start = utils::conditional_static_cast<u32>(j) * f.bin_size();
        pixels.emplace_back(PixelT{
            {bin_table, utils::conditional_static_cast<u32>(chrom_id), bin1_start, bin2_start}, n});

        if (pixels.size() == pixels.capacity()) {
          f.append_pixels(pixels.begin(), pixels.end(), utils::ndebug_not_defined(), chunk_size);
          pixels.clear();
        }
      }
    }
  }

  if (!pixels.empty()) {
    f.append_pixels(pixels.begin(), pixels.end(), utils::ndebug_not_defined(), pixels.size());
  }
}

template <class N>
inline void append_contact_matrix_to_cooler(coolerpp::File& f, std::string_view chrom_name,
                                            const ContactMatrixDense<N>& matrix) {
  append_contact_matrix_to_cooler(f, static_cast<usize>(f.chromosomes().get_id(chrom_name)),
                                  matrix);
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

    auto selector = f.template fetch<N>(chrom.name, static_cast<u32>(start_pos), chrom.name,
                                        static_cast<u32>(end_pos));
    for (const auto pixel : selector) {
      if (pixel.coords.bin2_start - pixel.coords.bin1_start < diagonal_width) {
        const auto i = (pixel.coords.bin1_start - static_cast<u32>(start_pos)) / f.bin_size();
        const auto j = (pixel.coords.bin2_start - static_cast<u32>(start_pos)) / f.bin_size();
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
  const auto match = f.template fetch<double>(chrom_name, static_cast<u32>(start_pos), chrom_name,
                                              static_cast<u32>(end_pos));

  return match.begin() == match.end();
}

}  // namespace modle::io
