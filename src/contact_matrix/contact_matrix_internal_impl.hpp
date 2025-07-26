// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>

#include "modle/common/common.hpp"
#include "modle/common/pixel.hpp"

namespace modle::internal {

constexpr PixelCoordinates transpose_coords(const std::size_t row, const std::size_t col) noexcept {
  if (row > col) {
    return {row - col, row};
  }
  return {col - row, col};
}

constexpr PixelCoordinates transpose_coords(const PixelCoordinates& coords) noexcept {
  return transpose_coords(coords.row, coords.col);
}

constexpr PixelCoordinates decode_idx(std::size_t i, std::size_t nrows) noexcept {
  assert(nrows != 0);
  const auto col = i / nrows;
  const auto row = i - (col * nrows);
  assert(row <= col);
  assert(row <= nrows);

  return {row, col};
}

constexpr std::size_t encode_idx(std::size_t row, std::size_t col, std::size_t nrows) noexcept {
  return (col * nrows) + row;
}

constexpr std::size_t encode_idx(const PixelCoordinates& coords, std::size_t nrows) noexcept {
  return encode_idx(coords.row, coords.col, nrows);
}

template <class ContactMatrix>
void bound_check_coords(const ContactMatrix& m, const std::size_t row, const std::size_t col) {
  if (row >= m.ncols() || col >= m.ncols()) [[unlikely]] {
    throw std::logic_error(
        fmt::format("detected an out-of-bound read: attempt to access "
                    "item at {}:{} of a matrix of shape {}x{} ({}x{})",
                    row, col, m.ncols(), m.ncols(), m.nrows(), m.ncols()));
  }
}

template <class T>
constexpr std::size_t compute_num_cols_per_chunk(std::size_t nrows,
                                                 std::size_t max_chunk_size_bytes) {
  const auto row_size_bytes = nrows * sizeof(T);
  return std::max(std::size_t(1), max_chunk_size_bytes / row_size_bytes);
}

}  // namespace modle::internal
