// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <xxhash.h>

#include <array>
#include <functional>
#include <type_traits>

#include "modle/common/common.hpp"

namespace modle {

struct PixelCoordinates {
  usize row;
  usize col;
  [[nodiscard]] constexpr bool operator==(const PixelCoordinates &other) const {
    return row == other.row && col == other.col;
  }
  [[nodiscard]] constexpr bool operator<(const PixelCoordinates &other) const {
    if (row == other.row) {
      return col < other.col;
    }
    return row < other.row;
  }
};

template <class N>
struct Pixel {
  static_assert(std::is_arithmetic_v<N>);
  using Coordinates = PixelCoordinates;

  Coordinates coords{0, 0};

  N count{0};

  constexpr Pixel() = default;
  constexpr Pixel(Coordinates coords_, N count_) noexcept : coords(coords_), count(count_) {}
  constexpr Pixel(usize row, usize col, N count_) noexcept : Pixel(Coordinates{row, col}, count_) {}

  [[nodiscard]] constexpr usize row() const noexcept { return coords.row; }
  [[nodiscard]] constexpr usize col() const noexcept { return coords.col; }
  [[nodiscard]] constexpr bool operator==(const Pixel &other) const {
    return coords == other.coords && count == other.count;
  }
};

}  // namespace modle

namespace std {

template <>
struct hash<modle::PixelCoordinates> {
  using usize = modle::usize;
  usize operator()(const typename modle::PixelCoordinates &coords) const {
    const std::array<usize, 2> buff{coords.row, coords.col};
    return modle::utils::conditional_static_cast<usize>(
        XXH3_64bits(buff.data(), sizeof(usize) * buff.size()));
  }
};

}  // namespace std
