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
  std::size_t row;
  std::size_t col;
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
  constexpr Pixel(std::size_t row, std::size_t col, N count_) noexcept
      : Pixel(Coordinates{row, col}, count_) {}

  [[nodiscard]] constexpr std::size_t row() const noexcept { return coords.row; }
  [[nodiscard]] constexpr std::size_t col() const noexcept { return coords.col; }
  [[nodiscard]] constexpr bool operator==(const Pixel &other) const {
    return coords == other.coords && count == other.count;
  }
};

}  // namespace modle

namespace std {

template <>
struct hash<modle::PixelCoordinates> {
  std::size_t operator()(const typename modle::PixelCoordinates &coords) const {
    const std::array<std::size_t, 2> buff{coords.row, coords.col};
    return modle::utils::conditional_static_cast<std::size_t>(
        XXH3_64bits(buff.data(), sizeof(std::size_t) * buff.size()));
  }
};

}  // namespace std
