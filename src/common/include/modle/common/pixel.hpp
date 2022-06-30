// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <xxhash.h>  // for XXH_INLINE_XXH3_64bits, XXH3_64bits

#include <array>        // for array
#include <functional>   // for hash
#include <type_traits>  // for is_arithmetic_v

#include "modle/common/common.hpp"  // for conditional_static_cast, usize

namespace modle {

struct PixelCoordinates {
  usize row;
  usize col;
  [[nodiscard]] constexpr bool operator==(const PixelCoordinates &other) const {
    return this->row == other.row && this->col == other.col;
  }
  [[nodiscard]] constexpr bool operator<(const PixelCoordinates &other) const {
    if (this->row == other.row) {
      return this->col < other.col;
    }
    return this->row < other.row;
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

  [[nodiscard]] constexpr usize row() const noexcept { return this->coords.row; }
  [[nodiscard]] constexpr usize col() const noexcept { return this->coords.col; }
  [[nodiscard]] constexpr bool operator==(const Pixel &other) const {
    return this->coords == other.coords && this->count == other.count;
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
