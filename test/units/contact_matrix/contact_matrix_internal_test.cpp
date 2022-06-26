// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/internal/contact_matrix_internal.hpp"

#include <catch2/catch_test_macros.hpp>

#include "modle/common/common.hpp"  // for usize

namespace modle::test::cmatrix {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix internal: transpose_coords", "[cmatrix][short]") {
  const auto c1t = internal::transpose_coords(1, 2);
  CHECK(c1t == PixelCoordinates{1, 2});

  const auto c2t = internal::transpose_coords(2, 1);
  CHECK(c2t == PixelCoordinates{1, 2});

  const auto c3t = internal::transpose_coords(0, 0);
  CHECK(c3t == PixelCoordinates{0, 0});
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix internal: encode_idx", "[cmatrix][short]") {
  constexpr usize nrows = 4;
  [[maybe_unused]] constexpr usize ncols = 4;

  const auto i1 = internal::encode_idx(PixelCoordinates{0, 0}, nrows);
  CHECK(i1 == 0);

  const auto i2 = internal::encode_idx(1, 2, nrows);
  CHECK(i2 == 9);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix internal: decode_idx", "[cmatrix][short]") {
  constexpr usize nrows = 4;
  [[maybe_unused]] constexpr usize ncols = 4;

  const auto c1 = internal::decode_idx(0, nrows);
  CHECK(c1 == PixelCoordinates{0, 0});

  const auto c2 = internal::decode_idx(9, nrows);
  CHECK(c2 == PixelCoordinates{1, 2});
}

}  // namespace modle::test::cmatrix
