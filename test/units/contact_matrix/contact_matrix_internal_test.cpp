// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/internal/contact_matrix_internal.hpp"

#include <catch2/catch_test_macros.hpp>

#include "modle/common/common.hpp"  // for usize
#include "modle/common/random.hpp"

namespace modle::cmatrix::test {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrix internal: transpose_coords", "[cmatrix][short]") {
  const auto c1t = internal::transpose_coords(1, 2);
  CHECK(c1t == PixelCoordinates{1, 2});

  const auto c2t = internal::transpose_coords(2, 1);
  CHECK(c2t == PixelCoordinates{1, 2});

  const auto c3t = internal::transpose_coords(0, 0);
  CHECK(c3t == PixelCoordinates{0, 0});
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrix internal: encode_idx", "[cmatrix][short]") {
  constexpr usize nrows = 4;
  [[maybe_unused]] constexpr usize ncols = 4;

  const auto i1 = internal::encode_idx(PixelCoordinates{0, 0}, nrows);
  CHECK(i1 == 0);

  const auto i2 = internal::encode_idx(1, 2, nrows);
  CHECK(i2 == 9);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrix internal: decode_idx", "[cmatrix][short]") {
  constexpr usize nrows = 4;
  [[maybe_unused]] constexpr usize ncols = 4;

  const auto c1 = internal::decode_idx(0, nrows);
  CHECK(c1 == PixelCoordinates{0, 0});

  const auto c2 = internal::decode_idx(9, nrows);
  CHECK(c2 == PixelCoordinates{1, 2});
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrix internal: encode/decode roundtrip", "[cmatrix][short]") {
  constexpr usize nrows = 100;
  constexpr usize ncols = 5000;

  auto rand_eng = random::PRNG(1234567890ULL);
  for (usize i = 0; i < 1'000'000; ++i) {
    const auto col = random::uniform_int_distribution<usize>{0, ncols - 1}(rand_eng);
    const auto row =
        std::min(col + random::uniform_int_distribution<usize>{0, nrows - 1}(rand_eng), ncols - 1);

    const auto coordst = internal::transpose_coords(row, col);
    const auto idx = internal::encode_idx(coordst, nrows);
    CHECK(idx < nrows * ncols);
    const auto coordst_decoded = internal::decode_idx(idx, nrows);
    CHECK(coordst_decoded == coordst);
  }
}

}  // namespace modle::cmatrix::test
