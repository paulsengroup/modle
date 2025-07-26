// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "modle/common/common.hpp"
#include "modle/common/pixel.hpp"

namespace modle::internal {
/// Map coordinates from the NxN to the NxM space of a ContactMatrixDense
//
// This function does not do any kind of validation on the input params.
// Passing coordinates that lie outside of the NxM space is UB
[[nodiscard]] constexpr PixelCoordinates transpose_coords(std::size_t row,
                                                          std::size_t col) noexcept;
[[nodiscard]] constexpr PixelCoordinates transpose_coords(const PixelCoordinates& coords) noexcept;

/// Encode 2D bin coordinates as 1D coordinates
//
// This function does not do any kind of validation on the input params.
// Passing indices not mapping to the NxM space represented by a ContactMatrixDense is UB
[[nodiscard]] constexpr std::size_t encode_idx(const PixelCoordinates& coords,
                                               std::size_t nrows) noexcept;
[[nodiscard]] constexpr std::size_t encode_idx(std::size_t row, std::size_t col,
                                               std::size_t nrows) noexcept;

/// Decode 1D coordinate to 2D bin coordinates
//
// This function does not do any kind of validation on the input params.
// Passing and index mapping outside to the NxM space represented by a ContactMatrixDense is UB
[[nodiscard]] constexpr PixelCoordinates decode_idx(std::size_t i, std::size_t nrows) noexcept;

template <class ContactMatrix>
inline void bound_check_coords(const ContactMatrix& m, std::size_t row, std::size_t col);

template <class T>
[[nodiscard]] constexpr std::size_t compute_num_cols_per_chunk(
    std::size_t nrows, std::size_t max_chunk_size_bytes = 4096ULL << 10ULL);
}  // namespace modle::internal

#include "../../../contact_matrix_internal_impl.hpp"  // IWYU pragma: export
// IWYU pragma: "../../../contact_matrix_internal_impl.hpp"
