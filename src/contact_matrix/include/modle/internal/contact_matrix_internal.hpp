// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <boost/serialization/access.hpp>
#include <ostream>
#include <type_traits>
#include <vector>

#include "modle/common/common.hpp"  // for utils::ndebug_defined
#include "modle/common/pixel.hpp"   // for PixelCoordinates

namespace modle::internal {
/// Map coordinates from the NxN to the NxM space of a ContactMatrixDense
//
// This function does not do any kind of validation on the input params.
// Passing coordinates that lie outside of the NxM space is UB
[[nodiscard]] constexpr PixelCoordinates transpose_coords(usize row, usize col) noexcept;
[[nodiscard]] constexpr PixelCoordinates transpose_coords(const PixelCoordinates& coords) noexcept;

/// Encode 2D bin coordinates as 1D coordinates
//
// This function does not do any kind of validation on the input params.
// Passing indices not mapping to the NxM space represented by a ContactMatrixDense is UB
[[nodiscard]] constexpr usize encode_idx(const PixelCoordinates& coords, usize nrows) noexcept;
[[nodiscard]] constexpr usize encode_idx(usize row, usize col, usize nrows) noexcept;

/// Decode 1D coordinate to 2D bin coordinates
//
// This function does not do any kind of validation on the input params.
// Passing and index mapping outside to the NxM space represented by a ContactMatrixDense is UB
[[nodiscard]] constexpr PixelCoordinates decode_idx(usize i, usize nrows) noexcept;

template <class ContactMatrix>
inline void bound_check_coords(const ContactMatrix& m, usize row, usize col);

template <class T>
[[nodiscard]] constexpr usize compute_num_cols_per_chunk(usize nrows,
                                                         usize max_chunk_size_bytes = 4096ULL
                                                                                      << 10ULL);
template <class SumT>
struct SerializedContactMatrixHeader {
  static_assert(std::is_arithmetic_v<SumT>);
  friend class boost::serialization::access;

  u64 nrows{0};
  u64 ncols{0};
  usize cols_per_chunk{0};
  SumT tot_contacts{0};
  usize nnz{0};
  usize updates_missed{0};

  std::vector<i64> chunk_offsets{};

  SerializedContactMatrixHeader() = default;
  template <class ContactMatrix>
  inline explicit SerializedContactMatrixHeader(const ContactMatrix& m);
  inline i64 serialize(std::ostream& out_stream) const;
  [[nodiscard]] inline usize num_chunks() const noexcept;

  static inline SerializedContactMatrixHeader<SumT> deserialize(std::istream& in_stream);

 private:
  template <class BoostArchive>
  inline void serialize(BoostArchive& ar, unsigned int version);
};

template <class N>
class SerializationTmpBuffers {
  static constexpr usize compute_size_from_bytes(usize bytes) noexcept;
  static constexpr usize default_size{compute_size_from_bytes(64ULL * (1024ULL << 10U))};

  std::vector<usize> _idx_buff;
  std::vector<N> _count_buff;
  std::vector<usize> _sorting_tmp_buff{};

 public:
  friend class boost::serialization::access;
  inline SerializationTmpBuffers();
  inline explicit SerializationTmpBuffers(usize size_bytes);

  [[nodiscard]] inline const std::vector<usize>& idx() const noexcept;
  [[nodiscard]] inline const std::vector<N>& counts() const noexcept;

  [[nodiscard]] inline usize size() const noexcept;
  [[nodiscard]] inline usize capacity() const noexcept;
  [[nodiscard]] inline bool empty() const noexcept;
  [[nodiscard]] inline bool full() const noexcept;
  inline void reserve(usize new_size);

  inline void sort();

  inline void clear() noexcept;
  inline void push_back(usize i, N count) noexcept;

  inline i64 serialize(std::ostream& out_stream, u32 compression_level = 1) const;
  static inline i64 deserialize(std::istream& in_stream, SerializationTmpBuffers& buff);

  inline void copy_to_buff(std::vector<N>& buff) const;
  template <class ContactMatrixSparse>
  inline void copy_to_buff(ContactMatrixSparse& matrix) const;

  template <class BoostArchive>
  inline void serialize(BoostArchive& ar, unsigned int version);
};

}  // namespace modle::internal

#include "../../../contact_matrix_internal_impl.hpp"  // IWYU pragma: export
// IWYU pragma: "../../../contact_matrix_internal_impl.hpp"
