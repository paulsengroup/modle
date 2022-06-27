// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cpp-sort/sorters/pdq_sorter.h>
#include <fmt/format.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <cassert>
#include <numeric>
#include <ostream>
#include <stdexcept>

#include "modle/common/common.hpp"  // for ndebug_defined, ndebug_not_defined
#include "modle/common/pixel.hpp"   // for PixelCoordinates

namespace modle::internal {

constexpr PixelCoordinates transpose_coords(const usize row, const usize col) noexcept {
  if (row > col) {
    return {row - col, row};
  }
  return {col - row, col};
}

constexpr PixelCoordinates transpose_coords(const PixelCoordinates& coords) noexcept {
  return transpose_coords(coords.row, coords.col);
}

constexpr PixelCoordinates decode_idx(usize i, usize nrows) noexcept {
  assert(nrows != 0);
  const auto col = i / nrows;
  const auto row = i - (col * nrows);
  assert(row <= col);
  assert(row <= nrows);

  return {row, col};
}

constexpr usize encode_idx(usize row, usize col, usize nrows) noexcept {
  return (col * nrows) + row;
}

constexpr usize encode_idx(const PixelCoordinates& coords, usize nrows) noexcept {
  return encode_idx(coords.row, coords.col, nrows);
}

template <class ContactMatrix>
void bound_check_coords(const ContactMatrix& m, const usize row, const usize col) {
  if (MODLE_UNLIKELY(row >= m.ncols()) || col >= m.ncols()) {
    throw std::logic_error(
        fmt::format(FMT_STRING("Detected an out-of-bound read: attempt to access "
                               "item at {}:{} of a matrix of shape {}x{} ({}x{})"),
                    row, col, m.ncols(), m.ncols(), m.nrows(), m.ncols()));
  }
}

template <class T>
constexpr usize compute_num_cols_per_chunk(usize nrows, usize max_chunk_size_bytes) {
  const auto row_size_bytes = nrows * sizeof(T);
  return std::max(usize(1), max_chunk_size_bytes / row_size_bytes);
}

template <class SumT>
template <class ContactMatrix>
SerializedContactMatrixHeader<SumT>::SerializedContactMatrixHeader(const ContactMatrix& m)
    : nrows(m.nrows()),
      ncols(m.ncols()),
      tot_contacts(m.unsafe_get_tot_contacts()),
      nnz(m.unsafe_get_nnz()),
      updates_missed(m.get_n_of_missed_updates()) {
  using N = typename ContactMatrix::value_type;

  cols_per_chunk = compute_num_cols_per_chunk<N>(nrows);
  const auto num_chunks = (this->ncols + cols_per_chunk - 1) / cols_per_chunk;
  chunk_offsets.resize(num_chunks);
}

template <class SumT>
i64 SerializedContactMatrixHeader<SumT>::serialize(std::ostream& out_stream) const {
  namespace bio = boost::iostreams;
  boost::iostreams::filtering_ostreambuf fos;
  fos.push(out_stream);

  boost::archive::binary_oarchive bo(fos);
  bo << *this;

  return out_stream.tellp();
}

template <class SumT>
usize SerializedContactMatrixHeader<SumT>::num_chunks() const noexcept {
  return this->chunk_offsets.size();
}

template <class SumT>
SerializedContactMatrixHeader<SumT> SerializedContactMatrixHeader<SumT>::deserialize(
    std::istream& in_stream) {
  SerializedContactMatrixHeader<SumT> header{};

  boost::iostreams::filtering_istreambuf fis;
  fis.push(in_stream);

  boost::archive::binary_iarchive bi(fis);
  bi >> header;

  return header;
}

template <class SumT>
template <class BoostArchive>
void SerializedContactMatrixHeader<SumT>::serialize(BoostArchive& ar,
                                                    [[maybe_unused]] const unsigned int version) {
  ar& this->nrows;
  ar& this->ncols;
  ar& this->cols_per_chunk;
  ar& this->tot_contacts;
  ar& this->nnz;
  ar& this->updates_missed;

  ar& this->chunk_offsets;
}

template <class N>
SerializationTmpBuffers<N>::SerializationTmpBuffers()
    : _idx_buff(default_size), _count_buff(default_size) {
  this->clear();
}

template <class N>
SerializationTmpBuffers<N>::SerializationTmpBuffers(usize size)
    : _idx_buff(compute_size_from_bytes(size)), _count_buff(compute_size_from_bytes(size)) {
  this->clear();
}

template <class N>
const std::vector<usize>& SerializationTmpBuffers<N>::idx() const noexcept {
  return this->_idx_buff;
}

template <class N>
const std::vector<N>& SerializationTmpBuffers<N>::counts() const noexcept {
  return this->_count_buff;
}

template <class N>
usize SerializationTmpBuffers<N>::size() const noexcept {
  assert(this->_idx_buff.size() == this->_count_buff.size());
  return this->_idx_buff.size();
}

template <class N>
usize SerializationTmpBuffers<N>::capacity() const noexcept {
  assert(this->_idx_buff.capacity() == this->_count_buff.capacity());
  return this->_idx_buff.capacity();
}

template <class N>
bool SerializationTmpBuffers<N>::empty() const noexcept {
  return this->size() == 0;
}

template <class N>
bool SerializationTmpBuffers<N>::full() const noexcept {
  return this->size() == this->capacity();
}

template <class N>
void SerializationTmpBuffers<N>::reserve(usize new_size) {
  this->_idx_buff.reserve(new_size);
  this->_count_buff.reserve(new_size);
}

template <class N>
void SerializationTmpBuffers<N>::sort() {
  this->_sorting_tmp_buff.resize(this->size());
  std::iota(this->_sorting_tmp_buff.begin(), this->_sorting_tmp_buff.end(), usize(0));

  cppsort::pdq_sort(
      this->_sorting_tmp_buff.begin(), this->_sorting_tmp_buff.end(),
      [&](const auto i1, const auto i2) { return this->_idx_buff[i1] < this->_idx_buff[i2]; });

  for (usize i = 0; i < this->size(); ++i) {
    // https://stackoverflow.com/a/22218699
    while (this->_sorting_tmp_buff[i] != i) {
      const auto j = this->_sorting_tmp_buff[i];
      const auto k = this->_sorting_tmp_buff[j];

      std::swap(this->_idx_buff[j], this->_idx_buff[k]);
      std::swap(this->_count_buff[j], this->_count_buff[k]);
      std::swap(this->_sorting_tmp_buff[i], this->_sorting_tmp_buff[j]);
    }
  }
}

template <class N>
void SerializationTmpBuffers<N>::clear() noexcept {
  this->_idx_buff.clear();
  this->_count_buff.clear();
}

template <class N>
void SerializationTmpBuffers<N>::push_back(usize i, N count) noexcept {
  assert(this->size() <= this->capacity());

  this->_idx_buff.push_back(i);
  this->_count_buff.push_back(count);

  assert(this->size() <= this->capacity());
}

template <class N>
i64 SerializationTmpBuffers<N>::serialize(std::ostream& out_stream, u32 compression_level) const {
  namespace bio = boost::iostreams;
  bio::filtering_ostreambuf fos;

  // push the ostream and the compressor
  fos.push(bio::zstd_compressor(bio::zstd_params(compression_level)));
  fos.push(out_stream);

  // start the archive on the filtering buffer:
  boost::archive::binary_oarchive bo(fos);
  bo << *this;

  return out_stream.tellp();
}

template <class N>
i64 SerializationTmpBuffers<N>::deserialize(std::istream& in_stream,
                                            SerializationTmpBuffers<N>& buff) {
  namespace bio = boost::iostreams;

  bio::filtering_istreambuf fis;
  fis.push(bio::zstd_decompressor());
  fis.push(in_stream);
  boost::archive::binary_iarchive bi(fis);
  buff.clear();

  bi >> buff;

  return in_stream.tellg();
}

template <class N>
void SerializationTmpBuffers<N>::copy_to_buff(std::vector<N>& buff) const {
  for (usize i = 0; i < this->size(); ++i) {
    const auto j = this->_idx_buff[i];
    buff[j] = this->_count_buff[i];
  }
}

template <class N>
template <class ContactMatrixSparse>
void SerializationTmpBuffers<N>::copy_to_buff(ContactMatrixSparse& matrix) const {
  for (usize i = 0; i < this->size(); ++i) {
    const auto [rowt, colt] = internal::decode_idx(this->_idx_buff[i], matrix.nrows());
    assert(colt >= rowt);
    const auto row = colt - rowt;
    const auto col = colt;
    matrix.set(row, col, this->_count_buff[i]);
    assert(matrix.get_n_of_missed_updates() == 0);
  }
}

template <class N>
template <class BoostArchive>
void SerializationTmpBuffers<N>::serialize(BoostArchive& ar,
                                           [[maybe_unused]] const unsigned int version) {
  ar& this->_idx_buff;
  ar& this->_count_buff;
}

template <class N>
constexpr usize SerializationTmpBuffers<N>::compute_size_from_bytes(usize bytes) noexcept {
  constexpr usize idx_size = sizeof(typename decltype(_idx_buff)::value_type);
  constexpr usize value_size = sizeof(N);

  return std::max(usize(1), bytes / (idx_size + value_size));
}
}  // namespace modle::internal
