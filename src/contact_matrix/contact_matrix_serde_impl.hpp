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
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#include <cassert>
#include <numeric>
#include <ostream>
#include <stdexcept>

#include "modle/common/common.hpp"  // for ndebug_defined, ndebug_not_defined
#include "modle/common/pixel.hpp"   // for PixelCoordinates

namespace modle {

template <class N>
constexpr bool ContactMatrixSerde<N>::ChunkMetadata::operator<(
    const ChunkMetadata& other) noexcept {
  assert(this->first_col < this->last_col);
  assert(other.first_col < other.last_col);

  return this->first_col < other.first_col;
}

template <class N>
constexpr bool ContactMatrixSerde<N>::ChunkMetadata::operator==(
    const ChunkMetadata& other) noexcept {
  assert(this->first_col < this->last_col);
  assert(other.first_col < other.last_col);

  return this->first_col == other.first_col && this->last_col == other.last_col;
}

template <class N>
template <class BoostArchive>
void ContactMatrixSerde<N>::ChunkMetadata::serialize(BoostArchive& ar,
                                                     [[maybe_unused]] unsigned int version) {
  ar& this->offset;
  ar& this->nnz;
  ar& this->sum;
}

template <class N>
constexpr usize ContactMatrixSerde<N>::Header::compute_num_cols_per_chunk(
    usize nrows, usize max_chunk_size_bytes) noexcept {
  const auto row_size_bytes = nrows * sizeof(N);
  return std::max(usize(1), max_chunk_size_bytes / row_size_bytes);
}

template <class N>
template <class ContactMatrix>
ContactMatrixSerde<N>::Header::Header(const ContactMatrix& m, std::string_view metadata_str)
    : metadata(metadata_str),
      nrows(m.nrows()),
      ncols(m.ncols()),
      tot_contacts(m.unsafe_get_tot_contacts()),
      nnz(m.unsafe_get_nnz()),
      updates_missed(m.get_n_of_missed_updates()) {
  if constexpr (std::is_same_v<ContactMatrixDense<N>, std::remove_cv_t<ContactMatrix>>) {
    _num_chunks = (nnz + ChunkBuffer::DEFAULT_CAPACITY - 1) / ChunkBuffer::DEFAULT_CAPACITY;
  } else {
    _num_chunks = m.num_blocks();
  }
  chunk_metadata.reserve(_num_chunks);
}

template <class N>
usize ContactMatrixSerde<N>::Header::num_chunks() const noexcept {
  return this->_num_chunks;
}

template <class N>
void ContactMatrixSerde<N>::Header::clear() noexcept {
  this->metadata.clear();
  this->nrows = 0;
  this->ncols = 0;
  this->tot_contacts = 0;
  this->nnz = 0;
  this->updates_missed = 0;
  this->_num_chunks = 0;
  this->chunk_metadata.clear();
}

template <class N>
template <class BoostArchive>
void ContactMatrixSerde<N>::Header::save(BoostArchive& ar,
                                         [[maybe_unused]] const unsigned int version) const {
  ar << this->nrows;
  ar << this->ncols;
  ar << this->tot_contacts;
  ar << this->nnz;
  ar << this->updates_missed;

  assert(this->chunk_metadata.is_BST());
  assert(this->_num_chunks == this->chunk_metadata.size());

  std::vector<usize> idx_buff(chunk_metadata.capacity());
  std::copy(this->chunk_metadata.starts_begin(), this->chunk_metadata.starts_end(),
            idx_buff.begin());
  ar << idx_buff;

  std::copy(this->chunk_metadata.ends_begin(), this->chunk_metadata.ends_end(), idx_buff.begin());
  ar << idx_buff;

  std::vector<ChunkMetadata> chunk_metadata_flat(chunk_metadata.capacity());
  std::copy(this->chunk_metadata.data_begin(), this->chunk_metadata.data_end(),
            chunk_metadata_flat.begin());
  ar << chunk_metadata_flat;
}

template <class N>
template <class BoostArchive>
void ContactMatrixSerde<N>::Header::load(BoostArchive& ar,
                                         [[maybe_unused]] const unsigned int version) {
  ar >> this->nrows;
  ar >> this->ncols;
  ar >> this->tot_contacts;
  ar >> this->nnz;
  ar >> this->updates_missed;

  std::vector<usize> idx_buff1(chunk_metadata.capacity());
  std::vector<usize> idx_buff2(chunk_metadata.capacity());
  std::vector<ChunkMetadata> chunk_metadata_flat(chunk_metadata.capacity());

  ar >> idx_buff1;
  ar >> idx_buff2;
  ar >> chunk_metadata_flat;

  assert(idx_buff1.size() == chunk_metadata.capacity());
  assert(idx_buff2.size() == chunk_metadata.capacity());
  assert(chunk_metadata_flat.size() == chunk_metadata.capacity());

  for (usize i = 0; i < chunk_metadata.capacity(); ++i) {
    this->chunk_metadata.insert(idx_buff1[i], idx_buff2[i], chunk_metadata_flat[i]);
  }

  this->chunk_metadata.make_BST();
  this->_num_chunks = this->chunk_metadata.size();
}

template <class N>
template <class BoostArchive>
void ContactMatrixSerde<N>::Header::serialize(BoostArchive& ar, const unsigned int version) {
  boost::serialization::split_member(ar, *this, version);
}

template <class N>
ContactMatrixSerde<N>::ChunkBuffer::ChunkBuffer()
    : _idx_buff(DEFAULT_CAPACITY), _count_buff(DEFAULT_CAPACITY) {
  this->clear();
}

template <class N>
usize ContactMatrixSerde<N>::ChunkBuffer::size() const noexcept {
  assert(this->_idx_buff.size() == this->_count_buff.size());
  return this->_idx_buff.size();
}

template <class N>
usize ContactMatrixSerde<N>::ChunkBuffer::capacity() const noexcept {
  assert(this->_idx_buff.capacity() == this->_count_buff.capacity());
  return this->_idx_buff.capacity();
}

template <class N>
bool ContactMatrixSerde<N>::ChunkBuffer::empty() const noexcept {
  return this->size() == 0;
}

template <class N>
bool ContactMatrixSerde<N>::ChunkBuffer::full() const noexcept {
  return this->size() == this->capacity();
}

template <class N>
void ContactMatrixSerde<N>::ChunkBuffer::reserve(usize new_size) {
  this->_idx_buff.reserve(new_size);
  this->_count_buff.reserve(new_size);
}

template <class N>
void ContactMatrixSerde<N>::ChunkBuffer::sort(std::vector<usize>& sorting_buff) {
  sorting_buff.resize(this->size());
  std::iota(sorting_buff.begin(), sorting_buff.end(), usize(0));

  cppsort::pdq_sort(sorting_buff.begin(), sorting_buff.end(), [&](const auto i1, const auto i2) {
    return this->_idx_buff[i1] < this->_idx_buff[i2];
  });

  for (usize i = 0; i < this->size(); ++i) {
    // https://stackoverflow.com/a/22218699
    while (sorting_buff[i] != i) {
      const auto j = sorting_buff[i];
      const auto k = sorting_buff[j];

      std::swap(this->_idx_buff[j], this->_idx_buff[k]);
      std::swap(this->_count_buff[j], this->_count_buff[k]);
      std::swap(sorting_buff[i], sorting_buff[j]);
    }
  }
}

template <class N>
auto ContactMatrixSerde<N>::ChunkBuffer::sum_contacts() const noexcept -> SumT {
  return std::accumulate(this->_count_buff.begin(), this->_count_buff.end(), SumT(0));
}

template <class N>
void ContactMatrixSerde<N>::ChunkBuffer::clear() noexcept {
  this->_idx_buff.clear();
  this->_count_buff.clear();
}

template <class N>
void ContactMatrixSerde<N>::ChunkBuffer::push_back(usize i, N count) noexcept {
  this->_idx_buff.push_back(i);
  this->_count_buff.push_back(count);
}

template <class N>
i64 ContactMatrixSerde<N>::ChunkBuffer::serialize(std::ostream& stream,
                                                  u32 compression_level) const {
  namespace bio = boost::iostreams;
  bio::filtering_ostreambuf fos;

  // push the ostream and the compressor
  fos.push(bio::zstd_compressor(bio::zstd_params(compression_level)));
  fos.push(stream);

  // start the archive on the filtering buffer:
  boost::archive::binary_oarchive bo(fos);
  bo << *this;

  return stream.tellp();
}

template <class N>
i64 ContactMatrixSerde<N>::ChunkBuffer::deserialize(std::istream& stream,
                                                    ContactMatrixSerde<N>::ChunkBuffer& buff) {
  namespace bio = boost::iostreams;

  bio::filtering_istreambuf fis;
  fis.push(bio::zstd_decompressor());
  fis.push(stream);
  boost::archive::binary_iarchive bi(fis);
  buff.clear();

  bi >> buff;

  return stream.tellg();
}

template <class N>
template <class BoostArchive>
void ContactMatrixSerde<N>::ChunkBuffer::serialize(BoostArchive& ar,
                                                   [[maybe_unused]] const unsigned int version) {
  ar& this->_idx_buff;
  ar& this->_count_buff;
}

template <class N>
constexpr usize ContactMatrixSerde<N>::ChunkBuffer::compute_size_from_bytes(usize bytes) noexcept {
  constexpr usize idx_size = sizeof(typename decltype(_idx_buff)::value_type);
  constexpr usize value_size = sizeof(N);

  return std::max(usize(1), bytes / (idx_size + value_size));
}

template <class N>
void ContactMatrixSerde<N>::write_header(std::ostream& stream) {
  stream.seekp(std::ios::beg);
  this->_header.chunk_metadata.make_BST();

  boost::iostreams::filtering_ostreambuf fos;
  fos.push(stream);

  boost::archive::binary_oarchive bo(fos);
  bo << this->_header;
}

template <class N>
void ContactMatrixSerde<N>::write_empty_header(std::ostream& stream) {
  stream.seekp(std::ios::beg);

  const auto num_chunks = this->_header.chunk_metadata.capacity();

  for (usize i = 0; i < num_chunks; ++i) {
    this->_header.chunk_metadata.insert(i, i + 1, ChunkMetadata{});
  }
  this->_header.chunk_metadata.make_BST();

  boost::iostreams::filtering_ostreambuf fos;
  fos.push(stream);

  boost::archive::binary_oarchive bo(fos);
  bo << this->_header;
  this->_header.chunk_metadata.clear();
}

template <class N>
void ContactMatrixSerde<N>::read_header(std::istream& stream) {
  stream.seekg(std::ios::beg);

  boost::iostreams::filtering_istreambuf fis;
  fis.push(stream);

  boost::archive::binary_iarchive bi(fis);
  bi >> this->_header;
}

template <class N>
template <class ContactMatrix>
void ContactMatrixSerde<N>::read(const std::filesystem::path& in_path, ContactMatrix& matrix) {
  std::ifstream stream(in_path, std::ios::binary);
  stream.exceptions(std::ios::failbit);
  this->read_internal(stream, matrix);
}

template <class N>
template <class ContactMatrix>
void ContactMatrixSerde<N>::read(std::istream& stream, ContactMatrix& matrix) {
  const auto except_policy = stream.exceptions();
  stream.exceptions(std::ios::failbit);
  this->read_internal(stream, matrix);
  stream.exceptions(except_policy);
}

template <class N>
template <class ContactMatrix>
void ContactMatrixSerde<N>::write(const std::filesystem::path& out_path,
                                  const ContactMatrix& matrix, std::string_view metadata) {
  std::ofstream stream(out_path, std::ios::binary);
  stream.exceptions(std::ios::failbit);
  this->write_internal(stream, matrix, metadata);
}

template <class N>
template <class ContactMatrix>
void ContactMatrixSerde<N>::write(std::ostream& stream, const ContactMatrix& matrix,
                                  std::string_view metadata) {
  const auto except_policy = stream.exceptions();
  stream.exceptions(std::ios::failbit);
  this->write_internal(stream, matrix, metadata);
  stream.exceptions(except_policy);
}

template <class N>
template <class ContactMatrix>
void ContactMatrixSerde<N>::write_destructive(const std::filesystem::path& out_path,
                                              ContactMatrix& matrix, std::string_view metadata) {
  std::ofstream stream(out_path, std::ios::binary);
  stream.exceptions(std::ios::failbit);
  this->write_internal_destructive(stream, matrix, metadata);
}

template <class N>
template <class ContactMatrix>
void ContactMatrixSerde<N>::write_destructive(std::ostream& stream, ContactMatrix& matrix,
                                              std::string_view metadata) {
  const auto except_policy = stream.exceptions();
  this->write_internal_destructive(stream, matrix, metadata);
  stream.exceptions(except_policy);
}

template <class N>
void ContactMatrixSerde<N>::read_internal(std::istream& stream, ContactMatrixDense<N>& matrix) {
  stream.seekg(std::ios::beg);
  this->read_header(stream);

  matrix.unsafe_reset();
  matrix.unsafe_resize(this->_header.nrows, this->_header.ncols);

  for (const auto& chunk_metadata : this->_header.chunk_metadata.data()) {
    assert(chunk_metadata.offset != -1);
    this->read_chunk(stream, chunk_metadata.offset);
    assert(chunk_metadata.nnz == this->_chunk.size());
    assert(chunk_metadata.sum == this->_chunk.sum_contacts());

    for (usize i = 0; i < this->_chunk.size(); ++i) {
      const auto j = this->_chunk._idx_buff[i];
      assert(j < matrix._contacts.size());
      matrix._contacts[j] = this->_chunk._count_buff[i];
    }
  }

  assert(matrix.get_n_of_missed_updates() == 0);
  matrix._updates_missed = this->_header.updates_missed;
  if constexpr (utils::ndebug_not_defined()) {
    matrix._global_stats_outdated = true;
    matrix.unsafe_update_global_stats();
    if constexpr (std::is_integral_v<SumT>) {
      assert(matrix._tot_contacts == this->_header.tot_contacts);
    }
    assert(matrix._nnz == this->_header.nnz);
  }

  matrix._tot_contacts = this->_header.tot_contacts;
  matrix._nnz = this->_header.nnz;

  this->reset();
}

template <class N>
void ContactMatrixSerde<N>::read_internal(std::istream& stream, ContactMatrixSparse<N>& matrix) {
  stream.seekg(std::ios::beg);
  this->read_header(stream);

  matrix = ContactMatrixSparse<N>(this->_header.nrows, this->_header.ncols);

  for (const auto& chunk_metadata : this->_header.chunk_metadata.data()) {
    assert(chunk_metadata.offset != -1);
    this->read_chunk(stream, chunk_metadata.offset);
    assert(chunk_metadata.nnz == this->_chunk.size());
    assert(chunk_metadata.sum == this->_chunk.sum_contacts());

    for (usize i = 0; i < this->_chunk.size(); ++i) {
      const auto [rowt, colt] = internal::decode_idx(this->_chunk._idx_buff[i], matrix.nrows());
      assert(colt >= rowt);
      const auto row = colt - rowt;
      const auto col = colt;
      matrix.set(row, col, this->_chunk._count_buff[i]);
    }
  }

  assert(matrix.get_n_of_missed_updates() == 0);
  matrix._updates_missed = this->_header.updates_missed;
  if constexpr (utils::ndebug_not_defined()) {
    matrix._global_stats_outdated = true;
    matrix.update_global_stats(matrix.lock_tables());
    if constexpr (std::is_integral_v<SumT>) {
      assert(matrix._tot_contacts == this->_header.tot_contacts);
    }
    assert(matrix._nnz == this->_header.nnz);
  }

  matrix._tot_contacts = this->_header.tot_contacts;
  matrix._nnz = this->_header.nnz;

  this->reset();
}

template <class N>
void ContactMatrixSerde<N>::write_internal(std::ostream& stream,
                                           const ContactMatrixDense<N>& matrix,
                                           std::string_view metadata) {
  stream.seekp(std::ios::beg);
  this->_header = Header{matrix, metadata};
  this->write_empty_header(stream);

  usize j = 0;
  for (usize i = 0; i < this->_header.chunk_metadata.capacity(); ++i) {
    this->_chunk.clear();
    const auto chunk_start = j;
    ChunkMetadata m{stream.tellp(), 0, 0};

    for (; !this->_chunk.full() && j < matrix.npixels(); ++j) {
      if (const auto count = matrix._contacts[j]; count != 0) {
        this->_chunk.push_back(j, count);
        m.sum += count;
        m.nnz++;
      }
    }
    const auto chunk_end = j;

    this->_header.chunk_metadata.insert(chunk_start, chunk_end, m);
    this->write_chunk(stream);
  }

  assert(this->_header.chunk_metadata.size() == this->_header.num_chunks());
  this->write_header(stream);
  this->reset();
}

template <class N>
void ContactMatrixSerde<N>::write_internal(std::ostream& stream,
                                           const ContactMatrixSparse<N>& matrix,
                                           std::string_view metadata) {
  stream.seekp(std::ios::beg);

  const auto tables = matrix.lock_tables();

  this->_header = Header{matrix, metadata};
  this->write_empty_header(stream);

  usize chunk_start = 0;
  usize chunk_end = 0;
  for (const auto& block : tables) {
    chunk_start = chunk_end;
    ChunkMetadata m{stream.tellp(), 0, 0};

    this->_chunk.reserve(block.size());
    this->_chunk.clear();

    for (const auto& [idx, count] : block) {
      assert(idx < matrix.npixels());
      assert(count != 0);
      this->_chunk.push_back(idx, count);
      if constexpr (std::is_integral_v<N>) {
        m.sum += count;
      }
      m.nnz++;
    }

    this->_chunk.sort(this->_sorting_buff);
    chunk_end = this->_chunk._idx_buff.back() + 1;
    if constexpr (std::is_floating_point_v<N>) {
      // Summing FP contacts afrer sorting should yield reproducible results
      m.sum = this->_chunk.sum_contacts();
    }

    this->_header.chunk_metadata.insert(chunk_start, chunk_end, m);
    this->write_chunk(stream);
  }

  assert(this->_header.chunk_metadata.size() == this->_header.num_chunks());
  this->write_header(stream);
  this->reset();
}

template <class N>
void ContactMatrixSerde<N>::reset() noexcept {
  this->_header.clear();
  this->_chunk.clear();
  this->_sorting_buff.clear();
}

template <class N>
void ContactMatrixSerde<N>::read_chunk(std::istream& stream, i64 offset) {
  namespace bio = boost::iostreams;

  if (offset >= 0) {
    stream.seekg(offset);
  }

  bio::filtering_istreambuf fis{};
  fis.push(bio::zstd_decompressor());
  fis.push(stream);

  boost::archive::binary_iarchive bi(fis);
  bi >> this->_chunk;
}

template <class N>
void ContactMatrixSerde<N>::write_internal_destructive(std::ostream& stream,
                                                       ContactMatrixDense<N>& matrix,
                                                       std::string_view metadata) {
  this->write_internal(stream, matrix, metadata);
  ContactMatrixDense<N> empty_matrix{};
  std::swap(matrix, empty_matrix);
}

template <class N>
void ContactMatrixSerde<N>::write_internal_destructive(std::ostream& stream,
                                                       ContactMatrixSparse<N>& matrix,
                                                       std::string_view metadata) {
  using BlockT = typename ContactMatrixSparse<N>::BlockT;

  stream.seekp(std::ios::beg);

  auto tables = matrix.lock_tables();
  std::reverse(tables.begin(), tables.end());

  this->_header = Header{matrix, metadata};
  this->write_empty_header(stream);

  usize chunk_start = 0;
  usize chunk_end = 0;
  const auto ntables = tables.size();
  for (usize i = 0; i < ntables; ++i) {
    chunk_start = chunk_end;
    ChunkMetadata m{stream.tellp(), 0, 0};

    this->_chunk.reserve(tables.back().size());
    this->_chunk.clear();

    for (const auto& [idx, count] : tables.back()) {
      assert(idx < matrix.npixels());
      assert(count != 0);
      this->_chunk.push_back(idx, count);
      if constexpr (std::is_integral_v<N>) {
        m.sum += count;
      }
      m.nnz++;
    }

    {
      tables.pop_back();
      BlockT empty_block{0};
      std::swap(matrix._contact_blocks[i], empty_block);
    }

    this->_chunk.sort(this->_sorting_buff);
    chunk_end = this->_chunk._idx_buff.back() + 1;
    if constexpr (std::is_floating_point_v<N>) {
      // Summing FP contacts afrer sorting should yield reproducible results
      m.sum = this->_chunk.sum_contacts();
    }

    this->_header.chunk_metadata.insert(chunk_start, chunk_end, m);
    this->write_chunk(stream);
  }

  matrix = ContactMatrixSparse<N>{};

  assert(this->_header.chunk_metadata.size() == this->_header.num_chunks());
  this->write_header(stream);
  this->reset();
}

template <class N>
void ContactMatrixSerde<N>::write_chunk(std::ostream& stream, i64 offset) {
  namespace bio = boost::iostreams;

  if (offset >= 0) {
    stream.seekp(offset);
  }

  bio::filtering_ostreambuf fos{};
  fos.push(bio::zstd_compressor(bio::zstd_params(this->_compression_lvl)));
  fos.push(stream);

  boost::archive::binary_oarchive bo(fos);
  bo << this->_chunk;
}

}  // namespace modle
