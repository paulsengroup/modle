// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/time/clock.h>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <fstream>
#include <numeric>
#include <ostream>

#include "modle/common/common.hpp"
#include "modle/internal/contact_matrix_internal.hpp"

namespace modle {

template <class N>
ContactMatrixSparse<N>::ContactMatrixSparse(const ContactMatrixSparse<N>& other)
    : _nrows(other.nrows()),
      _ncols(other.ncols()),
      _cols_per_chunk(other._cols_per_chunk),
      // _contact_blocks(other._contact_blocks),
      _tot_contacts(other._tot_contacts.load()),
      _nnz(other._nnz.load()),
      _global_stats_outdated(other._global_stats_outdated.load()),
      _updates_missed(other._updates_missed.load()) {
  // Temporary workaround for https://github.com/efficient/libcuckoo/issues/148
  // All of this is just to avoid calling libcuckoo's map and bucker_container copy
  // constructor/operator
  const auto tables = other.lock_tables();
  _contact_blocks.resize(tables.size());

  for (usize i = 0; i < _contact_blocks.size(); ++i) {
    const auto& table = tables[i];
    auto& map = _contact_blocks[i];
    map.reserve(table.size());
    for (const auto& [k, v] : table) {
      map.template insert(k, v);
    }
  }
}

template <class N>
ContactMatrixSparse<N>::ContactMatrixSparse(usize nrows, usize ncols, ChunkSize max_chunk_size)
    : _nrows(std::min(nrows, ncols)),
      _ncols(ncols),
      _cols_per_chunk(compute_cols_per_chunk(_nrows, max_chunk_size)),
      _contact_blocks(compute_num_chunks(_ncols, _cols_per_chunk)) {}

template <class N>
ContactMatrixSparse<N>::ContactMatrixSparse(bp_t length, bp_t diagonal_width, bp_t bin_size,
                                            ChunkSize max_chunk_size)
    : ContactMatrixSparse((diagonal_width + bin_size - 1) / bin_size,
                          (length + bin_size - 1) / bin_size, max_chunk_size) {}

template <class N>
ContactMatrixSparse<N>& ContactMatrixSparse<N>::operator=(const ContactMatrixSparse<N>& other) {
  if (this == &other) {
    return *this;
  }

  _nrows = other.nrows();
  _ncols = other.ncols();
  _cols_per_chunk = other._cols_per_chunk;
  // _contact_blocks = other._contact_blocks;
  _tot_contacts = other._tot_contacts.load();
  _nnz = other._nnz.load();
  _global_stats_outdated = other._global_stats_outdated.load();
  _updates_missed = other._updates_missed.load();

  // Temporary workaround for https://github.com/efficient/libcuckoo/issues/148
  // All of this is just to avoid calling libcuckoo's map and bucker_container copy
  // constructor/operator
  const auto tables = other.lock_tables();
  _contact_blocks.resize(tables.size());

  for (usize i = 0; i < _contact_blocks.size(); ++i) {
    const auto& table = tables[i];
    auto& map = _contact_blocks[i];
    map.reserve(table.size());
    for (const auto& [k, v] : table) {
      map.template insert(k, v);
    }
  }

  return *this;
}

template <class N>
constexpr usize ContactMatrixSparse<N>::compute_cols_per_chunk(usize nrows,
                                                               ChunkSize max_chunk_size) {
  return std::max(usize(1), max_chunk_size.size / nrows);
}

template <class N>
constexpr usize ContactMatrixSparse<N>::compute_num_chunks(usize ncols, usize cols_per_chunk) {
  return (ncols + cols_per_chunk - 1) / cols_per_chunk;
}

template <class N>
usize ContactMatrixSparse<N>::compute_block_idx(usize col) const noexcept {
  return col / this->_cols_per_chunk;
}

template <class N>
usize ContactMatrixSparse<N>::num_blocks() const noexcept {
  return this->_contact_blocks.size();
}

template <class N>
auto ContactMatrixSparse<N>::get_block(usize col) noexcept -> BlockT& {
  const auto i = this->compute_block_idx(col);
  assert(i < this->num_blocks());
  return this->_contact_blocks[i];
}

template <class N>
auto ContactMatrixSparse<N>::get_block(usize col) const noexcept -> const BlockT& {
  const auto i = this->compute_block_idx(col);
  assert(i < this->num_blocks());
  return this->_contact_blocks[i];
}

template <class N>
N ContactMatrixSparse<N>::get(usize row, usize col) const {
  const auto [rowt, colt] = internal::transpose_coords(row, col);
  this->bound_check_coords(rowt, colt);

  if (rowt > this->nrows()) {
    return 0;
  }

  const auto i = internal::encode_idx(rowt, colt, this->_nrows);

  N n{0};
  auto& block = this->get_block(colt);
  block.find(i, n);
  return n;
}

template <class N>
void ContactMatrixSparse<N>::set(usize row, usize col, N n) {
  const auto [rowt, colt] = internal::transpose_coords(row, col);
  this->bound_check_coords(rowt, colt);

  if (rowt > this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, usize(1), std::memory_order_relaxed);
    return;
  }

  const auto i = internal::encode_idx(rowt, colt, this->_nrows);

  auto& block = this->get_block(colt);
  block.uprase_fn(
      i,
      [n](auto& num) {
        if (MODLE_UNLIKELY(n == 0)) {
          return true;
        }
        num = n;
        return false;
      },
      n);
  this->_global_stats_outdated = true;
}

template <class N>
void ContactMatrixSparse<N>::add(usize row, usize col, N n) {
  const auto [rowt, colt] = internal::transpose_coords(row, col);
  this->bound_check_coords(rowt, colt);

  if (rowt > this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, usize(1), std::memory_order_relaxed);
    return;
  }

  const auto i = internal::encode_idx(rowt, colt, this->_nrows);

  auto& block = this->get_block(colt);
  block.upsert(
      i, [n](auto& num) { num += n; }, n);
  this->_global_stats_outdated = true;
}

template <class N>
void ContactMatrixSparse<N>::subtract(const usize row, const usize col, const N n) {
  assert(n >= 0);
  const auto [rowt, colt] = internal::transpose_coords(row, col);
  this->bound_check_coords(rowt, colt);

  if (rowt > this->nrows()) {
    std::atomic_fetch_add_explicit(&this->_updates_missed, usize(1), std::memory_order_relaxed);
    return;
  }

  const auto i = internal::encode_idx(rowt, colt, this->_nrows);

  auto& block = this->get_block(colt);
  block.uprase_fn(
      i,
      [n](auto& num) {
        num -= n;
        if (MODLE_UNLIKELY(num == 0)) {
          return true;
        }
        return false;
      },
      n);
  this->_global_stats_outdated = true;
}

template <class N>
void ContactMatrixSparse<N>::increment(usize row, usize col) {
  this->add(row, col, 1);
}

template <class N>
void ContactMatrixSparse<N>::decrement(usize row, usize col) {
  this->subtract(row, col, N(1));
}

template <class N>
constexpr usize ContactMatrixSparse<N>::ncols() const {
  return this->_ncols;
}

template <class N>
constexpr usize ContactMatrixSparse<N>::nrows() const {
  return this->_nrows;
}

template <class N>
constexpr usize ContactMatrixSparse<N>::npixels() const {
  return this->_nrows * this->_ncols;
}

template <class N>
constexpr usize ContactMatrixSparse<N>::get_n_of_missed_updates() const noexcept {
  return this->_updates_missed.load();
}

template <class N>
double ContactMatrixSparse<N>::get_fraction_of_missed_updates() const {
  return static_cast<double>(this->get_n_of_missed_updates()) /
         utils::conditional_static_cast<double>(this->get_tot_contacts());
}

template <class N>
double ContactMatrixSparse<N>::unsafe_get_fraction_of_missed_updates() const {
  return static_cast<double>(this->get_n_of_missed_updates()) /
         utils::conditional_static_cast<double>(this->unsafe_get_tot_contacts());
}

template <class N>
double ContactMatrixSparse<N>::get_avg_contact_density() const {
  return utils::conditional_static_cast<double>(this->get_tot_contacts()) /
         static_cast<double>(this->npixels());
}

template <class N>
double ContactMatrixSparse<N>::unsafe_get_avg_contact_density() const {
  return utils::conditional_static_cast<double>(this->unsafe_get_tot_contacts()) /
         static_cast<double>(this->npixels());
}

template <class N>
void ContactMatrixSparse<N>::update_global_stats(
    const std::vector<typename BlockT::locked_table>& tables) const {
  if (!this->_global_stats_outdated) {
    return;
  }

  usize nnz{0};
  SumT tot_contacts{0};
  for (auto& table : tables) {
    tot_contacts +=
        std::accumulate(table.begin(), table.end(), SumT(0),
                        [](auto accumulator, const auto& it) { return accumulator + it.second; });
    nnz += table.size();
  }
  this->_nnz = nnz;
  this->_tot_contacts = tot_contacts;
  this->_global_stats_outdated = false;
}

template <class N>
auto ContactMatrixSparse<N>::unsafe_get_tot_contacts() const -> SumT {
  if (this->_global_stats_outdated) {
    this->update_global_stats(this->lock_tables());
  }
  return this->_tot_contacts.load();
}

template <class N>
auto ContactMatrixSparse<N>::get_tot_contacts() const -> SumT {
  auto table_locks = this->lock_tables();
  if (this->_global_stats_outdated) {
    this->update_global_stats(table_locks);
  }
  return this->_tot_contacts.load();
}

template <class N>
usize ContactMatrixSparse<N>::unsafe_get_nnz() const {
  if (this->_global_stats_outdated) {
    this->update_global_stats(this->lock_tables());
  }
  return this->_nnz.load();
}

template <class N>
usize ContactMatrixSparse<N>::get_nnz() const {
  auto table_locks = this->lock_tables();
  if (this->_global_stats_outdated) {
    this->update_global_stats(table_locks);
  }
  return this->_nnz.load();
}

template <class N>
auto ContactMatrixSparse<N>::lock_tables() const -> std::vector<typename BlockT::locked_table> {
  using LockT = typename BlockT::locked_table;
  std::vector<LockT> locks;
  std::transform(this->_contact_blocks.begin(), this->_contact_blocks.end(),
                 std::back_inserter(locks), [&](auto& blk) { return blk.lock_table(); });
  return locks;
}

template <class N>
void ContactMatrixSparse<N>::bound_check_coords([[maybe_unused]] const usize row,
                                                [[maybe_unused]] const usize col) const {
  if constexpr (utils::ndebug_not_defined()) {
    internal::bound_check_coords(*this, row, col);
  }
}

template <class N>
i64 ContactMatrixSparse<N>::serialize(const std::filesystem::path& path) const {
  std::ofstream ofs(path, std::ios::binary);
  ofs.exceptions(std::ios::failbit);
  return this->serialize(ofs);
}

template <class N>
i64 ContactMatrixSparse<N>::serialize(std::ostream& out_stream) const {
  using HeaderT = internal::SerializedContactMatrixHeader<SumT>;
  using BuffT = internal::SerializationTmpBuffers<N>;

  try {
    const auto t0 = absl::Now();
    const auto tables = this->lock_tables();

    HeaderT header{*this};
    header.chunk_offsets.resize(this->num_blocks());
    assert(tables.size() == header.chunk_offsets.size());

    // Reserve space for header
    header.serialize(out_stream);

    BuffT buff{};
    for (size_t i = 0; i < this->num_blocks(); ++i) {
      const auto& block = tables[i];
      buff.reserve(block.size());
      buff.clear();

      header.chunk_offsets[i] = out_stream.tellp();
      for (const auto& [idx, count] : block) {
        assert(idx < this->npixels());
        buff.push_back(idx, count);
      }

      buff.sort();
      buff.serialize(out_stream);
    }

    // Write header
    out_stream.seekp(std::ios::beg);
    header.serialize(out_stream);

    const auto t1 = absl::Now();
    spdlog::debug(FMT_STRING("Serialized {} contacts in {}"), header.nnz,
                  absl::FormatDuration(t1 - t0));

  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while serializing an instance of "
                               "ContactMatrixSparse: {}"),
                    e.what()));
  }

  return out_stream.seekp(std::ios::end).tellp();
}

template <class N>
ContactMatrixSparse<N> ContactMatrixSparse<N>::deserialize(const std::filesystem::path& path) {
  std::ifstream ifs(path, std::ios::binary);
  ifs.exceptions(std::ios::failbit);

  return deserialize(ifs);
}

template <class N>
ContactMatrixSparse<N> ContactMatrixSparse<N>::deserialize(std::istream& in_stream) {
  using HeaderT = internal::SerializedContactMatrixHeader<SumT>;
  using BuffT = internal::SerializationTmpBuffers<N>;

  const auto header = HeaderT::deserialize(in_stream);
  BuffT tmp_buff{};

  ContactMatrixSparse<N> m(header.nrows, header.ncols);

  for (const auto& chunk_offset : header.chunk_offsets) {
    in_stream.seekg(chunk_offset);
    BuffT::deserialize(in_stream, tmp_buff);
    tmp_buff.copy_to_buff(m);
  }

  assert(m.get_n_of_missed_updates() == 0);
  m._updates_missed = header.updates_missed;
  if constexpr (utils::ndebug_not_defined()) {
    m._global_stats_outdated = true;
    m.update_global_stats(m.lock_tables());
    assert(m._tot_contacts == header.tot_contacts);
    assert(m._nnz == header.nnz);
  }

  m._tot_contacts = header.tot_contacts;
  m._nnz = header.nnz;

  return m;
}

}  // namespace modle
