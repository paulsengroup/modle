// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
// Copyright (C) 2019 Dana-Farber Cancer Institute
//
// SPDX-License-Identifier: MIT

#pragma once
#include <absl/types/span.h>
#include <cpp-sort/sorters/pdq_sorter.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"

namespace modle {

template <class N, class T>
IITree<N, T>::StackCell::StackCell(i64 _level, usize _node_idx, bool _left_child_alread_processed)
    : level(_level),
      node_idx(_node_idx),
      left_child_already_processed(_left_child_alread_processed) {}

template <class N, class T>
IITree<N, T>::IITree(N start_pos_, N end_pos_) : _start_pos(start_pos_), _end_pos(end_pos_) {
  assert(start_pos_ <= end_pos_);
}

template <class N, class T>
void IITree<N, T>::insert(const N start, const N end, const T &data) {
  assert(start_pos() <= start);
  assert(end_pos() >= end);

  _indexed = false;
  _start.emplace_back(N(start));
  _end.emplace_back(N(end));
  _max.emplace_back(N(end));
  _data.emplace_back(data);
}

template <class N, class T>
void IITree<N, T>::emplace(const N start, const N end, T &&data) {
  assert(start_pos() <= start);
  assert(end_pos() >= end);

  _indexed = false;
  _start.emplace_back(N(start));
  _end.emplace_back(N(end));
  _max.emplace_back(N(end));
  _data.emplace_back(std::move(data));
}

template <class N, class T>
void IITree<N, T>::make_BST() {
  if (_data.empty() || _indexed) {
    return;
  }
  std::vector<usize> ranks(size());
  std::iota(ranks.begin(), ranks.end(), 0);
  cppsort::pdq_sort(ranks.begin(), ranks.end(), [this](const auto i1, const auto i2) {
    assert(i1 < size());
    assert(i2 < size());
    return _start[i1] < _start[i2];
  });

  // https://stackoverflow.com/a/22218699
  for (usize i = 0; i < size(); i++) {
    while (ranks[i] != i) {
      const auto &j = ranks[i];
      const auto &k = ranks[j];
      std::swap(_start[j], _start[k]);
      std::swap(_end[j], _end[k]);
      std::swap(_max[j], _max[k]);
      std::swap(_data[j], _data[k]);

      std::swap(ranks[i], ranks[j]);
    }
  }

  // assert(std::is_sorted(_start.begin(), _start.end()));

  _max_level = [&]() {
    if (empty()) {
      return -1LL;
    }

    for (usize i = 0; i < size(); i += 2) {  // leaves (i.e. at level 0)
      _max[i] = _end[i];
    }

    // last_i points to the rightmost node in the tree
    auto [last_i, last] = [&]() {
      auto i = size() - 1;
      i -= i % 2;
      return std::make_pair(i, static_cast<N>(_end[i]));
    }();

    auto k = 1LL;
    const auto k1 = static_cast<i64>(size());
    for (; 1LL << k <= k1; ++k) {  // process internal nodes in the bottom-up order
      const auto x = static_cast<usize>(1L << (k - 1));
      const auto i0 = (x << 1) - 1;
      const auto step = x << 2;

      for (auto i = i0; i < size(); i += step) {  // i0 is the first node
        // traverse all nodes at level level
        N el = _max[i - x];                          // _max value of the left child
        N er = i + x < size() ? _max[i + x] : last;  // of the right child
        N e = _end[i];
        e = e > el ? e : el;
        e = e > er ? e : er;
        _max[i] = e;  // set the _max value for node i
      }
      last_i = last_i >> k & 1
                   ? last_i - x
                   : last_i + x;  // last_i now points to the parent of the original last_i
      if (last_i < size() && _max[last_i] > last) {  // update last accordingly
        last = _max[last_i];
      }
    }
    return k - 1LL;
  }();

  _indexed = true;
}

// GCC chokes if we use an using directive here
#define T_iterator       typename IITree<N, T>::T_iterator
#define I_iterator       typename IITree<N, T>::I_iterator
#define T_iterator_const typename IITree<N, T>::T_iterator_const
#define I_iterator_const typename IITree<N, T>::I_iterator_const

template <class N, class T>
template <typename IITree<N, T>::QueryType Mode>
std::pair<T_iterator_const, T_iterator_const> IITree<N, T>::internal_equal_range(
    const N start, const N end) const noexcept {
  assert(start <= end);
  std::array<StackCell, 64> stack;
  if (_max_level < 0) {
    return std::make_pair(_data.end(), _data.end());
  }

  assert(_indexed);
  T_iterator_const overlap_begin = _data.end();
  T_iterator_const overlap_end = _data.begin();

  auto update_overlap_range = [&](auto idx) {  // The idea here is to have a single implementation
                                               // for lower/upper bound and equal range queries
    using ptr_difference_t = typename std::iterator_traits<T_iterator_const>::difference_type;
    if constexpr ((Mode & QueryType::LOWER_BOUND) == 0) {
      overlap_begin = std::min(overlap_begin, _data.begin() + static_cast<ptr_difference_t>(idx));
    }
    if constexpr ((Mode & QueryType::UPPER_BOUND) == 0) {
      overlap_end = std::max(overlap_end, _data.begin() + static_cast<ptr_difference_t>(idx + 1));
    }
    if constexpr ((Mode & QueryType::CONTAINS) || (Mode & QueryType::UPPER_BOUND) != 0) {
      return true;
    }
    return false;
  };

  // push the root; this is a top down traversal
  stack.front() = StackCell{_max_level, (1ULL << _max_level) - 1, false};

  // The following guarantees that we visit overlapping intervals in increasing order
  for (usize t = 1; t > 0;) {
    assert(t <= 64);
    const auto cell = stack[--t];
    if (cell.level < 4) {  // we are in a small subtree; traverse every node in this subtree
      const auto i0 = cell.node_idx >> cell.level << cell.level;
      const auto i1 = std::min(i0 + (1UL << (cell.level + 1)) - 1, size());
      for (auto i = i0; i < i1 && _start[i] < end; ++i) {
        if (start < _end[i]) {  // found an overlap
          if (update_overlap_range(i)) {
            goto return_stmt;
          }
        }
      }
    } else if (!cell.left_child_already_processed) {
      const auto lchild_idx =
          cell.node_idx -
          (1ULL << (cell.level - 1));  // the left child of cell.node_idx
                                       // NB: lchild_idx may be overlapping_intervals of range
                                       // (i.e. lchild_idx >= _data.size())

      // re-add node cell.node_idx, but mark the left child as having been processed
      stack[t++] = StackCell{cell.level, cell.node_idx, true};
      if (lchild_idx >= size() || _max[lchild_idx] > start) {
        // push the left child if lchild_idx is out of range or may overlap with the query
        stack[t++] = StackCell{cell.level - 1, lchild_idx, false};
      }
    } else if (cell.node_idx < size() && _start[cell.node_idx] < end) {
      // need to push the right child
      // test if cell.node_idx overlaps the query; if yes, append to overlapping_intervals[]
      if (start < _end[cell.node_idx]) {  // found an overlap
        if (update_overlap_range(cell.node_idx)) {
          goto return_stmt;
        }
      }
      assert(cell.level > 0);
      stack[t++] =
          StackCell{cell.level - 1, cell.node_idx + (1ULL << static_cast<usize>(cell.level - 1)),
                    false};  // push the right child
    }
  }

return_stmt:
  if (overlap_begin == _data.end()) {
    overlap_end = overlap_begin;
  }
  return std::make_pair(overlap_begin, overlap_end);
}

template <class N, class T>
std::pair<T_iterator_const, T_iterator_const> IITree<N, T>::equal_range(
    const N start, const N end) const noexcept {
  return internal_equal_range<QueryType::EQUAL_RANGE>(start, end);
}

template <class N, class T>
std::pair<T_iterator_const, T_iterator_const> IITree<N, T>::find_overlaps(
    const N start, const N end) const noexcept {
  return equal_range(start, end);
}

template <class N, class T>
T_iterator_const IITree<N, T>::lower_bound(const N start, const N end) const noexcept {
  const auto [overlap_begin, _] = internal_equal_range<QueryType::LOWER_BOUND>(start, end);
  return overlap_begin;
}

template <class N, class T>
T_iterator_const IITree<N, T>::upper_bound(const N start, const N end) const noexcept {
  const auto [_, overlap_end] = internal_equal_range<QueryType::UPPER_BOUND>(start, end);
  return overlap_end;
}

template <class N, class T>
std::pair<usize, usize> IITree<N, T>::equal_range_idx(const N start, const N end) const noexcept {
  const auto [overlap_begin, overlap_end] =
      internal_equal_range<QueryType::EQUAL_RANGE>(start, end);

  return std::make_pair(overlap_begin - _data.begin(), overlap_end - _data.begin());
}

template <class N, class T>
std::pair<usize, usize> IITree<N, T>::find_overlaps_idx(const N start, const N end) const noexcept {
  return equal_range_idx(start, end);
}

template <class N, class T>
usize IITree<N, T>::lower_bound_idx(const N start, const N end) const noexcept {
  [[maybe_unused]] const auto [overlap_begin, _] =
      internal_equal_range<QueryType::LOWER_BOUND>(start, end);
  return _data.begin() - overlap_begin;
}

template <class N, class T>
usize IITree<N, T>::upper_bound_idx(const N start, const N end) const noexcept {
  [[maybe_unused]] const auto [_, overlap_end] =
      internal_equal_range<QueryType::UPPER_BOUND>(start, end);
  return _data.begin() - overlap_end;
}

template <class N, class T>
bool IITree<N, T>::overlaps_with(const N start, const N end) const noexcept {
  const auto overlap_begin = internal_equal_range<QueryType::CONTAINS>(start, end).first;
  return overlap_begin != _data.end();
}

template <class N, class T>
usize IITree<N, T>::count(const N start, const N end) const noexcept {
  const auto [overlap_begin, overlap_end] =
      internal_equal_range<QueryType::EQUAL_RANGE>(start, end);
  if (overlap_begin == _data.end()) {
    return 0UL;
  }
  assert(overlap_begin <= overlap_end);
  return static_cast<usize>(overlap_end - overlap_begin);
}

template <class N, class T>
constexpr usize IITree<N, T>::size() const noexcept {
  assert(_data.size() == _start.size());
  assert(_data.size() == _end.size());
  assert(_data.size() == _max.size());
  return _data.size();
}

template <class N, class T>
constexpr bool IITree<N, T>::is_BST() const noexcept {
  return _indexed;
}

template <class N, class T>
constexpr usize IITree<N, T>::span() const noexcept {
  assert(start_pos() <= end_pos());
  return end_pos() - start_pos();
}

template <class N, class T>
constexpr N IITree<N, T>::start_pos() const noexcept {
  return _start_pos;
}

template <class N, class T>
constexpr N IITree<N, T>::end_pos() const noexcept {
  return _end_pos;
}

template <class N, class T>
constexpr usize IITree<N, T>::capacity() const noexcept {
  return _data.capacity();
}

template <class N, class T>
constexpr bool IITree<N, T>::empty() const noexcept {
  return _data.empty();
}

template <class N, class T>
void IITree<N, T>::reserve(const usize new_capacity) {
  _start.reserve(new_capacity);
  _end.reserve(new_capacity);
  _max.reserve(new_capacity);
  _data.reserve(new_capacity);
}

template <class N, class T>
void IITree<N, T>::clear() {
  _start.clear();
  _end.clear();
  _max.clear();
  _data.clear();

  _max_level = 0;
  _indexed = false;
}

template <class N, class T>
N IITree<N, T>::get_overlap_start(const usize i) const {
  return _start[i];
}

template <class N, class T>
N IITree<N, T>::get_overlap_end(const usize i) const {
  return _end[i];
}

template <class N, class T>
const T &IITree<N, T>::get_overlap_data(const usize i) const {
  return _data[i];
}

template <class N, class T>
const absl::Span<const N> IITree<N, T>::starts() const {
  return _start;
}

template <class N, class T>
const absl::Span<const N> IITree<N, T>::ends() const {
  return _data;
}

template <class N, class T>
const absl::Span<const T> IITree<N, T>::data() const {
  return _data;
}

template <class N, class T>
I_iterator IITree<N, T>::starts_begin() {
  return _start.begin();
}

template <class N, class T>
I_iterator IITree<N, T>::starts_end() {
  return _start.end();
}

template <class N, class T>
I_iterator IITree<N, T>::ends_begin() {
  return _end.begin();
}

template <class N, class T>
I_iterator IITree<N, T>::ends_end() {
  return _end.end();
}

template <class N, class T>
T_iterator IITree<N, T>::data_begin() {
  return _data.begin();
}

template <class N, class T>
T_iterator IITree<N, T>::data_end() {
  return _data.end();
}

template <class N, class T>
I_iterator_const IITree<N, T>::starts_begin() const {
  return _start.begin();
}

template <class N, class T>
I_iterator_const IITree<N, T>::starts_end() const {
  return _start.end();
}

template <class N, class T>
I_iterator_const IITree<N, T>::ends_begin() const {
  return _end.begin();
}

template <class N, class T>
I_iterator_const IITree<N, T>::ends_end() const {
  return _end.end();
}

template <class N, class T>
T_iterator_const IITree<N, T>::data_begin() const {
  return _data.begin();
}

template <class N, class T>
T_iterator_const IITree<N, T>::data_end() const {
  return _data.end();
}

#undef T_iterator
#undef I_iterator
#undef T_iterator_const
#undef I_iterator_const

}  // namespace modle
