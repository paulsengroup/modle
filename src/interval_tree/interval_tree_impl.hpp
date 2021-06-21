//   The MIT License
//   Copyright (c) 2019     Dana-Farber Cancer Institute
//   Copyright (c) 2021     modified by Roberto Rossini
//   Permission is hereby granted, free of charge, to any person obtaining
//   a copy of this software and associated documentation files (the
//   "Software"), to deal in the Software without restriction, including
//   without limitation the rights to use, copy, modify, merge, publish,
//   distribute, sublicense, and/or sell copies of the Software, and to
//   permit persons to whom the Software is furnished to do so, subject to
//   the following conditions:
//   The above copyright notice and this permission notice shall be
//   included in all copies or substantial portions of the Software.
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
//   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
//   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//   SOFTWARE.

#pragma once
#include <cpp-sort/sorters/pdq_sorter.h>  // for pdq_sorter
#include <fmt/format.h>

#include <array>    // for array
#include <cassert>  // for assert
#include <cstddef>  // for size_t
#include <stack>    // for stack
#include <vector>   // for vector

namespace modle {

template <typename I, typename T>
IITree<I, T>::StackCell::StackCell(int64_t _level, size_t _node_idx,
                                   bool _left_child_alread_processed)
    : level(_level),
      node_idx(_node_idx),
      left_child_already_processed(_left_child_alread_processed){};

template <typename I, typename T>
IITree<I, T>::IITree(size_t start_pos_, size_t end_pos_)
    : _start_pos(start_pos_), _end_pos(end_pos_) {
  assert(start_pos_ <= end_pos_);  // NOLINT
}

template <typename I, typename T>
template <typename I2, typename>
void IITree<I, T>::insert(const I2 start, const I2 end, const T &data) {
  assert(this->start_pos() <= start);  // NOLINT
  assert(this->end_pos() > end);       // NOLINT

  this->_indexed = false;
  this->_start.emplace_back(I2(start));
  this->_end.emplace_back(I2(end));
  this->_max.emplace_back(I2(end));
  this->_data.emplace_back(data);
}

template <typename I, typename T>
template <typename I2, typename>
void IITree<I, T>::emplace(const I2 start, const I2 end, T &&data) {
  assert(this->start_pos() <= start);  // NOLINT
  assert(this->end_pos() > end);       // NOLINT

  this->_indexed = false;
  this->_start.emplace_back(I2(start));
  this->_end.emplace_back(I2(end));
  this->_max.emplace_back(I2(end));
  this->_data.emplace_back(std::move(data));
}

template <typename I, typename T>
void IITree<I, T>::make_BST() {
  if (this->_data.empty() || this->_indexed) {
    return;
  }
  std::vector<size_t> ranks(this->size());
  std::iota(ranks.begin(), ranks.end(), 0UL);
  cppsort::pdq_sort(ranks.begin(), ranks.end(), [this](const auto i1, const auto i2) {
    assert(i1 < this->size());  // NOLINT
    assert(i2 < this->size());  // NOLINT
    return this->_start[i1] < this->_start[i2];
  });

  // https://stackoverflow.com/a/22218699
  for (auto i = 0UL; i < this->size(); i++) {
    while (ranks[i] != i) {
      const auto &j = ranks[i];
      const auto &k = ranks[j];
      std::swap(this->_start[j], this->_start[k]);
      std::swap(this->_end[j], this->_end[k]);
      std::swap(this->_max[j], this->_max[k]);
      std::swap(this->_data[j], this->_data[k]);

      std::swap(ranks[i], ranks[j]);
    }
  }

  // assert(std::is_sorted(this->_start.begin(), this->_start.end()));  // NOLINT

  this->_max_level = [&]() {
    if (this->empty()) {
      return -1LL;
    }

    for (auto i = 0UL; i < this->size(); i += 2) {  // leaves (i.e. at level 0)
      this->_max[i] = this->_end[i];
    }

    // last_i points to the rightmost node in the tree
    auto [last_i, last] = [&]() {
      auto i = this->size() - 1;
      i -= i % 2;
      return std::make_pair(i, static_cast<I>(this->_end[i]));
    }();

    auto k = 1LL;
    const auto k1 = static_cast<int64_t>(this->size());
    for (; 1LL << k <= k1; ++k) {  // process internal nodes in the bottom-up order
      const auto x = 1ULL << (k - 1);
      const auto i0 = (x << 1) - 1;
      const auto step = x << 2;

      for (auto i = i0; i < this->size(); i += step) {  // i0 is the first node
        // traverse all nodes at level level
        I el = this->_max[i - x];                                // _max value of the left child
        I er = i + x < this->size() ? this->_max[i + x] : last;  // of the right child
        I e = this->_end[i];
        e = e > el ? e : el;
        e = e > er ? e : er;
        this->_max[i] = e;  // set the _max value for node i
      }
      last_i = last_i >> k & 1
                   ? last_i - x
                   : last_i + x;  // last_i now points to the parent of the original last_i
      if (last_i < this->size() && this->_max[last_i] > last)  // update last accordingly
        last = this->_max[last_i];
    }
    return k - 1LL;
  }();

  this->_indexed = true;
}

template <typename I, typename T>
bool IITree<I, T>::find_overlaps(const I start, const I end,
                                 std::vector<size_t> &overlapping_intervals) const {
  std::array<StackCell, 64> stack{};  // NOLINT

  return this->find_overlaps(start, end, overlapping_intervals, absl::MakeSpan(stack));
}

template <typename I, typename T>
bool IITree<I, T>::find_overlaps(const I start, const I end,
                                 std::vector<size_t> &overlapping_intervals) {
  return this->find_overlaps(start, end, overlapping_intervals, absl::MakeSpan(this->_stack));
}

template <typename I, typename T>
bool IITree<I, T>::find_overlaps(const I start, const I end,
                                 std::vector<size_t> &overlapping_intervals,
                                 const absl::Span<StackCell> stack) const {
  assert(start <= end);  // NOLINT
  overlapping_intervals.clear();
  if (this->_max_level < 0) {
    return false;
  }

  assert(this->_indexed);  // NOLINT

  // push the root; this is a top down traversal
  stack.front() = StackCell{this->_max_level, (1ULL << this->_max_level) - 1, false};

  for (auto t = 1UL; t;) {  // the following guarantees that numbers in
                            // overlapping_intervals[] are always sorted
    const auto cell = stack[--t];
    if (cell.level < 4) {  // we are in a small subtree; traverse every node in this subtree
      const auto i0 = cell.node_idx >> cell.level << cell.level;
      const auto i1 = std::min(i0 + (1UL << (cell.level + 1)) - 1, this->size());
      for (auto i = i0; i < i1 && this->_start[i] < end; ++i)
        if (start < this->_end[i]) {  // if find_overlaps, append to overlapping_intervals[]
          overlapping_intervals.push_back(i);
        }
      return !overlapping_intervals.empty();
    }

    if (!cell.left_child_already_processed) {
      const auto lchild_idx =
          cell.node_idx -
          (1ULL << (cell.level - 1));  // the left child of cell.node_idx
                                       // NB: lchild_idx may be overlapping_intervals of range
                                       // (i.e. lchild_idx >= _data.size())

      // re-add node cell.node_idx, but mark the left child as having been processed
      stack[t++] = StackCell{cell.level, cell.node_idx, true};
      if (lchild_idx >= this->size() || this->_max[lchild_idx] > start) {
        // push the left child if lchild_idx is overlapping_intervals of range or may
        // find_overlaps with the query
        stack[t++] = StackCell{cell.level - 1, lchild_idx, false};
      }
    } else if (cell.node_idx < this->size() && this->_start[cell.node_idx] < end) {
      // need to push the right child
      // test if cell.node_idx overlaps the query; if yes, append to overlapping_intervals[]
      if (start < this->_end[cell.node_idx]) {
        overlapping_intervals.push_back(cell.node_idx);
      }
      stack[t++] = StackCell{cell.level - 1, cell.node_idx + (1LL << (cell.level - 1)),
                             false};  // push the right child
    }
  }
  return !overlapping_intervals.empty();
}

template <typename I, typename T>
bool IITree<I, T>::overlaps_with(const I start, const I end) noexcept {
  return this->overlaps_with(start, end, absl::MakeSpan(this->_stack));
}

template <typename I, typename T>
bool IITree<I, T>::overlaps_with(const I start, const I end) const {
  std::array<StackCell, 64> stack;  // NOLINT
  return this->overlaps_with(start, end, absl::MakeSpan(stack));
}

template <typename I, typename T>
bool IITree<I, T>::overlaps_with(const I start, const I end,
                                 const absl::Span<StackCell> stack) const noexcept {
  assert(start <= end);  // NOLINT
  if (this->_max_level < 0) {
    return false;
  }

  assert(this->_indexed);  // NOLINT

  // push the root; this is a top down traversal
  stack.front() = StackCell{this->_max_level, (1ULL << this->_max_level) - 1, false};

  for (auto t = 1UL; t;) {  // the following guarantees that numbers in out[] are always sorted
    const auto cell = stack[--t];
    if (cell.level < 4) {  // we are in a small subtree; traverse every node in this subtree
      const auto i0 = cell.node_idx >> cell.level << cell.level;
      const auto i1 = std::min(i0 + (1UL << (cell.level + 1)) - 1, this->size());
      for (auto i = i0; i < i1 && this->_start[i] < end; ++i)
        if (start < this->_end[i]) {  // found an overlap
          return true;
        }
      return false;
    }

    if (!cell.left_child_already_processed) {  // if left child not processed
      const auto lchild_idx =
          cell.node_idx - (1ULL << (cell.level - 1));  // the left child of cell.node_idx
                                                       // NB: lchild_idx may be out of range
                                                       // (i.e. lchild_idx >= _data.size())

      // re-add node cell.node_idx, but mark the left child as having been processed
      stack[t++] = StackCell{cell.level, cell.node_idx, true};
      if (lchild_idx >= this->size() || this->_max[lchild_idx] > start) {
        // push the left child if lchild_idx is out of range or may overlap with the query
        assert(t < stack.size());  // NOLINT
        stack[t++] = StackCell{cell.level - 1, lchild_idx, false};
      }
    } else if (cell.node_idx < this->size() && this->_start[cell.node_idx] < end) {
      // need to push the right child
      if (start < this->_end[cell.node_idx]) {
        return true;  // test if cell.node_idx overlaps the query; if yes, append to out[]
      }
      stack[t++] = StackCell{cell.level - 1, cell.node_idx + (1LL << (cell.level - 1)),
                             false};  // push the right child
    }
  }
  return false;
}

template <typename I, typename T>
constexpr size_t IITree<I, T>::size() const noexcept {
  assert(this->_data.size() == this->_start.size());  // NOLINT
  assert(this->_data.size() == this->_end.size());    // NOLINT
  assert(this->_data.size() == this->_max.size());    // NOLINT
  return this->_data.size();
}

template <typename I, typename T>
constexpr bool IITree<I, T>::is_BST() const noexcept {
  return this->_indexed;
}

template <typename I, typename T>
constexpr size_t IITree<I, T>::span() const noexcept {
  assert(this->start_pos() <= this->end_pos());  // NOLINT
  return this->end_pos() - this->start_pos();
}

template <typename I, typename T>
constexpr size_t IITree<I, T>::start_pos() const noexcept {
  return this->_start_pos;
}

template <typename I, typename T>
constexpr size_t IITree<I, T>::end_pos() const noexcept {
  return this->_end_pos;
}

template <typename I, typename T>
constexpr size_t IITree<I, T>::capacity() const noexcept {
  return this->_data.capacity();
}

template <typename I, typename T>
constexpr bool IITree<I, T>::empty() const noexcept {
  return this->_data.empty();
}

template <typename I, typename T>
void IITree<I, T>::reserve(const size_t new_capacity) {
  this->_start.reserve(new_capacity);
  this->_end.reserve(new_capacity);
  this->_max.reserve(new_capacity);
  this->_data.reserve(new_capacity);
}

template <typename I, typename T>
void IITree<I, T>::clear() {
  this->_start.clear();
  this->_end.clear();
  this->_max.clear();
  this->_data.clear();

  this->_max_level = 0;
  this->_indexed = false;
}

template <typename I, typename T>
I IITree<I, T>::get_overlap_start(const size_t i) const {
  return this->_start[i];
}

template <typename I, typename T>
I IITree<I, T>::get_overlap_end(const size_t i) const {
  return this->_end[i];
}

template <typename I, typename T>
const T &IITree<I, T>::get_overlap_data(const size_t i) const {
  return this->_data[i];
}

template <typename I, typename T>
const absl::Span<const I> IITree<I, T>::starts() const {
  return this->_start;
}

template <typename I, typename T>
const absl::Span<const I> IITree<I, T>::ends() const {
  return this->_data;
}

template <typename I, typename T>
const absl::Span<const T> IITree<I, T>::data() const {
  return this->_data;
}

template <typename I, typename T>
typename IITree<I, T>::I_iterator IITree<I, T>::starts_begin() {
  return this->_start.begin();
}

template <typename I, typename T>
typename IITree<I, T>::I_iterator IITree<I, T>::starts_end() {
  return this->_start.end();
}

template <typename I, typename T>
typename IITree<I, T>::I_iterator IITree<I, T>::ends_begin() {
  return this->_end.begin();
}

template <typename I, typename T>
typename IITree<I, T>::I_iterator IITree<I, T>::ends_end() {
  return this->_end.end();
}

template <typename I, typename T>
typename IITree<I, T>::T_iterator IITree<I, T>::data_begin() {
  return this->_data.begin();
}

template <typename I, typename T>
typename IITree<I, T>::T_iterator IITree<I, T>::data_end() {
  return this->_data.end();
}

template <typename I, typename T>
typename IITree<I, T>::I_iterator_const IITree<I, T>::starts_begin() const {
  return this->_starts.begin();
}

template <typename I, typename T>
typename IITree<I, T>::I_iterator_const IITree<I, T>::starts_end() const {
  return this->_start.end();
}

template <typename I, typename T>
typename IITree<I, T>::I_iterator_const IITree<I, T>::ends_begin() const {
  return this->_end.begin();
}

template <typename I, typename T>
typename IITree<I, T>::I_iterator_const IITree<I, T>::ends_end() const {
  return this->_end.end();
}

template <typename I, typename T>
typename IITree<I, T>::T_iterator_const IITree<I, T>::data_begin() const {
  return this->_data.begin();
}

template <typename I, typename T>
typename IITree<I, T>::T_iterator_const IITree<I, T>::data_end() const {
  return this->_data.end();
}

}  // namespace modle
