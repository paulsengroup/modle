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

#include <array>    // for array
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
IITree<I, T>::Interval::Interval(I start_, I end_, const T &data_)
    : _start(start_), _end(end_), _max(_end), _data(data_){};

template <typename I, typename T>
IITree<I, T>::Interval::Interval(I start_, I end_, T &&data_)
    : _start(start_), _end(end_), _max(_end), _data(data_){};

template <typename I, typename T>
constexpr bool IITree<I, T>::Interval::operator<(
    const IITree<I, T>::Interval &other) const noexcept {
  if (this->_start != other._start) {
    return this->_start < other._start;
  }
  return this->_end < other._end;
}

template <typename I, typename T>
constexpr I IITree<I, T>::Interval::start() const noexcept {
  return this->_start;
}

template <typename I, typename T>
constexpr I IITree<I, T>::Interval::end() const noexcept {
  return this->_end;
}

template <typename I, typename T>
constexpr T *IITree<I, T>::Interval::data() const noexcept {
  return &this->_data;
}

template <typename I, typename T>
void IITree<I, T>::insert(const I start, const I end, const T &data) {
  this->_indexed = false;
  this->_data.emplace_back(start, end, data);
}

template <typename I, typename T>
void IITree<I, T>::emplace(const I start, const I end, T &&data) {
  this->_indexed = false;
  this->_data.emplace_back(start, end, std::move(data));
}

template <typename I, typename T>
void IITree<I, T>::make_BST() {
  cppsort::pdq_sort(this->_data.begin(), this->_data.end());
  this->_max_level = [&]() {
    if (this->empty()) {
      return -1LL;
    }

    for (auto i = 0UL; i < this->size(); i += 2) {  // leaves (i.e. at level 0)
      this->_data[i]._max = this->_data[i]._end;
    }

    // last_i points to the rightmost node in the tree
    auto [last_i, last] = [&]() {
      auto i = this->size() - 1;
      i -= i % 2;
      return std::make_pair(i, static_cast<I>(this->_data[i]._end));
    }();

    auto k = 1LL;
    const auto k1 = static_cast<int64_t>(this->size());
    for (; 1LL << k <= k1; ++k) {  // process internal nodes in the bottom-up order
      const auto x = 1ULL << (k - 1);
      const auto i0 = (x << 1) - 1;
      const auto step = x << 2;

      for (auto i = i0; i < this->size(); i += step) {  // i0 is the first node
        // traverse all nodes at level level
        I el = this->_data[i - x]._max;  // _max value of the left child
        I er = i + x < this->size() ? this->_data[i + x]._max : last;  // of the right child
        I e = this->_data[i]._end;
        e = e > el ? e : el;
        e = e > er ? e : er;
        this->_data[i]._max = e;  // set the _max value for node i
      }
      last_i = last_i >> k & 1
                   ? last_i - x
                   : last_i + x;  // last_i now points to the parent of the original last_i
      if (last_i < this->size() && this->_data[last_i]._max > last)  // update last accordingly
        last = this->_data[last_i]._max;
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
      for (auto i = i0; i < i1 && this->_data[i]._start < end; ++i)
        if (start < this->_data[i]._end) {  // if find_overlaps, append to overlapping_intervals[]
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
      if (lchild_idx >= this->size() || this->_data[lchild_idx]._max > start) {
        // push the left child if lchild_idx is overlapping_intervals of range or may find_overlaps
        // with the query
        stack[t++] = StackCell{cell.level - 1, lchild_idx, false};
      }
    } else if (cell.node_idx < this->size() && this->_data[cell.node_idx]._start < end) {
      // need to push the right child
      // test if cell.node_idx overlaps the query; if yes, append to overlapping_intervals[]
      if (start < this->_data[cell.node_idx]._end) {
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
      for (auto i = i0; i < i1 && this->_data[i]._start < end; ++i)
        if (start < this->_data[i]._end) {  // found an overlap
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
      if (lchild_idx >= this->size() || this->_data[lchild_idx]._max > start) {
        // push the left child if lchild_idx is out of range or may overlap with the query
        assert(t < stack.size());  // NOLINT
        stack[t++] = StackCell{cell.level - 1, lchild_idx, false};
      }
    } else if (cell.node_idx < this->size() && this->_data[cell.node_idx]._start < end) {
      // need to push the right child
      if (start < this->_data[cell.node_idx]._end) {
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
  return this->_data.size();
}

template <typename I, typename T>
constexpr bool IITree<I, T>::is_BST() const noexcept {
  return this->_indexed;
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
  return this->_data.reserve(new_capacity);
}

template <typename I, typename T>
I IITree<I, T>::get_overlap_start(const size_t i) const {
  return this->_data[i]._start;
}

template <typename I, typename T>
I IITree<I, T>::get_overlap_end(const size_t i) const {
  return this->_data[i]._end;
}

template <typename I, typename T>
const T& IITree<I, T>::get_overlap_data(const size_t i) const {
  return this->_data[i]._data;
}

template <typename I, typename T>
typename IITree<I, T>::iterator IITree<I, T>::begin() {
  return this->_data.begin();
}
template <typename I, typename T>
typename IITree<I, T>::iterator IITree<I, T>::end() {
  return this->_data.end();
}

template <typename I, typename T>
typename IITree<I, T>::const_iterator IITree<I, T>::cbegin() const {
  return this->_data.cbegin();
}
template <typename I, typename T>
typename IITree<I, T>::const_iterator IITree<I, T>::cend() const {
  return this->_data.cend();
}

}  // namespace modle
