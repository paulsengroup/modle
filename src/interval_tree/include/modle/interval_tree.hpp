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

#include <array>        // for array
#include <cstddef>      // for size_t, std::ptrdiff_t
#include <cstdint>      // for uint64_t
#include <type_traits>  // for is_integral_v
#include <vector>       // for vector

namespace modle {

/* Suppose there are N=2^(K+1)-1 sorted numbers in an array a[]. They
 * implicitly form a complete binary tree of height K+1. We consider leaves to
 * be at level 0. The binary tree has the following properties:
 *
 * 1. The lowest k-1 bits of nodes at level k are all 1. The k-th bit is 0.
 *    The first node at level k is indexed by 2^k-1. The root of the tree is
 *    indexed by 2^K-1.
 *
 * 2. For a node x at level k, its left child is x-2^(k-1) and the right child
 *    is x+2^(k-1).
 *
 * 3. For a node x at level k, it is a left child if its (k+1)-th bit is 0. Its
 *    parent node is x+2^k. Similarly, if the (k+1)-th bit is 1, x is a right
 *    child and its parent is x-2^k.
 *
 * 4. For a node x at level k, there are 2^(k+1)-1 nodes in the subtree
 *    descending from x, including x. The left-most leaf is x&~(2^k-1) (masking
 *    the lowest k bits to 0).
 *
 * When numbers can't fill a complete binary tree, the parent of a node may not
 * be present in the array. The implementation here still mimics a complete
 * tree, though getting the special casing right is a little complex. There may
 * be alternative solutions.
 *
 * As a sorted array can be considered as a binary search tree, we can
 * implement an interval tree on top of the idea. We only need to record, for
 * each node, the maximum value in the subtree descending from the node.
 */

template <typename I, typename T>
class IITree {
  static_assert(std::is_integral_v<I>, "I should be an integral type.");
  struct StackCell;

  using I_iterator = typename std::vector<I>::iterator;
  using T_iterator = typename std::vector<T>::iterator;

  using I_iterator_const = typename std::vector<I>::const_iterator;
  using T_iterator_const = typename std::vector<T>::const_iterator;

 public:
  inline explicit IITree(size_t start_pos_ = 0,
                         size_t end_pos_ = std::numeric_limits<size_t>::max());

  template <typename I2, typename = std::enable_if_t<std::is_integral_v<I2>>>
  inline void insert(I2 start, I2 end, const T &data);
  template <typename I2, typename = std::enable_if_t<std::is_integral_v<I2>>>
  inline void emplace(I2 start, I2 end, T &&data);

  [[nodiscard]] inline bool overlaps_with(I start, I end) noexcept;
  [[nodiscard]] inline bool overlaps_with(I start, I end) const;
  inline bool find_overlaps(I start, I end, std::vector<size_t> &overlapping_intervals);
  inline bool find_overlaps(I start, I end, std::vector<size_t> &overlapping_intervals) const;

  [[nodiscard]] inline constexpr size_t capacity() const noexcept;
  [[nodiscard]] inline constexpr bool empty() const noexcept;
  [[nodiscard]] inline constexpr size_t size() const noexcept;
  [[nodiscard]] inline constexpr bool is_BST() const noexcept;

  [[nodiscard]] inline constexpr size_t span() const noexcept;
  [[nodiscard]] inline constexpr size_t start_pos() const noexcept;
  [[nodiscard]] inline constexpr size_t end_pos() const noexcept;

  [[nodiscard]] inline I get_overlap_start(size_t i) const;
  [[nodiscard]] inline I get_overlap_end(size_t i) const;
  [[nodiscard]] inline const T &get_overlap_data(size_t i) const;

  [[nodiscard]] inline const absl::Span<const I> starts() const;
  [[nodiscard]] inline const absl::Span<const I> ends() const;
  [[nodiscard]] inline const absl::Span<const T> data() const;

  [[nodiscard]] inline I_iterator starts_begin();
  [[nodiscard]] inline I_iterator starts_end();

  [[nodiscard]] inline I_iterator ends_begin();
  [[nodiscard]] inline I_iterator ends_end();

  [[nodiscard]] inline T_iterator data_begin();
  [[nodiscard]] inline T_iterator data_end();

  [[nodiscard]] inline I_iterator_const starts_begin() const;
  [[nodiscard]] inline I_iterator_const starts_end() const;

  [[nodiscard]] inline I_iterator_const ends_begin() const;
  [[nodiscard]] inline I_iterator_const ends_end() const;

  [[nodiscard]] inline T_iterator_const data_begin() const;
  [[nodiscard]] inline T_iterator_const data_end() const;

  inline void clear();

  inline void reserve(size_t new_capacity);
  inline void make_BST();

 private:
  struct StackCell {
    inline StackCell() = default;
    inline StackCell(int64_t _level, size_t _node_idx, bool _left_child_alread_processed);

    int64_t level;
    size_t node_idx;
    bool left_child_already_processed;
  };

  size_t _start_pos;
  size_t _end_pos;
  int64_t _max_level{0};
  bool _indexed{false};

  std::vector<I> _start;
  std::vector<I> _end;
  std::vector<I> _max;
  std::vector<T> _data;

  std::array<StackCell, 64> _stack{};  // NOLINT

  [[nodiscard]] inline bool overlaps_with(I start, I end,
                                          absl::Span<StackCell> stack) const noexcept;
  inline bool find_overlaps(I start, I end, std::vector<size_t> &overlapping_intervals,
                            absl::Span<StackCell> stack) const;
};
}  // namespace modle

#include "../../interval_tree_impl.hpp"
