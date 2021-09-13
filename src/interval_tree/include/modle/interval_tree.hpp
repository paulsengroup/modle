// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
// Copyright (c) 2019     Dana-Farber Cancer Institute
// Copyright (c) 2021     modified by Roberto Rossini
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span

#include <array>        // for array
#include <cstddef>      // for size_t
#include <cstdint>      // for int64_t, uint8_t
#include <iterator>     // for pair
#include <limits>       // for numeric_limits
#include <type_traits>  // for enable_if_t
#include <utility>      // for pair
#include <vector>       // for vector

#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_PADDED, DISABLE_W...

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

DISABLE_WARNING_PUSH
DISABLE_WARNING_PADDED
template <typename I, typename T>
class IITree {
  DISABLE_WARNING_POP
  static_assert(std::is_integral_v<I>, "I should be an integral type.");
  struct StackCell;

  using I_iterator = typename std::vector<I>::iterator;
  using T_iterator = typename std::vector<T>::iterator;

  using I_iterator_const = typename std::vector<I>::const_iterator;
  using T_iterator_const = typename std::vector<T>::const_iterator;

 public:
  inline explicit IITree(size_t start_pos_ = 0,
                         size_t end_pos_ = (std::numeric_limits<size_t>::max)());

  template <typename I2, typename = std::enable_if_t<std::is_integral_v<I2>>>
  inline void insert(I2 start, I2 end, const T &data);
  template <typename I2, typename = std::enable_if_t<std::is_integral_v<I2>>>
  inline void emplace(I2 start, I2 end, T &&data);

  [[nodiscard]] inline bool overlaps_with(I start, I end) const noexcept;
  [[nodiscard]] inline std::pair<T_iterator_const, T_iterator_const> equal_range(
      I start, I end) const noexcept;
  [[nodiscard]] inline std::pair<T_iterator_const, T_iterator_const> find_overlaps(
      I start, I end) const noexcept;
  [[nodiscard]] inline T_iterator_const lower_bound(I start, I end) const noexcept;
  [[nodiscard]] inline T_iterator_const upper_bound(I start, I end) const noexcept;

  [[nodiscard]] inline std::pair<size_t, size_t> equal_range_idx(I start, I end) const noexcept;
  [[nodiscard]] inline std::pair<size_t, size_t> find_overlaps_idx(I start, I end) const noexcept;
  [[nodiscard]] inline size_t lower_bound_idx(I start, I end) const noexcept;
  [[nodiscard]] inline size_t upper_bound_idx(I start, I end) const noexcept;

  [[nodiscard]] inline size_t count(I start, I end) const noexcept;

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
  bool _indexed{true};

  std::vector<I> _start;
  std::vector<I> _end;
  std::vector<I> _max;
  std::vector<T> _data;

  // clang-format off
  enum QueryType : uint8_t {
    LOWER_BOUND = 0x01,
    UPPER_BOUND = 0x02,
    EQUAL_RANGE = 0x04,
    CONTAINS =    0x08
  };
  // clang-format on

  template <QueryType Mode>
  inline std::pair<T_iterator_const, T_iterator_const> internal_equal_range(I start,
                                                                            I end) const noexcept;
};
}  // namespace modle

#include "../../interval_tree_impl.hpp"  // IWYU pragma: export

// IWYU pragma: no_include "../../interval_tree_impl.hpp"
