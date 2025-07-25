// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
// Copyright (C) 2019     Dana-Farber Cancer Institute
// Copyright (C) 2021     modified by Roberto Rossini
//
// SPDX-License-Identifier: MIT

#pragma once

#include <iterator>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"

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
template <class N, class T>
class IITree {
  DISABLE_WARNING_POP
  static_assert(std::is_arithmetic_v<N>);
  struct StackCell;

  using I_iterator = typename std::vector<N>::iterator;
  using T_iterator = typename std::vector<T>::iterator;

  using I_iterator_const = typename std::vector<N>::const_iterator;
  using T_iterator_const = typename std::vector<T>::const_iterator;

 public:  // NOLINTNEXTLINE(google-explicit-constructor)
  inline IITree(N start_pos_ = std::numeric_limits<N>::lowest(),
                N end_pos_ = (std::numeric_limits<N>::max)());

  inline void insert(N start, N end, const T &data);
  inline void emplace(N start, N end, T &&data);

  [[nodiscard]] inline bool overlaps_with(N start, N end) const noexcept;
  [[nodiscard]] inline std::pair<T_iterator_const, T_iterator_const> equal_range(
      N start, N end) const noexcept;
  [[nodiscard]] inline std::pair<T_iterator_const, T_iterator_const> find_overlaps(
      N start, N end) const noexcept;
  [[nodiscard]] inline T_iterator_const lower_bound(N start, N end) const noexcept;
  [[nodiscard]] inline T_iterator_const upper_bound(N start, N end) const noexcept;

  [[nodiscard]] inline std::pair<std::size_t, std::size_t> equal_range_idx(N start,
                                                                           N end) const noexcept;
  [[nodiscard]] inline std::pair<std::size_t, std::size_t> find_overlaps_idx(N start,
                                                                             N end) const noexcept;
  [[nodiscard]] inline std::size_t lower_bound_idx(N start, N end) const noexcept;
  [[nodiscard]] inline std::size_t upper_bound_idx(N start, N end) const noexcept;

  [[nodiscard]] inline std::size_t count(N start, N end) const noexcept;

  [[nodiscard]] inline constexpr std::size_t capacity() const noexcept;
  [[nodiscard]] inline constexpr bool empty() const noexcept;
  [[nodiscard]] inline constexpr std::size_t size() const noexcept;
  [[nodiscard]] inline constexpr bool is_BST() const noexcept;

  [[nodiscard]] inline constexpr std::size_t span() const noexcept;
  [[nodiscard]] inline constexpr N start_pos() const noexcept;
  [[nodiscard]] inline constexpr N end_pos() const noexcept;

  [[nodiscard]] inline N get_overlap_start(std::size_t i) const;
  [[nodiscard]] inline N get_overlap_end(std::size_t i) const;
  [[nodiscard]] inline const T &get_overlap_data(std::size_t i) const;

  [[nodiscard]] inline const std::span<const N> starts() const;
  [[nodiscard]] inline const std::span<const N> ends() const;
  [[nodiscard]] inline const std::span<const T> data() const;

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

  inline void reserve(std::size_t new_capacity);
  inline void make_BST();

 private:
  struct StackCell {
    inline StackCell() = default;
    inline StackCell(std::int64_t _level, std::size_t _node_idx, bool _left_child_alread_processed);

    std::int64_t level;
    std::size_t node_idx;
    bool left_child_already_processed;
  };

  N _start_pos;
  N _end_pos;
  std::int64_t _max_level{0};
  bool _indexed{true};

  std::vector<N> _start;
  std::vector<N> _end;
  std::vector<N> _max;
  std::vector<T> _data;

  // clang-format off
  enum QueryType : std::uint8_t {
    LOWER_BOUND = 0x01,
    UPPER_BOUND = 0x02,
    EQUAL_RANGE = 0x04,
    CONTAINS =    0x08
  };
  // clang-format on

  template <QueryType Mode>
  inline std::pair<T_iterator_const, T_iterator_const> internal_equal_range(N start,
                                                                            N end) const noexcept;
};
}  // namespace modle

#include "../../interval_tree_impl.hpp"  // IWYU pragma: export

// IWYU pragma: no_include "../../interval_tree_impl.hpp"
