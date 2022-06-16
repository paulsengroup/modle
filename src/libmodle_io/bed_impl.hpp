// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span, MakeConstSpan
#include <fmt/format.h>       // for format_parse_context, FMT_STRING, join, format_error

#include <algorithm>    // for min
#include <array>        // for array
#include <cassert>      // for assert
#include <filesystem>   // for path
#include <memory>       // for unique_ptr
#include <numeric>      // for accumulate
#include <string>       // for string
#include <string_view>  // for string_view, basic_string_view
#include <tuple>        // for ignore
#include <type_traits>  // for add_const<>::type
#include <utility>      // for move, pair, make_pair
#include <vector>       // for vector

#include "modle/common/common.hpp"  // for bp_t, u64, u32, u8, usize

namespace modle::bed {

template <typename K, typename I>
BED_tree<K, I>::BED_tree(const std::filesystem::path& path_to_bed, BED::Dialect dialect)
    : BED_tree(path_to_bed.empty()
                   ? BED_tree<>{}
                   : bed::Parser(path_to_bed, dialect).parse_all_in_interval_tree()) {}

// For some reason GCC 11 is confused by this alias...
// template <typename K, typename I>
// using value_type_t = typename BED_tree<K, I>::value_type;
template <typename K, typename I>
using iterator_t = typename BED_tree<K, I>::iterator;

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::insert(BED interval) {
  [[maybe_unused]] auto [node, inserted] = this->_trees.try_emplace(interval.chrom, IITree_t{});
  node->second.insert(I(interval.chrom_start), I(interval.chrom_end), std::move(interval));
  return std::make_pair(node, true);
}

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::insert(const K& chrom_name, value_type tree) {
  return this->_trees.emplace(chrom_name, std::move(tree));
}

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::insert(const K& chrom_name, I chrom_start,
                                                         I chrom_end) {
  auto node = this->_trees.try_emplace(std::move(chrom_name), IITree_t{}).first;
  node->second.insert(chrom_start, chrom_end, BED{chrom_name, bp_t(chrom_start), bp_t(chrom_end)});
  return std::make_pair(node, true);
}

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::emplace(BED&& interval) {
#if __GNUC__ > 7
  [[maybe_unused]] auto [node, inserted] = this->_trees.try_emplace(interval.chrom, IITree_t{});
  node->second.emplace(I(interval.chrom_start), I(interval.chrom_end), std::move(interval));
  return std::make_pair(node, true);
#else
  // This inefficient, ugly piece of code is required to workaround a mysterious
  // -Wduplicated-branches warning that is raised on GCC 7.5
  if (this->_trees.contains(interval.chrom)) {
    this->_trees.at(interval.chrom)
        .insert(I(interval.chrom_start), I(interval.chrom_end), interval);
    return std::make_pair(this->_trees.find(interval.chrom), true);
  }
  this->_trees.emplace(interval.chrom, IITree_t{});
  this->_trees.at(interval.chrom).insert(I(interval.chrom_start), I(interval.chrom_end), interval);
  return std::make_pair(this->_trees.find(interval.chrom), true);
#endif
}

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::emplace(const K& chrom_name, value_type&& tree) {
  return this->_trees.emplace(chrom_name, tree);
}

template <typename K, typename I>
void BED_tree<K, I>::insert(const absl::Span<const BED> intervals) {
  for (const auto& interval : intervals) {
    [[maybe_unused]] auto [node, inserted] = this->_trees.try_emplace(interval.chrom, IITree_t{});
    node->second.insert(I(interval.chrom_start), I(interval.chrom_end), interval);
  }
}

template <typename K, typename I>
void BED_tree<K, I>::emplace(std::vector<BED>&& intervals) {
  for (auto&& interval : intervals) {
    [[maybe_unused]] auto [node, inserted] = this->_trees.try_emplace(interval.chrom, IITree_t{});
    node->second.emplace(I(interval.chrom_start), I(interval.chrom_end), std::move(interval));
  }
}

template <typename K, typename I>
const typename BED_tree<K, I>::value_type& BED_tree<K, I>::at(const K& chrom) const {
  return this->_trees.at(chrom);
}

template <typename K, typename I>
void BED_tree<K, I>::index() {
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_DUPLICATED_BRANCHES
  std::for_each(this->_trees.begin(), this->_trees.end(),
                [](auto& node) { node.second.make_BST(); });
  DISABLE_WARNING_POP
}

template <typename K, typename I>
void BED_tree<K, I>::index(const K& chrom_name) {
  this->_trees.at(chrom_name).make_BST();
}

template <typename K, typename I>
bool BED_tree<K, I>::contains(const K& chrom_name) const {
  return this->_trees.contains(chrom_name);
}

template <typename K, typename I>
bool BED_tree<K, I>::contains_overlap(const BED& interval) const {
  return this->contains_overlap(interval.chrom, I(interval.chrom_start), I(interval.chrom_end));
}

template <typename K, typename I>
bool BED_tree<K, I>::contains_overlap(const K& chrom_name, u64 chrom_start, u64 chrom_end) const {
  auto it = this->_trees.find(chrom_name);
  if (it == this->_trees.end()) {
    return false;
  }
  assert(it->second.is_BST());

  return it->second.overlaps_with(chrom_start, chrom_end);
}

template <typename K, typename I>
usize BED_tree<K, I>::count_overlaps(const BED& interval) const {
  return this->count_overlaps(interval.chrom, I(interval.chrom_start), I(interval.chrom_end));
}

template <typename K, typename I>
usize BED_tree<K, I>::count_overlaps(const K& chrom_name, u64 chrom_start, u64 chrom_end) const {
  auto it = this->_trees.find(chrom_name);
  if (it == this->_trees.end()) {
    return 0UL;
  }
  assert(it->second.is_BST());

  return it->second.count(chrom_start, chrom_end);
}

template <typename K, typename I>
absl::Span<const BED> BED_tree<K, I>::find_overlaps(const K& chrom_name, u64 chrom_start,
                                                    u64 chrom_end) const {
  auto it = this->_trees.find(chrom_name);
  if (it == this->_trees.end()) {
    return absl::Span<const BED>{};
  }
  assert(it->second.is_BST());

  const auto [overlap_begin, overlap_end] = it->second.find_overlaps(chrom_start, chrom_end);
  if (overlap_begin == it->second.data_end()) {
    return absl::Span<const BED>{};
  }
  return absl::MakeConstSpan(&(*overlap_begin), static_cast<usize>(overlap_end - overlap_begin));
}

template <typename K, typename I>
absl::Span<const BED> BED_tree<K, I>::find_overlaps(const BED& interval) const {
  return this->find_overlaps(interval.chrom, interval.chrom_start, interval.chrom_end);
}

template <typename K, typename I>
bool BED_tree<K, I>::empty() const {
  return this->_trees.empty();
}

template <typename K, typename I>
usize BED_tree<K, I>::size() const {
  return std::accumulate(this->_trees.begin(), this->_trees.end(), usize(0),
                         [](const auto accumulator, const auto& node) {
                           const auto& tree = node.second;
                           return accumulator + tree.size();
                         });
}

template <typename K, typename I>
usize BED_tree<K, I>::size(const K& chrom_name) const {
  return this->at(chrom_name).size();
}

template <typename K, typename I>
void BED_tree<K, I>::clear() {
  this->_trees.clear();
}

template <typename K, typename I>
void BED_tree<K, I>::clear(const K& chrom_name) {
  this->_trees.at(chrom_name).clear();
}

template <typename K, typename I>
typename BED_tree<K, I>::iterator BED_tree<K, I>::begin() {
  return this->_trees.begin();
}

template <typename K, typename I>
typename BED_tree<K, I>::iterator BED_tree<K, I>::end() {
  return this->_trees.end();
}

template <typename K, typename I>
typename BED_tree<K, I>::const_iterator BED_tree<K, I>::begin() const {
  return this->_trees.begin();
}

template <typename K, typename I>
typename BED_tree<K, I>::const_iterator BED_tree<K, I>::end() const {
  return this->_trees.end();
}

template <typename K, typename I>
typename BED_tree<K, I>::const_iterator BED_tree<K, I>::cbegin() const {
  return this->_trees.cbegin();
}

template <typename K, typename I>
typename BED_tree<K, I>::const_iterator BED_tree<K, I>::cend() const {
  return this->_trees.cend();
}
}  // namespace modle::bed

constexpr auto fmt::formatter<modle::bed::RGB>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.end();
}

template <typename FormatContext>
auto fmt::formatter<modle::bed::RGB>::format(const modle::bed::RGB& rgb, FormatContext& ctx)
    -> decltype(ctx.out()) {
  return fmt::format_to(ctx.out(), FMT_STRING("{},{},{}"), rgb.r, rgb.g, rgb.b);
}

constexpr auto fmt::formatter<modle::bed::BED>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  const auto* it = ctx.begin();
  const auto* end = ctx.end();
  const auto fmt_string =
      std::string_view{&(*ctx.begin()), static_cast<modle::usize>(ctx.end() - ctx.begin())};
  if (it != end) {
    for (const auto& [k, v] : this->presentation_mappings) {
      if (const auto pos = fmt_string.find(k); pos != std::string_view::npos) {
        this->presentation = v;
        it += k.size();
        break;
      }
    }
  }
  // Check if reached the end of the range:
  if (it != end && *it != '}') {
    throw fmt::format_error("invalid format");
  }

  // Return an iterator past the end of the parsed range:
  return it;
}

// Formats the point p using the parsed format specification (presentation)
// stored in this formatter.
template <typename FormatContext>
auto fmt::formatter<modle::bed::BED>::format(const modle::bed::BED& b, FormatContext& ctx)
    -> decltype(ctx.out()) {
  assert(!b.chrom.empty());
  const auto n = std::min(b.num_fields(), static_cast<modle::usize>(this->presentation));
  auto it = ctx.out();
  if (n == 3U) {
    it = fmt::format_to(it, FMT_STRING("{}\t{}\t{}"), b.chrom, b.chrom_start, b.chrom_end);
  } else if (n == 4U) {
    it = fmt::format_to(it, FMT_STRING("{}\t{}\t{}\t{}"), b.chrom, b.chrom_start, b.chrom_end,
                        b.name);
  } else if (n == 5U) {
    it = fmt::format_to(it, FMT_STRING("{}\t{}\t{}\t{}\t{}"), b.chrom, b.chrom_start, b.chrom_end,
                        b.name, b.score);
  } else if (n == 6U) {
    it = fmt::format_to(it, FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}"), b.chrom, b.chrom_start,
                        b.chrom_end, b.name, b.score, b.strand);
  } else if (n == 9U) {
    assert(b.rgb);
    it =
        fmt::format_to(it, FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"), b.chrom, b.chrom_start,
                       b.chrom_end, b.name, b.score, b.strand, b.thick_start, b.thick_end, *b.rgb);
  } else if (n >= 12U) {
    assert(b.rgb);
    it = fmt::format_to(it, FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"), b.chrom,
                        b.chrom_start, b.chrom_end, b.name, b.score, b.strand, b.thick_start,
                        b.thick_end, *b.rgb, fmt::join(b.block_sizes, ","),
                        fmt::join(b.block_starts, ","));
  }

  if (!b.extra_tokens.empty()) {
    it = fmt::format_to(it, FMT_STRING("\t{}"), b.extra_tokens);
  }

  return it;
}
