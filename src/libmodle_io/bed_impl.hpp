// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>
#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <filesystem>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"

namespace modle::bed {

constexpr BED::operator bool() const noexcept { return id() != null_id; }

constexpr BED::Dialect BED::get_standard() const noexcept { return _standard; }

constexpr std::size_t BED::id() const noexcept { return _id; }

constexpr std::size_t BED::size() const noexcept { return chrom_end - chrom_start; }

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
  [[maybe_unused]] auto [node, inserted] = _trees.try_emplace(interval.chrom, IITree_t{});
  node->second.insert(I(interval.chrom_start), I(interval.chrom_end), std::move(interval));
  return std::make_pair(node, true);
}

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::insert(const K& chrom_name, value_type tree) {
  return _trees.emplace(chrom_name, std::move(tree));
}

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::insert(const K& chrom_name, I chrom_start,
                                                         I chrom_end) {
  auto node = _trees.try_emplace(std::move(chrom_name), IITree_t{}).first;
  node->second.insert(chrom_start, chrom_end, BED{chrom_name, bp_t(chrom_start), bp_t(chrom_end)});
  return std::make_pair(node, true);
}

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::emplace(BED&& interval) {
#if __GNUC__ > 7
  [[maybe_unused]] auto [node, inserted] = _trees.try_emplace(interval.chrom, IITree_t{});
  node->second.emplace(I(interval.chrom_start), I(interval.chrom_end), std::move(interval));
  return std::make_pair(node, true);
#else
  // This inefficient, ugly piece of code is required to workaround a mysterious
  // -Wduplicated-branches warning that is raised on GCC 7.5
  if (_trees.contains(interval.chrom)) {
    _trees.at(interval.chrom).insert(I(interval.chrom_start), I(interval.chrom_end), interval);
    return std::make_pair(_trees.find(interval.chrom), true);
  }
  _trees.emplace(interval.chrom, IITree_t{});
  _trees.at(interval.chrom).insert(I(interval.chrom_start), I(interval.chrom_end), interval);
  return std::make_pair(_trees.find(interval.chrom), true);
#endif
}

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::emplace(const K& chrom_name, value_type&& tree) {
  return _trees.emplace(chrom_name, tree);
}

template <typename K, typename I>
void BED_tree<K, I>::insert(const absl::Span<const BED> intervals) {
  for (const auto& interval : intervals) {
    [[maybe_unused]] auto [node, inserted] = _trees.try_emplace(interval.chrom, IITree_t{});
    node->second.insert(I(interval.chrom_start), I(interval.chrom_end), interval);
  }
}

template <typename K, typename I>
void BED_tree<K, I>::emplace(std::vector<BED>&& intervals) {
  for (auto&& interval : intervals) {
    [[maybe_unused]] auto [node, inserted] = _trees.try_emplace(interval.chrom, IITree_t{});
    node->second.emplace(I(interval.chrom_start), I(interval.chrom_end), std::move(interval));
  }
}

template <typename K, typename I>
const typename BED_tree<K, I>::value_type& BED_tree<K, I>::at(const K& chrom) const {
  return _trees.at(chrom);
}

template <typename K, typename I>
void BED_tree<K, I>::index() {
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_DUPLICATED_BRANCHES
  std::for_each(_trees.begin(), _trees.end(), [](auto& node) { node.second.make_BST(); });
  DISABLE_WARNING_POP
}

template <typename K, typename I>
void BED_tree<K, I>::index(const K& chrom_name) {
  _trees.at(chrom_name).make_BST();
}

template <typename K, typename I>
bool BED_tree<K, I>::contains(const K& chrom_name) const {
  return _trees.contains(chrom_name);
}

template <typename K, typename I>
bool BED_tree<K, I>::contains_overlap(const BED& interval) const {
  return contains_overlap(interval.chrom, I(interval.chrom_start), I(interval.chrom_end));
}

template <typename K, typename I>
bool BED_tree<K, I>::contains_overlap(const K& chrom_name, std::uint64_t chrom_start,
                                      std::uint64_t chrom_end) const {
  auto it = _trees.find(chrom_name);
  if (it == _trees.end()) {
    return false;
  }
  assert(it->second.is_BST());

  return it->second.overlaps_with(chrom_start, chrom_end);
}

template <typename K, typename I>
std::size_t BED_tree<K, I>::count_overlaps(const BED& interval) const {
  return count_overlaps(interval.chrom, I(interval.chrom_start), I(interval.chrom_end));
}

template <typename K, typename I>
std::size_t BED_tree<K, I>::count_overlaps(const K& chrom_name, std::uint64_t chrom_start,
                                           std::uint64_t chrom_end) const {
  auto it = _trees.find(chrom_name);
  if (it == _trees.end()) {
    return 0UL;
  }
  assert(it->second.is_BST());

  return it->second.count(chrom_start, chrom_end);
}

template <typename K, typename I>
absl::Span<const BED> BED_tree<K, I>::find_overlaps(const K& chrom_name, std::uint64_t chrom_start,
                                                    std::uint64_t chrom_end) const {
  auto it = _trees.find(chrom_name);
  if (it == _trees.end()) {
    return absl::Span<const BED>{};
  }
  assert(it->second.is_BST());

  const auto [overlap_begin, overlap_end] = it->second.find_overlaps(chrom_start, chrom_end);
  if (overlap_begin == it->second.data_end()) {
    return absl::Span<const BED>{};
  }
  return absl::MakeConstSpan(&(*overlap_begin),
                             static_cast<std::size_t>(overlap_end - overlap_begin));
}

template <typename K, typename I>
absl::Span<const BED> BED_tree<K, I>::find_overlaps(const BED& interval) const {
  return find_overlaps(interval.chrom, interval.chrom_start, interval.chrom_end);
}

template <typename K, typename I>
bool BED_tree<K, I>::empty() const {
  return _trees.empty();
}

template <typename K, typename I>
std::size_t BED_tree<K, I>::size() const {
  return std::accumulate(_trees.begin(), _trees.end(), std::size_t(0),
                         [](const auto accumulator, const auto& node) {
                           const auto& tree = node.second;
                           return accumulator + tree.size();
                         });
}

template <typename K, typename I>
std::size_t BED_tree<K, I>::size(const K& chrom_name) const {
  return at(chrom_name).size();
}

template <typename K, typename I>
void BED_tree<K, I>::clear() {
  _trees.clear();
}

template <typename K, typename I>
void BED_tree<K, I>::clear(const K& chrom_name) {
  _trees.at(chrom_name).clear();
}

template <typename K, typename I>
typename BED_tree<K, I>::iterator BED_tree<K, I>::begin() {
  return _trees.begin();
}

template <typename K, typename I>
typename BED_tree<K, I>::iterator BED_tree<K, I>::end() {
  return _trees.end();
}

template <typename K, typename I>
typename BED_tree<K, I>::const_iterator BED_tree<K, I>::begin() const {
  return _trees.begin();
}

template <typename K, typename I>
typename BED_tree<K, I>::const_iterator BED_tree<K, I>::end() const {
  return _trees.end();
}

template <typename K, typename I>
typename BED_tree<K, I>::const_iterator BED_tree<K, I>::cbegin() const {
  return _trees.cbegin();
}

template <typename K, typename I>
typename BED_tree<K, I>::const_iterator BED_tree<K, I>::cend() const {
  return _trees.cend();
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
auto fmt::formatter<modle::bed::RGB>::format(const modle::bed::RGB& rgb, FormatContext& ctx) const
    -> decltype(ctx.out()) {
  return fmt::format_to(ctx.out(), "{},{},{}", rgb.r, rgb.g, rgb.b);
}

constexpr auto fmt::formatter<modle::bed::BED>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  const auto* it = ctx.begin();
  const auto* end = ctx.end();
  const auto fmt_string =
      std::string_view{&(*ctx.begin()), static_cast<std::size_t>(ctx.end() - ctx.begin())};
  if (it != end) {
    for (const auto& [k, v] : presentation_mappings) {
      if (const auto pos = fmt_string.find(k); pos != std::string_view::npos) {
        presentation = v;
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
auto fmt::formatter<modle::bed::BED>::format(const modle::bed::BED& b, FormatContext& ctx) const
    -> decltype(ctx.out()) {
  assert(!b.chrom.empty());
  const auto n = std::min(b.num_fields(), static_cast<std::size_t>(presentation));
  auto it = ctx.out();
  if (n == 3U) {
    it = fmt::format_to(it, "{}\t{}\t{}", b.chrom, b.chrom_start, b.chrom_end);
  } else if (n == 4U) {
    it = fmt::format_to(it, "{}\t{}\t{}\t{}", b.chrom, b.chrom_start, b.chrom_end, b.name);
  } else if (n == 5U) {
    it = fmt::format_to(it, "{}\t{}\t{}\t{}\t{}", b.chrom, b.chrom_start, b.chrom_end, b.name,
                        b.score);
  } else if (n == 6U) {
    it = fmt::format_to(it, "{}\t{}\t{}\t{}\t{}\t{}", b.chrom, b.chrom_start, b.chrom_end, b.name,
                        b.score, b.strand);
  } else if (n == 9U) {
    assert(b.rgb);
    it = fmt::format_to(it, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", b.chrom, b.chrom_start,
                        b.chrom_end, b.name, b.score, b.strand, b.thick_start, b.thick_end, *b.rgb);
  } else if (n >= 12U) {
    assert(b.rgb);
    it = fmt::format_to(it, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", b.chrom, b.chrom_start,
                        b.chrom_end, b.name, b.score, b.strand, b.thick_start, b.thick_end, *b.rgb,
                        fmt::join(b.block_sizes, ","), fmt::join(b.block_starts, ","));
  }

  if (!b.extra_tokens.empty()) {
    it = fmt::format_to(it, "\t{}", b.extra_tokens);
  }

  return it;
}
