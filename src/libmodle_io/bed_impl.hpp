#pragma once

#include <absl/container/btree_map.h>  // for btree_map

#include <cassert>  // for assert
#include <cstddef>  // for size_t

#include "modle/interval_tree.hpp"  // for IITree

namespace modle::bed {

template <typename K, typename I>
BED_tree<K, I>::BED_tree(const boost::filesystem::path& path_to_bed, BED::Dialect dialect)
    : BED_tree(bed::Parser(path_to_bed, dialect).parse_all_in_interval_tree()) {}

// For some reason GCC 11 is confused by this alias...
// template <typename K, typename I>
// using value_type_t = typename BED_tree<K, I>::value_type;
template <typename K, typename I>
using iterator_t = typename BED_tree<K, I>::iterator;

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::insert(BED interval) {
  auto [node, inserted] = this->_trees.try_emplace(interval.chrom, IITree_t{});
  (void)inserted;
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
  auto [node, inserted] = this->_trees.try_emplace(chrom_name, IITree_t{});
  (void)inserted;
  node->second.insert(chrom_start, chrom_end, BED{chrom_name, bp_t(chrom_start), bp_t(chrom_end)});
  return std::make_pair(node, true);
}

template <typename K, typename I>
std::pair<iterator_t<K, I>, bool> BED_tree<K, I>::emplace(BED&& interval) {
#if __GNUC__ > 7
  auto [node, inserted] = this->_trees.try_emplace(interval.chrom, IITree_t{});
  (void)inserted;
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
  return this->_trees.emplace(chrom_name, std::move(tree));
}

template <typename K, typename I>
void BED_tree<K, I>::insert(const absl::Span<const BED> intervals) {
  for (const auto& interval : intervals) {
    auto [node, inserted] = this->_trees.try_emplace(interval.chrom, IITree_t{});
    (void)inserted;
    node->second.insert(I(interval.chrom_start), I(interval.chrom_end), interval);
  }
}

template <typename K, typename I>
void BED_tree<K, I>::emplace(std::vector<BED>&& intervals) {
  for (auto&& interval : intervals) {
    auto [node, inserted] = this->_trees.try_emplace(interval.chrom, IITree_t{});
    (void)inserted;
    node->second.emplace(I(interval.chrom_start), I(interval.chrom_end), std::move(interval));
  }
}

template <typename K, typename I>
const typename BED_tree<K, I>::value_type& BED_tree<K, I>::at(const K& chrom) const {
  return this->_trees.at(chrom);
}

template <typename K, typename I>
void BED_tree<K, I>::index() {
  // TODO Workaround the -Wduplicated-branches produced by GCC 7.5
  std::for_each(this->_trees.begin(), this->_trees.end(),
                [](auto& node) { node.second.make_BST(); });
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
bool BED_tree<K, I>::contains_overlap(const K& chrom_name, uint64_t chrom_start,
                                      uint64_t chrom_end) const {
  auto it = this->_trees.find(chrom_name);
  if (it == this->_trees.end()) {
    return false;
  }
  assert(it->second.is_BST());  // NOLINT You forgot to call index/make_BST()!

  return it->second.overlaps_with(chrom_start, chrom_end);
}

template <typename K, typename I>
size_t BED_tree<K, I>::count_overlaps(const BED& interval) const {
  return this->count_overlaps(interval.chrom, I(interval.chrom_start), I(interval.chrom_end));
}

template <typename K, typename I>
size_t BED_tree<K, I>::count_overlaps(const K& chrom_name, uint64_t chrom_start,
                                      uint64_t chrom_end) const {
  auto it = this->_trees.find(chrom_name);
  if (it == this->_trees.end()) {
    return 0UL;
  }
  assert(it->second.is_BST());  // NOLINT You forgot to call index/make_BST()!

  return it->second.count(chrom_start, chrom_end);
}

template <typename K, typename I>
absl::Span<const BED> BED_tree<K, I>::find_overlaps(const K& chrom_name, uint64_t chrom_start,
                                                    uint64_t chrom_end) const {
  auto it = this->_trees.find(chrom_name);
  if (it == this->_trees.end()) {
    return absl::Span<const BED>{};
  }
  assert(it->second.is_BST());  // NOLINT You forgot to call index/make_BST()!

  const auto [overlap_begin, overlap_end] = it->second.find_overlaps(chrom_start, chrom_end);
  if (overlap_begin == it->second.data_end()) {
    return absl::Span<const BED>{};
  }
  return absl::MakeConstSpan(&(*overlap_begin), static_cast<size_t>(overlap_end - overlap_begin));
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
size_t BED_tree<K, I>::size() const {
  return std::accumulate(this->_trees.begin(), this->_trees.end(), 0UL,
                         [](const auto accumulator, const auto& node) {
                           const auto& tree = node.second;
                           return accumulator + tree.size();
                         });
}

template <typename K, typename I>
size_t BED_tree<K, I>::size(const K& chrom_name) const {
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
