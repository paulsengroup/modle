// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <initializer_list>
#include <stdexcept>
#include <utility>

#include "modle/common/common.hpp"

namespace modle::utils {
template <class Key, class Value, usize Size>
constexpr ConstMap<Key, Value, Size>::ConstMap(
    const std::initializer_list<std::pair<Key, Value>> &args) {
  if (args.size() < Size) {
    throw std::logic_error("Not enough arguments to initialize ConstMap");
  }
  if (args.size() > Size) {
    throw std::out_of_range("Too many arguments to initialize ConstMap");
  }

  usize i = 0;
  for (const auto &[k, v] : args) {
    _keys[i] = k;
    _vals[i] = v;
    ++i;
  }

  detect_duplicate_keys();
}

template <class Key, class Value, usize Size>
template <class KeyInputIt, class ValueInputIt>
constexpr ConstMap<Key, Value, Size>::ConstMap(KeyInputIt first_key, KeyInputIt last_key,
                                               ValueInputIt first_value) {
  const auto dist = static_cast<usize>(std::distance(first_key, last_key));
  if (dist < Size) {
    throw std::logic_error("Not enough values to initialize ConstMap");
  }
  if (dist > Size) {
    throw std::out_of_range("Too many values to initialize ConstMap");
  }

  std::move(first_key, last_key, _keys.begin());
  std::move(first_value, first_value + static_cast<isize>(dist), _vals.begin());

  detect_duplicate_keys();
}

template <class Key, class Value, usize Size>
constexpr void ConstMap<Key, Value, Size>::detect_duplicate_keys() const {
  if constexpr (utils::ndebug_not_defined()) {
    for (const auto &k : _keys) {
      usize i = 0;
      for (const auto &kk : _keys) {
        i += k == kk;
      }
      if (i != 1) {
        throw std::logic_error("Detected a duplicate key");
      }
    }
  }
}

template <class Key, class Value, usize Size>
constexpr bool ConstMap<Key, Value, Size>::empty() const noexcept {
  return _keys.empty();
}

template <class Key, class Value, usize Size>
constexpr usize ConstMap<Key, Value, Size>::size() const noexcept {
  return _keys.size();
}

template <class Key, class Value, usize Size>
constexpr const Value &ConstMap<Key, Value, Size>::at(const Key &key) const {
  const auto itr = find(key);
  if (itr != end()) {
    assert(itr.second);
    return *itr.second;
  }
  throw std::range_error("unable to find key");
}

template <class Key, class Value, usize Size>
constexpr const Value &ConstMap<Key, Value, Size>::operator[](const Key &key) const {
  return *find(key)->second;
}

template <class Key, class Value, usize Size>
constexpr auto ConstMap<Key, Value, Size>::find(const Key &key) const noexcept -> const_iterator {
  auto it = std::find(_keys.begin(), _keys.end(), key);
  return iterator(*this, static_cast<usize>(std::distance(_keys.begin(), it)));
}

template <class Key, class Value, usize Size>
constexpr bool ConstMap<Key, Value, Size>::contains(const Key &key) const noexcept {
  return find(key) != end();
}

template <class Key, class Value, usize Size>
constexpr auto ConstMap<Key, Value, Size>::begin() const noexcept -> const_iterator {
  return iterator(*this, 0);
}

template <class Key, class Value, usize Size>
constexpr auto ConstMap<Key, Value, Size>::end() const noexcept -> const_iterator {
  return iterator(*this, size());
}

template <class Key, class Value, usize Size>
constexpr const std::array<Key, Size> &ConstMap<Key, Value, Size>::keys() const noexcept {
  return _keys;
}

template <class Key, class Value, usize Size>
constexpr const std::array<Value, Size> &ConstMap<Key, Value, Size>::values() const noexcept {
  return _vals;
}

template <class Key, class Value, usize Size>
constexpr ConstMap<Key, Value, Size>::iterator::iterator(const ConstMap<Key, Value, Size> &map,
                                                         usize i) noexcept
    : first(map._keys.data() + i), second(map._vals.data() + i) {}

template <class Key, class Value, usize Size>
constexpr bool ConstMap<Key, Value, Size>::iterator::operator==(
    const ConstMap<Key, Value, Size>::iterator &other) const noexcept {
  return first == other.first;
}

template <class Key, class Value, usize Size>
constexpr bool ConstMap<Key, Value, Size>::iterator::operator!=(
    const ConstMap<Key, Value, Size>::iterator &other) const noexcept {
  return !(*this == other);
}

template <class Key, class Value, usize Size>
constexpr auto ConstMap<Key, Value, Size>::iterator::operator*() const -> const value_type {
  assert(first);
  assert(second);
  return std::make_pair(first, second);
}

template <class Key, class Value, usize Size>
constexpr auto ConstMap<Key, Value, Size>::iterator::operator++() noexcept -> iterator & {
  assert(first);
  assert(second);
  ++first;
  ++second;
  return *this;
}

template <class Key, class Value, usize Size>
constexpr auto ConstMap<Key, Value, Size>::iterator::operator++(int) const noexcept -> iterator {
  auto it = *this;
  std::ignore = ++it;
  return it;
}
}  // namespace modle::utils
