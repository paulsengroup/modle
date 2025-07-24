// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>
#include <xxhash.h>

#include <cassert>
#include <exception>
#include <functional>
#include <future>
#include <initializer_list>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"

namespace modle::utils {

template <class T>
constexpr auto get_printable_type_name() noexcept {
  std::string_view name = "Error: unsupported compiler";
  std::string_view prefix;
  std::string_view suffix;
#ifdef __clang__
  name = __PRETTY_FUNCTION__;
  prefix = "auto get_printable_type_name() [T = ";
  suffix = "]";
#elif defined(__GNUC__)
  name = __PRETTY_FUNCTION__;
  prefix = "constexpr auto type_name() [with T = ";
  suffix = "]";
#elif defined(_MSC_VER)
  name = __FUNCSIG__;
  prefix = "auto __cdecl type_name<";
  suffix = ">(void) noexcept";
#endif
  name.remove_prefix(prefix.size());
  name.remove_suffix(suffix.size());
  return name;
}

void XXH3_Deleter::operator()(XXH3_state_t *state) noexcept { XXH3_freeState(state); }

template <class T>
constexpr T &&identity::operator()(T &&a) const noexcept {
  return std::forward<T>(a);
}

template <class T>
RepeatIterator<T>::RepeatIterator(T value) : _value(std::move(value)) {}

template <class T>
constexpr const T &RepeatIterator<T>::operator*() const {
  return _value;
}

template <class T>
constexpr const T &RepeatIterator<T>::operator[]([[maybe_unused]] usize i) const {
  return *this;
}

template <class T>
constexpr const RepeatIterator<T> &RepeatIterator<T>::operator++() const {
  return *this;
}

template <class T>
constexpr const RepeatIterator<T> RepeatIterator<T>::operator++(int) const {
  return *this;
}

template <class T>
template <class I, class>
constexpr const RepeatIterator<T> &RepeatIterator<T>::operator+=([[maybe_unused]] I i) const {
  return *this;
}

template <class T>
template <class I, class>
constexpr const RepeatIterator<T> &RepeatIterator<T>::operator+([[maybe_unused]] I i) const {
  return *this;
}

template <class T>
constexpr const RepeatIterator<T> &RepeatIterator<T>::operator--() const {
  return *this;
}

template <class T>
constexpr const RepeatIterator<T> RepeatIterator<T>::operator--(int) const {
  return *this;
}

template <class T>
template <class I, class>
constexpr const RepeatIterator<T> &RepeatIterator<T>::operator-=([[maybe_unused]] I i) const {
  return *this;
}

template <class T>
template <class I, class>
constexpr const RepeatIterator<T> &RepeatIterator<T>::operator-([[maybe_unused]] I i) const {
  return *this;
}

template <class T>
constexpr bool RepeatIterator<T>::operator==([[maybe_unused]] const RepeatIterator &other) const {
  return false;
}

template <class T>
constexpr bool RepeatIterator<T>::operator!=([[maybe_unused]] const RepeatIterator &other) const {
  return !(*this == other);
}

template <class I, class>
constexpr I next_pow2(const I n) noexcept {
  using ull = unsigned long long;
  if constexpr (std::is_signed_v<I>) {
    assert(n >= 0);
    return utils::conditional_static_cast<I>(next_pow2(static_cast<ull>(n)));
  } else {
    auto m = utils::conditional_static_cast<ull>(n);
#ifndef __GNUC__
    // https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
    --m;
    m |= m >> 1;
    m |= m >> 2;
    m |= m >> 4;
    m |= m >> 8;
    m |= m >> 16;
    m |= m >> 32;
    return utils::conditional_static_cast<I>(m + 1);
#else
    // https://jameshfisher.com/2018/03/30/round-up-power-2/
    // https://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html

    return utils::conditional_static_cast<I>(
        m <= 1 ? m : u64(1) << (u64(64) - u64(__builtin_clzll(m - 1))));
#endif
  }
}

template <class T>
constexpr std::future<T> make_ready_future(T &&v) {
  std::promise<T> p;
  p.set_value(std::forward<T>(v));
  return p.get_future();
}

template <class MutexT>
LockRangeExclusive<MutexT>::LockRangeExclusive(const absl::Span<MutexT> mutexes) {
  usize i = 0;
  try {
    for (; i < mutexes.size(); ++i) {
      mutexes[i].lock();
    }
  } catch ([[maybe_unused]] const std::exception &e) {
    for (i = i - 1; i != 0; --i) {
      mutexes[i].unlock();
    }
    throw;
  }
  _mutexes = mutexes;
}

template <class MutexT>
LockRangeExclusive<MutexT>::LockRangeExclusive(std::vector<MutexT> &mutexes)
    : LockRangeExclusive(absl::MakeSpan(mutexes)) {}

template <class MutexT>
LockRangeExclusive<MutexT>::~LockRangeExclusive() noexcept {
  std::for_each(_mutexes.begin(), _mutexes.end(), [](auto &m) { m.unlock(); });
}

template <class MutexT>
LockRangeShared<MutexT>::LockRangeShared(const absl::Span<MutexT> mutexes) {
  usize i = 0;
  try {
    for (; i < mutexes.size(); ++i) {
      mutexes[i].lock();
    }
  } catch ([[maybe_unused]] const std::exception &e) {
    for (i = i - 1; i != 0; --i) {
      mutexes[i].unlock();
    }
    throw;
  }
  _mutexes = mutexes;
}

template <class MutexT>
LockRangeShared<MutexT>::LockRangeShared(std::vector<MutexT> &mutexes)
    : LockRangeShared(absl::MakeSpan(mutexes)) {}

template <class MutexT>
LockRangeShared<MutexT>::~LockRangeShared() noexcept {
  std::for_each(_mutexes.begin(), _mutexes.end(), [](auto &m) { m.unlock_shared(); });
}

constexpr std::string_view strip_quote_pairs(std::string_view s) noexcept {
  if (s.size() < 2) {
    return s;
  }
  const auto str_begins_with_quote = s.front() == '\'' || s.front() == '"';
  const auto str_ends_with_quote = s.back() == '\'' || s.back() == '"';
  if (str_begins_with_quote && str_ends_with_quote) {
    return s.substr(1, s.size() - 2);
  }
  return s;
}
}  // namespace modle::utils

// IWYU pragma: private, include "modle/utils.hpp"
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <cstdlib>
// IWYU pragma: no_include <memory>
