// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/strings/match.h>  // for StartsWithIgnoreCase
#include <absl/types/span.h>     // for MakeSpan, Span
#include <xxhash.h>              // for XXH_INLINE_XXH3_freeState, XXH3_f...

#include <cassert>           // for assert
#include <exception>         // for exception
#include <functional>        // for reference_wrapper
#include <future>            // for future, promise
#include <initializer_list>  // for initializer_list
#include <string_view>       // for string_view, basic_string_view
#include <type_traits>
#include <utility>  // for pair, make_pair, forward
#include <vector>   // for vector

#include "modle/common/common.hpp"                      // for usize, i64, u64
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...

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
  return this->_value;
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
  this->_mutexes = mutexes;
}

template <class MutexT>
LockRangeExclusive<MutexT>::LockRangeExclusive(std::vector<MutexT> &mutexes)
    : LockRangeExclusive(absl::MakeSpan(mutexes)) {}

template <class MutexT>
LockRangeExclusive<MutexT>::~LockRangeExclusive() noexcept {
  std::for_each(this->_mutexes.begin(), this->_mutexes.end(), [](auto &m) { m.unlock(); });
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
  this->_mutexes = mutexes;
}

template <class MutexT>
LockRangeShared<MutexT>::LockRangeShared(std::vector<MutexT> &mutexes)
    : LockRangeShared(absl::MakeSpan(mutexes)) {}

template <class MutexT>
LockRangeShared<MutexT>::~LockRangeShared() noexcept {
  std::for_each(this->_mutexes.begin(), this->_mutexes.end(), [](auto &m) { m.unlock_shared(); });
}

constexpr std::string_view strip_quotes(std::string_view s) noexcept {
  assert(!s.empty());
  if (s.front() == '\'' || s.front() == '"') {
    s = s.substr(1);
  }
  assert(!s.empty());
  if (s.back() == '\'' || s.back() == '"') {
    s = s.substr(0, s.size() - 1);
  }
  assert(!s.empty());
  return s;
}
}  // namespace modle::utils

// IWYU pragma: private, include "modle/utils.hpp"
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <cstdlib>
// IWYU pragma: no_include <memory>
