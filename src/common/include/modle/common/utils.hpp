// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>
#include <xxhash.h>

#include <array>
#include <filesystem>
#include <functional>
#include <future>
#include <initializer_list>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"

namespace modle::utils {

struct identity {
  template <class T>
  [[nodiscard]] constexpr T&& operator()(T&& a) const noexcept;
  using is_transparent = void;
};

// Typetraits stuff
template <class T>
[[nodiscard]] inline std::string get_printable_type_name(const T& var);

// https://stackoverflow.com/a/56766138
template <class T>
constexpr auto get_printable_type_name() noexcept;

// Various

struct XXH3_Deleter {
  inline void operator()(XXH3_state_t* state) noexcept;
};

#ifdef __has_cpp_attribute
#if __has_cpp_attribute(clang::no_destroy)
#define MODLE_NO_DESTROY [[clang::no_destroy]]
#else
#define MODLE_NO_DESTROY
#endif
#else
#define MODLE_NO_DESTROY
#endif

template <class T>
class RepeatIterator {
  T _value;

 public:
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using pointer = T*;
  using reference = T&;
  using iterator_category = std::bidirectional_iterator_tag;

  RepeatIterator() = delete;
  explicit RepeatIterator(T value);

  [[nodiscard]] constexpr const T& operator*() const;
  [[nodiscard]] constexpr const T& operator[](std::size_t i) const;

  constexpr const RepeatIterator& operator++() const;
  constexpr const RepeatIterator operator++(int) const;
  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  constexpr const RepeatIterator& operator+=(I i) const;
  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  constexpr const RepeatIterator& operator+(I i) const;
  constexpr const RepeatIterator& operator--() const;
  constexpr const RepeatIterator operator--(int) const;
  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  constexpr const RepeatIterator& operator-=(I i) const;
  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  constexpr const RepeatIterator& operator-(I i) const;

  constexpr bool operator==(const RepeatIterator& other) const;
  constexpr bool operator!=(const RepeatIterator& other) const;
};

template <class I, class = std::enable_if_t<std::is_integral_v<I>>>
[[nodiscard]] constexpr I next_pow2(I n) noexcept;

template <class T>
[[nodiscard]] constexpr std::future<T> make_ready_future(T&& v);

template <bool remove_source_files = false>
void concatenate_files(const std::filesystem::path& path_to_dest,
                       const std::filesystem::path& path_to_sources...);

template <class MutexT>
class LockRangeExclusive {
  absl::Span<MutexT> _mutexes{};

 public:
  LockRangeExclusive() = default;
  explicit LockRangeExclusive(absl::Span<MutexT> mutexes);
  explicit LockRangeExclusive(std::vector<MutexT>& mutexes);
  ~LockRangeExclusive() noexcept;

  LockRangeExclusive(const LockRangeExclusive& other) = delete;
  LockRangeExclusive(LockRangeExclusive&& other) noexcept = delete;
  LockRangeExclusive& operator=(const LockRangeExclusive& other) = delete;
  LockRangeExclusive& operator=(LockRangeExclusive&& other) noexcept = delete;
};

template <class MutexT>
class LockRangeShared {
  absl::Span<MutexT> _mutexes{};

 public:
  LockRangeShared() = default;
  explicit LockRangeShared(absl::Span<MutexT> mutexes);
  explicit LockRangeShared(std::vector<MutexT>& mutexes);
  ~LockRangeShared() noexcept;

  LockRangeShared(const LockRangeShared& other) = delete;
  LockRangeShared(LockRangeShared&& other) noexcept = delete;
  LockRangeShared& operator=(const LockRangeShared& other) = delete;
  LockRangeShared& operator=(LockRangeShared&& other) noexcept = delete;
};

[[nodiscard]] constexpr std::string_view strip_quote_pairs(std::string_view s) noexcept;

}  // namespace modle::utils

#include "../../../utils_impl.hpp"  // IWYU pragma: export
// IWYU pragma: no_include <iterator>
// IWYU pragma: no_include <boost/container/detail/std_fwd.hpp>
