// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span
#include <xxhash.h>           // for XXH3_state_t, XXH_INLINE_XXH3_state_t

#include <array>                      // for array
#include <boost/filesystem/path.hpp>  // for path
#include <functional>                 // for reference_wrapper
#include <future>                     // for future
#include <initializer_list>           // for initializer_list
#include <string>                     // for string
#include <string_view>                // for string_view
#include <type_traits>                // for enable_if, enable_if_t
#include <utility>                    // for pair
#include <vector>                     // for vector

#include "modle/common/common.hpp"  // for i64, usize, isize

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

[[maybe_unused]] [[nodiscard]] constexpr bool ndebug_defined() noexcept;
[[maybe_unused]] [[nodiscard]] constexpr bool ndebug_not_defined() noexcept;

#ifdef __has_cpp_attribute
#if __has_cpp_attribute(clang::no_destroy)
#define MODLE_NO_DESTROY [[clang::no_destroy]]
#else
#define MODLE_NO_DESTROY
#endif
#else
#define MODLE_NO_DESTROY
#endif

// Source: https://www.youtube.com/watch?v=INn3xa4pMfg
template <class Key, class Value, usize Size>
class ConstMap {
  std::array<std::pair<Key, Value>, Size> _buff;

 public:
  template <class... Args>
  constexpr ConstMap(Args&&... args) noexcept;

  using const_iterator = typename std::array<std::pair<Key, Value>, Size>::const_iterator;
  using key_type = Key;
  using mapped_type = Value;
  using value_type = std::pair<const Key, Value>;
  using first_type = const Key;
  using second_type = Value;

  [[nodiscard]] constexpr const Value& at(const Key& key) const;
  [[nodiscard]] constexpr const Value& operator[](const Key& key) const;
  [[nodiscard]] constexpr const_iterator find(const Key& key) const noexcept;
  [[nodiscard]] constexpr bool contains(const Key& key) const noexcept;

  [[nodiscard]] constexpr const_iterator begin() const noexcept;
  [[nodiscard]] constexpr const_iterator end() const noexcept;
};

template <class T>
class RepeatIterator {
  T _value;

 public:
  using value_type = T;
  using difference_type = isize;
  using pointer = T*;
  using reference = T&;
  using iterator_category = std::bidirectional_iterator_tag;

  RepeatIterator() = delete;
  explicit RepeatIterator(T value);

  [[nodiscard]] constexpr const T& operator*() const;
  [[nodiscard]] constexpr const T& operator[](usize i) const;

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

template <class InputIt1, class InputIt2, class N = double,
          class = std::enable_if_t<std::is_arithmetic_v<N>>>
[[nodiscard]] constexpr N convolve(InputIt1 kernel_first, InputIt1 kernel_last,
                                   InputIt2 buff_first);

template <class Rng1, class Rng2, class N = double,
          class = std::enable_if_t<std::is_arithmetic_v<N>>>
[[nodiscard]] constexpr N convolve(const Rng1& kernel, const Rng2& buff);

template <class I, class = std::enable_if_t<std::is_integral_v<I>>>
[[nodiscard]] constexpr I next_pow2(I n) noexcept;

template <class T>
[[nodiscard]] constexpr std::future<T> make_ready_future(T&& v);

template <bool remove_source_files = false>
inline void concatenate_files(
    const boost::filesystem::path& path_to_dest,
    std::initializer_list<std::reference_wrapper<const boost::filesystem::path>> path_to_sources);

template <bool remove_source_files = false>
inline void concatenate_files(const boost::filesystem::path& path_to_dest,
                              const boost::filesystem::path& path_to_src);

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

}  // namespace modle::utils

#include "../../../utils_impl.hpp"  // IWYU pragma: export
// IWYU pragma: no_include <iterator>
// IWYU pragma: no_include <boost/container/detail/std_fwd.hpp>
