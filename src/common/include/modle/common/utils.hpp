// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <xxh3.h>  // for XXH3_state_t, XXH_INLINE_XXH3_state_t

#include <array>                             // for array
#include <boost/filesystem/file_status.hpp>  // for regular_file, file_type
#include <boost/filesystem/path.hpp>         // for path
#include <cstdio>                            // for FILE
#include <future>                            // for future, get_future, promise, set_value
#include <string>                            // for string
#include <string_view>                       // for string_view
#include <system_error>                      // for errc
#include <type_traits>                       // for declval
#include <utility>                           // for pair
#include <vector>                            // for vector

#include "modle/common/common.hpp"  // for usize, i64, u64

namespace modle::utils {

template <class N>
inline auto from_chars(const char* first, const char* last, N& value) noexcept;

template <class N>
inline void parse_numeric_or_throw(std::string_view tok, N& field);

template <class N>
inline N parse_numeric_or_throw(std::string_view tok);

template <class N>
inline void parse_numeric_or_throw(const std::vector<std::string_view>& toks, usize idx, N& field);

template <class N>
inline void parse_vect_of_numbers_or_throw(const std::vector<std::string_view>& toks, usize idx,
                                           std::vector<N>& fields, u64 expected_size);

template <class N>
inline void throw_except_from_errc(std::string_view tok, usize idx, const N& field, const char* c,
                                   std::errc e);

[[nodiscard]] inline bool chrom_equal_operator(std::string_view chr1, std::string_view chr2);
[[nodiscard]] inline bool chrom_equal_operator(const std::pair<std::string_view, i64>& chr1,
                                               const std::pair<std::string_view, i64>& chr2);

[[nodiscard]] inline bool chrom_less_than_operator(std::string_view chr1, std::string_view chr2);
[[nodiscard]] inline bool chrom_less_than_operator(const std::pair<std::string_view, i64>& chr1,
                                                   const std::pair<std::string_view, i64>& chr2);

// Typetraits stuff

// Usage: static_assert(std::experimental::is_detected_convertible_v<const char*,
// utils::has_data_member_func, T>, "error msg");
// This static_assert will fail for types that don't have a member function .data() that returns a
// const char*

template <class T>
using has_data_member_func = decltype(std::declval<T&>().data());

template <class T>
using has_size_member_func = decltype(std::declval<T&>().size());

template <class T>
[[nodiscard]] inline std::string get_printable_type_name(const T& var);

// https://stackoverflow.com/a/56766138
template <class T>
constexpr auto get_printable_type_name() noexcept;

// Various

struct XXH3_Deleter {  // NOLINT
  inline void operator()(XXH3_state_t* state) noexcept;
};

inline void fclose(FILE* fp) noexcept(false);

// Try to convert str representations like "1.0" or "1.000000" to "1"
[[nodiscard]] inline std::string str_float_to_str_int(const std::string& s);

// Returns false in case a collision was detected
inline bool detect_path_collision(
    const boost::filesystem::path& p, std::string& error_msg, bool force_overwrite = false,
    boost::filesystem::file_type expected_type = boost::filesystem::regular_file);
[[nodiscard]] inline std::string detect_path_collision(
    const boost::filesystem::path& p, bool force_overwrite = false,
    boost::filesystem::file_type expected_type = boost::filesystem::regular_file);

struct identity {
  template <class T>
  [[nodiscard]] constexpr T&& operator()(T&& a) const noexcept;
  using is_transparent = void;
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
}  // namespace modle::utils

template <bool remove_source_files = false>
inline void concatenate_files(
    const boost::filesystem::path& path_to_dest,
    std::initializer_list<std::reference_wrapper<const boost::filesystem::path>> path_to_sources);

template <bool remove_source_files = false>
inline void concatenate_files(const boost::filesystem::path& path_to_dest,
                              const boost::filesystem::path& path_to_sources);

#include "../../../utils_impl.hpp"  // IWYU pragma: export
// IWYU pragma: no_include <iterator>
// IWYU pragma: no_include <boost/container/detail/std_fwd.hpp>
