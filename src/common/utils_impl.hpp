// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// clang-format off
#include "modle/common/suppress_compiler_warnings.hpp"
#if defined(MODLE_CHARCONV_FP_AVAILABLE) && defined(MODLE_CHARCONV_INT_AVAILABLE)
#include <charconv>                          // for from_chars
#else
DISABLE_WARNING_PUSH
DISABLE_WARNING_CONVERSION
DISABLE_WARNING_DOUBLE_PROMOTION
DISABLE_WARNING_SHORTEN_64_TO_32
DISABLE_WARNING_SIGN_CONVERSION
#include <msstl/charconv.hpp>

DISABLE_WARNING_POP
#endif
// clang-format on

#include <absl/strings/match.h>       // for StartsWithIgnoreCase
#include <absl/strings/str_format.h>  // for StrAppendFormat
#include <absl/strings/str_split.h>   // for StrSplit, Splitter
#include <fmt/format.h>               // for format, FMT_STRING, system_error
#include <xxh3.h>                     // for XXH_INLINE_XXH3_freeState, XXH3_freeState

#include <boost/filesystem.hpp>              // for path, status
#include <boost/filesystem/file_status.hpp>  // for file_type, regular_file, directory_file, fil...
#include <boost/filesystem/operations.hpp>   // for status
#include <boost/filesystem/path.hpp>         // for path
#include <cassert>                           // for assert
#include <cerrno>                            // for errno
#include <cmath>                             // for trunc
#include <cstddef>                           // for size_t
#include <cstdint>                           // for int64_t, uint64_t
#include <cstdio>                            // for fclose, FILE, stderr, stdout
#include <exception>                         // for exception
#include <limits>                            // for numeric_limits
#include <stdexcept>                         // for runtime_error, logic_error, range_error
#include <string>                            // for string
#include <string_view>                       // for string_view, basic_string_view, operator==
#include <system_error>                      // for errc, make_error_code, errc::invalid_argument
#include <type_traits>                       // for __strip_reference_wrapper<>::__type, is_arit...
#include <utility>                           // for pair, make_pair, forward
#include <vector>                            // for vector

namespace modle::utils {

template <typename N>
inline auto from_chars(const char *first, const char *last, N &value) noexcept {
#if defined(MODLE_CHARCONV_FP_AVAILABLE) && defined(MODLE_CHARCONV_INT_AVAILABLE)
  return std::from_chars(first, last, value);
#else
  return msstl::from_chars(first, last, value);
#endif
}

template <typename N>
void parse_numeric_or_throw(std::string_view tok, N &field) {
  auto [ptr, err] = utils::from_chars(tok.data(), tok.end(), field);
  if (ptr != tok.end() && err != std::errc{}) {
    throw_except_from_errc(tok, (std::numeric_limits<size_t>::max)(), field, ptr, err);
  }
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_MAYBE_UNINITIALIZED
template <typename N>
N parse_numeric_or_throw(std::string_view tok) {
  N field;
  utils::parse_numeric_or_throw(tok, field);
  return field;
}
DISABLE_WARNING_POP

template <typename N>
void parse_numeric_or_throw(const std::vector<std::string_view> &toks, size_t idx, N &field) {
  parse_numeric_or_throw(toks[idx], field);
}

template <typename N>
void parse_vect_of_numbers_or_throw(const std::vector<std::string_view> &toks, size_t idx,
                                    std::vector<N> &fields, uint64_t expected_size) {
  static_assert(std::is_arithmetic<N>());
  std::vector<std::string_view> ns = absl::StrSplit(toks[idx], ',');
  if (ns.size() != expected_size) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Expected {} fields, got {}."), expected_size, ns.size()));
  }
  fields.resize(ns.size());
  for (size_t i = 0; i < expected_size; ++i) {
    parse_numeric_or_throw(ns, i, fields[i]);
  }
}

template <typename N>
void throw_except_from_errc(std::string_view tok, size_t idx, const N &field, const char *c,
                            std::errc e) {
  (void)field;
  static_assert(std::is_arithmetic<N>());
  std::string base_error;
  if (idx != (std::numeric_limits<size_t>::max)()) {
    base_error = fmt::format(FMT_STRING("Unable to convert field {} ('{}') to a "), idx, tok);
  } else {
    base_error = fmt::format(FMT_STRING("Unable to convert field '{}' to"), tok);
  }
  if (std::is_integral<N>()) {
    if (std::is_unsigned<N>()) {
      base_error += " a positive integral number";
    } else {
      base_error += " an integral number";
    }
  } else {
    base_error += " a real number";
  }
  if (e == std::errc::invalid_argument) {
    if (c != nullptr) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("{}. Reason: found an invalid character '{}'"), base_error, *c));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("{}. Reason: found an invalid character"), base_error));
  }
  if (e == std::errc::result_out_of_range) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("{}. Reason: number {} is outside the range of representable numbers [{}, {}]."),
        base_error, tok, (std::numeric_limits<N>::min)(), (std::numeric_limits<N>::max)()));
  }
  throw std::logic_error(
      fmt::format(FMT_STRING("{}. If you see this error, report it to the developers on "
                             "GitHub.\nBED::throw_except_from_errc "
                             "called with an invalid std::errc '{}'. This should not be possible!"),
                  base_error, std::make_error_code(e).message()));
}

bool chrom_equal_operator(std::string_view chr1, std::string_view chr2) {
  return chrom_equal_operator(std::make_pair(chr1, 0), std::make_pair(chr2, 0));
}

bool chrom_equal_operator(const std::pair<std::string_view, int64_t> &chr1,
                          const std::pair<std::string_view, int64_t> &chr2) {
  if (chr1.second != chr2.second) {
    return false;
  }
  size_t offset1 = 0;
  size_t offset2 = 0;
  if (absl::StartsWithIgnoreCase(chr1.first, "chrom")) {
    offset1 = 3;
  }
  if (absl::StartsWithIgnoreCase(chr2.first, "chrom")) {
    offset2 = 3;
  }
  return chr1.first.substr(offset1) == chr2.first.substr(offset2);
}

bool chrom_less_than_operator(std::string_view chr1, std::string_view chr2) {
  return chrom_less_than_operator(std::make_pair(chr1, 0), std::make_pair(chr2, 0));
}

bool chrom_less_than_operator(const std::pair<std::string_view, int64_t> &chr1,
                              const std::pair<std::string_view, int64_t> &chr2) {
  size_t offset1 = 0;
  size_t offset2 = 0;
  if (absl::StartsWithIgnoreCase(chr1.first, "chrom")) {
    offset1 = 3;
  }
  if (absl::StartsWithIgnoreCase(chr2.first, "chrom")) {
    offset2 = 3;
  }
  if (chr1.first.substr(offset1) < chr2.first.substr(offset2)) {
    return true;
  }

  if (chr1.first.substr(offset1) == chr2.first.substr(offset2)) {
    return chr1.second < chr2.second;
  }
  return false;
}

template <typename T>
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

constexpr bool ndebug_defined() noexcept {
#ifdef NDEBUG
  return true;
#else
  return false;
#endif
}

constexpr bool ndebug_not_defined() noexcept { return !ndebug_defined(); }

void fclose(FILE *fp) noexcept(false) {
  if (!fp || fp == stdout || fp == stderr) {
    return;
  }
  if (std::fclose(fp) != 0) {  // NOLINT
    throw fmt::system_error(errno, FMT_STRING("Failed to close a file handle"));
  }
}

// Try to convert str representations like "1.0" or "1.000000" to "1"
std::string str_float_to_str_int(const std::string &s) {
  try {
    double n = parse_numeric_or_throw<double>(s);
    if (std::trunc(n) == n) {
      return fmt::format(FMT_STRING("{:.0f}"), n);
    }
  } catch (const std::exception &e) {  // Let CLI deal with invalid numbers
    return s;
  }
  return s;
}

std::string detect_path_collision(const boost::filesystem::path &p, bool force_overwrite,
                                  boost::filesystem::file_type expected_type) {
  std::string error_msg;
  detect_path_collision(p, error_msg, force_overwrite, expected_type);
  return error_msg;
}

bool detect_path_collision(const boost::filesystem::path &p, std::string &error_msg,
                           bool force_overwrite, boost::filesystem::file_type expected_type) {
  const auto path_type = boost::filesystem::status(p).type();
  if (force_overwrite && path_type == boost::filesystem::file_not_found) {
    return true;
  }

  if (expected_type != path_type) {
    switch (path_type) {
      case boost::filesystem::regular_file:
        absl::StrAppendFormat(&error_msg,
                              "Path \"%s\" already exists and is actually a file. Please remove "
                              "the file and try again",
                              p.string());
        return false;
      case boost::filesystem::directory_file:
        absl::StrAppendFormat(
            &error_msg,
            "Path \"%s\" already exists and is actually a directory. Please remove "
            "the directory and try again",
            p.string());
        return false;
      default:  // For the time being we only handle regular files and folders
        return true;
    }
  }

  if (force_overwrite) {
    return true;
  }

  if (path_type == boost::filesystem::regular_file) {
    absl::StrAppendFormat(&error_msg, "File \"%s\" already exists. Pass --force to overwrite",
                          p.string());
    return false;
  }
  return true;
}

void XXH3_Deleter::operator()(XXH3_state_t *state) noexcept { XXH3_freeState(state); }

template <class T>
constexpr T &&identity::operator()(T &&a) const noexcept {
  return std::forward<T>(a);
}

template <class Key, class Value, size_t Size>
template <class... Args>
constexpr ConstMap<Key, Value, Size>::ConstMap(Args &&...args) noexcept
    : _buff{{std::forward<Args>(args)...}} {}

template <class Key, class Value, size_t Size>
constexpr const Value &ConstMap<Key, Value, Size>::at(const Key &key) const {
  const auto itr = this->find(key);
  if (itr != this->end()) {
    return itr->second;
  }
  throw std::range_error(fmt::format(FMT_STRING("Unable to find key \"{}\""), key));
}

template <class Key, class Value, size_t Size>
constexpr const Value &ConstMap<Key, Value, Size>::operator[](const Key &key) const {
  return *this->find(key);
}

template <class Key, class Value, size_t Size>
constexpr typename ConstMap<Key, Value, Size>::const_iterator ConstMap<Key, Value, Size>::find(
    const Key &key) const noexcept {
  return std::find_if(this->begin(), this->end(), [&key](const auto &v) { return v.first == key; });
}

template <class Key, class Value, size_t Size>
constexpr bool ConstMap<Key, Value, Size>::contains(const Key &key) const noexcept {
  return this->find(key) != this->_buff.end();
}

template <class Key, class Value, size_t Size>
constexpr typename ConstMap<Key, Value, Size>::const_iterator ConstMap<Key, Value, Size>::begin()
    const noexcept {
  return this->_buff.begin();
}

template <class Key, class Value, size_t Size>
constexpr typename ConstMap<Key, Value, Size>::const_iterator ConstMap<Key, Value, Size>::end()
    const noexcept {
  return this->_buff.end();
}

}  // namespace modle::utils

// IWYU pragma: private, include "modle/utils.hpp"
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <cstdlib>
// IWYU pragma: no_include <memory>
