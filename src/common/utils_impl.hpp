#pragma once

// IWYU pragma: private, include "modle/utils.hpp"
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <modle/src/utils/utils_impl.hpp>

#if defined(MODLE_CHARCONV_FP_AVAILABLE) || defined(MODLE_CHARCONV_INT_AVAILABLE)
#include <charconv>  // for from_chars
#endif
#ifndef MODLE_CHARCONV_FP_AVAILABLE
#include <absl/strings/charconv.h>  // for from_chars (FP)
#endif

#include <absl/strings/match.h>       // for StartsWithIgnoreCase, EndsWith
#include <absl/strings/str_format.h>  // for absl::StrAppendFormat
#include <absl/strings/str_split.h>   // for StrSplit, Splitter
#include <fmt/format.h>               // for format

#include <boost/filesystem.hpp>  // for path, status
#include <cassert>               // for assert
#include <cmath>                 // for HUGE_VAL
#include <cstddef>               // IWYU pragma: keep for size_t
#include <cstdint>               // for int64_t, SIZE_MAX, uint64_t
#include <cstdio>                // for fclose, FILE
#include <cstdlib>               // for strtod
#include <limits>                // for numeric_limits
#include <memory>                // for allocator_traits<>::value_type
#include <stdexcept>             // for runtime_error, logic_error
#include <string>                // for operator+, string, allocator, basi...
#include <string_view>           // for string_view, basic_string_view
#include <system_error>          // for errc, errc::invalid_argument, errc...
#include <type_traits>           // for __strip_reference_wrapper<>::__type
#include <utility>               // for pair, make_pair
#include <vector>                // for vector

#include "modle/common/suppress_compiler_warnings.hpp"

namespace modle::utils {

template <typename N, typename>
void parse_numeric_or_throw(std::string_view tok, N &field) {
  if constexpr (std::is_integral_v<N>) {
#ifdef MODLE_CHARCONV_INT_AVAILABLE  // str -> integral
    auto [ptr, err] = std::from_chars(tok.data(), tok.end(), field);
    if (ptr != tok.end() && err != std::errc{}) {
      throw_except_from_errc(tok, SIZE_MAX, field, ptr, err);
    }
#else
    if constexpr (std::is_same_v<N, int32_t>) {  // str -> int32
      const auto tmp = static_cast<N>(std::stol(tok.data(), nullptr));
      if (tmp == LONG_MAX) {
        throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::result_out_of_range);
      } else if (tmp == 0 && tok != "0") {
        throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::invalid_argument);
      }
      field = tmp;
    } else if constexpr (std::is_same_v<N, int64_t>) {  // str -> int64
      const auto tmp = static_cast<N>(std::stoll(tok.data(), nullptr));
      if (tmp == LONG_LONG_MAX) {
        throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::result_out_of_range);
      } else if (tmp == 0 && tok != "0") {
        throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::invalid_argument);
      }
      field = tmp;
    } else if constexpr (std::is_same_v<N, uint32_t>) {  // str -> uint32
      const auto tmp = static_cast<N>(std::stoul(tok.data(), nullptr));
      if (tmp == ULONG_MAX) {
        throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::result_out_of_range);
      } else if (tmp == 0 && tok != "0") {
        throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::invalid_argument);
      }
      field = tmp;
    } else if constexpr (std::is_same_v<N, uint64_t>) {  // str -> uint64
      const auto tmp = static_cast<N>(std::stoull(tok.data(), nullptr));
      if (tmp == ULONG_LONG_MAX) {
        throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::result_out_of_range);
      } else if (tmp == 0 && tok != "0") {
        throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::invalid_argument);
      }
      field = tmp;
    } else {  // str -> other int
      const auto tmp = std::stoll(tok.data(), nullptr);
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_USELESS_CAST
      if (tmp == LONG_LONG_MAX || tmp < static_cast<int64_t>((std::numeric_limits<N>::min)()) ||
          tmp > static_cast<int64_t>((std::numeric_limits<N>::max)())) {
        throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::result_out_of_range);
      } else if (tmp == 0 && tok != "0") {
        throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::invalid_argument);
      }
      DISABLE_WARNING_POP
      field = static_cast<N>(tmp);
    }
#endif
    // str -> floating point
  } else if constexpr (std::is_floating_point_v<N>) {
#ifdef MODLE_CHARCONV_FP_AVAILABLE
    auto [ptr, err] = std::from_chars(tok.data(), tok.end(), field);
#else
    auto [ptr, err] = absl::from_chars(tok.data(), tok.end(), field);
#endif
    if (ptr != tok.end() && err != std::errc{}) {
      throw_except_from_errc(tok, SIZE_MAX, field, ptr, err);
    }
  }
}

template <typename N>
void parse_numeric_or_throw(const std::vector<std::string_view> &toks, size_t idx, N &field) {
  parse_numeric_or_throw(toks[idx], field);
}

template <typename N>
void parse_vect_of_numbers_or_throw(const std::vector<std::string_view> &toks, size_t idx,
                                    std::vector<N> &field, uint64_t expected_size) {
  static_assert(std::is_arithmetic<N>());
  std::vector<std::string_view> ns = absl::StrSplit(toks[idx], ',');
  if (ns.size() != expected_size) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Expected {} fields, got {}."), expected_size, ns.size()));
  }
  field = std::vector<N>(ns.size());
  for (size_t i = 0; i < expected_size; ++i) {
    parse_numeric_or_throw(ns, i, field[i]);
  }
}

template <typename N>
void throw_except_from_errc(std::string_view tok, size_t idx, const N &field, const char *c,
                            std::errc e) {
  (void)field;
  static_assert(std::is_arithmetic<N>());
  std::string base_error;
  if (idx != SIZE_MAX) {
    base_error = fmt::format(FMT_STRING("Unable to convert field {} ('{}') to a "), idx, tok);
  } else {
    base_error = fmt::format(FMT_STRING("Unable to convert field '{}' to"), tok);
  }
  if (std::is_integral<N>()) {
    if (std::is_unsigned<N>()) {
      base_error += " a positive integral number";
    } else {
      base_error += "an integral number";
    }
  } else {
    base_error += "a real number";
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
    double n;
    parse_numeric_or_throw(s, n);
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
}  // namespace modle::utils
