#pragma once

// IWYU pragma: private, include "modle/utils.hpp"

#include <absl/strings/match.h>      // for StartsWithIgnoreCase, EndsWith
#include <absl/strings/str_split.h>  // for StrSplit, Splitter
#include <fmt/format.h>              // for format

#include <boost/exception/exception.hpp>    // for enable_error_info, error_info_base...
#include <boost/exception/info.hpp>         // for error_info::name_value_string
#include <boost/stacktrace/stacktrace.hpp>  // for stacktrace, operator<<, to_string
#include <charconv>                         // for from_chars
#include <cmath>                            // for HUGE_VAL
#include <cstddef>                          // IWYU pragma: keep for size_t
#include <cstdint>                          // for int64_t, SIZE_MAX, uint64_t
#include <cstdlib>                          // for strtod
#include <limits>                           // for numeric_limits
#include <memory>                           // for allocator_traits<>::value_type
#include <stdexcept>                        // for runtime_error, logic_error
#include <string>                           // for operator+, string, allocator, basi...
#include <string_view>                      // for string_view, basic_string_view
#include <system_error>                     // for errc, errc::invalid_argument, errc...
#include <type_traits>                      // for __strip_reference_wrapper<>::__type
#include <utility>                          // for pair, make_pair
#include <vector>                           // for vector

namespace modle::utils {

template <typename N>
void parse_numeric_or_throw(std::string_view tok, N &field) {
  static_assert(std::is_arithmetic_v<N>);
  if constexpr (std::is_integral_v<N>) {
    parse_int_or_throw(tok, field);
  }
  if constexpr (std::is_floating_point_v<N>) {
    parse_real_or_throw(tok, field);
  }
}

template <typename I>
void parse_int_or_throw(std::string_view tok, I &field) {
  static_assert(std::is_integral_v<I>);
  auto [ptr, err] = std::from_chars(tok.data(), tok.end(), field);
  if (ptr != tok.end() && err != std::errc{}) {
    throw_except_from_errc(tok, SIZE_MAX, field, ptr, err);
  }
}

template <typename R>
void parse_real_or_throw(std::string_view tok, R &field) {
  static_assert(std::is_floating_point_v<R>);
  char *end = nullptr;
  R tmp = std::strtod(tok.data(), &end);
  if (tmp == HUGE_VAL)
    throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::result_out_of_range);
  else if (tmp == 0 && end == tok.data())
    throw_except_from_errc(tok, SIZE_MAX, tmp, nullptr, std::errc::invalid_argument);

  field = tmp;
}

template <typename N>
void parse_numeric_or_throw(const std::vector<std::string_view> &toks, std::size_t idx, N &field) {
  static_assert(std::is_arithmetic_v<N>);
  if constexpr (std::is_floating_point_v<N>) {
    parse_real_or_throw(toks, idx, field);
  } else {
    parse_int_or_throw(toks, idx, field);
  }
}

template <typename I>
void parse_int_or_throw(const std::vector<std::string_view> &toks, std::size_t idx, I &field) {
  static_assert(std::is_integral<I>());
  auto [ptr, err] = std::from_chars(toks[idx].data(), toks[idx].data() + toks[idx].size(), field);
  if (ptr != toks[idx].end() && err != std::errc{}) {
    throw_except_from_errc(toks[idx], idx, field, ptr, err);
  }
}

template <typename R>
void parse_real_or_throw(const std::vector<std::string_view> &toks, std::size_t idx, R &field) {
  static_assert(std::is_floating_point<R>());
  const std::string tok(toks[idx].begin(), toks[idx].end());
  char *end = nullptr;
  R tmp = std::strtod(tok.data(), &end);
  if (tmp == HUGE_VAL) {
    throw_except_from_errc(toks[idx], idx, tmp, nullptr, std::errc::result_out_of_range);
  } else if (tmp == 0 && end == tok.data()) {
    throw_except_from_errc(tok, idx, tmp, nullptr, std::errc::invalid_argument);
  }

  field = tmp;
}

template <typename N>
void parse_vect_of_numbers_or_throw(const std::vector<std::string_view> &toks, std::size_t idx,
                                    std::vector<N> &field, uint64_t expected_size) {
  static_assert(std::is_arithmetic<N>());
  std::vector<std::string_view> ns = absl::StrSplit(toks[idx], ',');
  if (ns.size() != expected_size) {
    throw std::runtime_error(
        fmt::format("Expected %lu fields, got %lu.", expected_size, ns.size()));
  }
  field = std::vector<N>(ns.size());
  for (std::size_t i = 0; i < expected_size; ++i) {
    parse_numeric_or_throw(ns, i, field[i]);
  }
}

template <typename N>
void throw_except_from_errc(std::string_view tok, std::size_t idx, const N &field, const char *c,
                            std::errc e) {
  (void)field;
  static_assert(std::is_arithmetic<N>());
  std::string base_error;
  if (idx != SIZE_MAX) {
    base_error = fmt::format("Unable to convert field {} ('{}') to a ", idx, tok);
  } else {
    base_error = fmt::format("Unable to convert field '{}' to", tok);
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
          fmt::format("{}. Reason: found an invalid character '{}'", base_error, *c));
    }
    throw std::runtime_error(fmt::format("{}. Reason: found an invalid character", base_error));
  }
  if (e == std::errc::result_out_of_range) {
    throw std::runtime_error(
        fmt::format("{}. Reason: number {} is outside the range of representable numbers [{}, {}].",
                    base_error, tok, std::numeric_limits<N>::min(), std::numeric_limits<N>::max()));
  }
  throw std::logic_error(
      fmt::format("{}. If you see this error, report it to the developers on "
                  "GitHub.\nBED::throw_except_from_errc "
                  "called with an invalid std::errc '{}'. This should not be possible!",
                  base_error, std::make_error_code(e).message()));
}

bool chr_equal_operator(std::string_view chr1, std::string_view chr2) {
  return chr_equal_operator(std::make_pair(chr1, 0), std::make_pair(chr2, 0));
}

bool chr_equal_operator(const std::pair<std::string_view, int64_t> &chr1,
                        const std::pair<std::string_view, int64_t> &chr2) {
  if (chr1.second != chr2.second) {
    return false;
  }
  std::size_t offset1 = 0;
  std::size_t offset2 = 0;
  if (absl::StartsWithIgnoreCase(chr1.first, "chr")) {
    offset1 = 3;
  }
  if (absl::StartsWithIgnoreCase(chr2.first, "chr")) {
    offset2 = 3;
  }
  return chr1.first.substr(offset1) == chr2.first.substr(offset2);
}

bool chr_less_than_operator(std::string_view chr1, std::string_view chr2) {
  return chr_less_than_operator(std::make_pair(chr1, 0), std::make_pair(chr2, 0));
}

bool chr_less_than_operator(const std::pair<std::string_view, int64_t> &chr1,
                            const std::pair<std::string_view, int64_t> &chr2) {
  std::size_t offset1 = 0;
  std::size_t offset2 = 0;
  if (absl::StartsWithIgnoreCase(chr1.first, "chr")) {
    offset1 = 3;
  }
  if (absl::StartsWithIgnoreCase(chr2.first, "chr")) {
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

typedef boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace> traced;  // NOLINT
template <class Except>
void throw_with_trace(const Except &e) {
  throw boost::enable_error_info(e) << traced(boost::stacktrace::stacktrace());
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

}  // namespace modle::utils
