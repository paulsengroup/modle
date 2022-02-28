// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// clang-format off
#if defined(MODLE_CHARCONV_FP_AVAILABLE) && defined(MODLE_CHARCONV_INT_AVAILABLE)
#include <charconv>                                     // for from_chars
#else
#include "modle/common/suppress_compiler_warnings.hpp"  // IWYU pragma: keep
DISABLE_WARNING_PUSH
DISABLE_WARNING_CONVERSION
DISABLE_WARNING_DOUBLE_PROMOTION
DISABLE_WARNING_SHORTEN_64_TO_32
DISABLE_WARNING_SIGN_CONVERSION
#include <msstl/charconv.hpp>

DISABLE_WARNING_POP
#endif
// clang-format on

#include <absl/strings/str_split.h>  // for StrSplit, Splitter
#include <fmt/format.h>              // for compile_string_to_view, FMT_STRING

#include <limits>        // for numeric_limits
#include <stdexcept>     // for runtime_error, logic_error
#include <string>        // for string
#include <string_view>   // for string_view
#include <system_error>  // for errc, make_error_code, errc::invalid_argument, errc:...
#include <type_traits>   // for is_arithmetic, is_integral, is_unsigned
#include <vector>        // for vector

#include "modle/common/common.hpp"  // for usize, u64

namespace modle::utils {

template <class N>
inline auto from_chars(const char *first, const char *last, N &value) noexcept {
#if defined(MODLE_CHARCONV_FP_AVAILABLE) && defined(MODLE_CHARCONV_INT_AVAILABLE)
  return std::from_chars(first, last, value);
#else
  return msstl::from_chars(first, last, value);
#endif
}

template <class N>
void parse_numeric_or_throw(std::string_view tok, N &field) {
  auto [ptr, err] = utils::from_chars(tok.data(), tok.end(), field);
  if (ptr != tok.end() && err != std::errc{}) {
    utils::detail::throw_except_from_errc(tok, (std::numeric_limits<usize>::max)(), field, ptr,
                                          err);
  }
}

template <class N>
N parse_numeric_or_throw(std::string_view tok) {
  N field{};
  utils::parse_numeric_or_throw(tok, field);
  return field;
}

template <class N>
void parse_numeric_or_throw(const std::vector<std::string_view> &toks, usize idx, N &field) {
  parse_numeric_or_throw(toks[idx], field);
}

template <class N>
void parse_vect_of_numbers_or_throw(const std::vector<std::string_view> &toks, usize idx,
                                    std::vector<N> &fields, u64 expected_size) {
  static_assert(std::is_arithmetic<N>());
  std::vector<std::string_view> ns = absl::StrSplit(toks[idx], ',');
  if (ns.size() != expected_size) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Expected {} fields, got {}."), expected_size, ns.size()));
  }
  fields.resize(ns.size());
  for (usize i = 0; i < expected_size; ++i) {
    parse_numeric_or_throw(ns, i, fields[i]);
  }
}

namespace detail {
template <class N>
void throw_except_from_errc(std::string_view tok, usize idx, [[maybe_unused]] const N &field,
                            const char *c, std::errc e) {
  static_assert(std::is_arithmetic<N>());
  std::string base_error;
  if (idx != (std::numeric_limits<usize>::max)()) {
    base_error = fmt::format(FMT_STRING("Unable to convert field {} (\"{}\") to a "), idx, tok);
  } else {
    base_error = fmt::format(FMT_STRING("Unable to convert field \"{}\" to"), tok);
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
          fmt::format(FMT_STRING("{}. Reason: found an invalid character \"{}\""), base_error, *c));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("{}. Reason: found an invalid character"), base_error));
  }
  if (e == std::errc::result_out_of_range) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("{}. Reason: number {} is outside the range of representable numbers [{}, {}]."),
        base_error, tok, (std::numeric_limits<N>::min)(), (std::numeric_limits<N>::max)()));
  }

  throw std::logic_error(fmt::format(
      FMT_STRING("{}. If you see this error, report it to the developers on "
                 "GitHub.\n throw_except_from_errc "
                 "called with an invalid std::errc \"{}\". This should not be possible!"),
      base_error, std::make_error_code(e).message()));
}
}  // namespace detail

}  // namespace modle::utils

// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <cstdlib>
// IWYU pragma: no_include <memory>
