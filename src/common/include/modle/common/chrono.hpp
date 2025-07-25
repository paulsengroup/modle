// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/chrono.h>
#include <fmt/format.h>

#include <cassert>
#include <chrono>
#include <string>

namespace modle {
namespace internal {

#ifdef __cpp_lib_constexpr_string
[[nodiscard]] constexpr std::string strip_leading_zero(std::string s) {
#else
[[nodiscard]] inline std::string strip_leading_zero(std::string s) {
#endif
  assert(!s.empty());
  if (s.front() != '0') {
    return s;
  }
  return s.substr(1);
}
}  // namespace internal

template <typename Duration>
#ifdef __cpp_lib_constexpr_string
[[nodiscard]] constexpr std::string format_duration(const Duration& duration) {
#else
[[nodiscard]] inline std::string format_duration(const Duration& duration) {
#endif
  if (duration < std::chrono::microseconds(1)) {
    return fmt::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(duration));
  }

  if (duration < std::chrono::milliseconds(1)) {
    return fmt::to_string(std::chrono::duration_cast<std::chrono::microseconds>(duration));
  }

  if (duration < std::chrono::seconds(1)) {
    return fmt::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(duration));
  }

  if (duration < std::chrono::minutes(1)) {
    return internal::strip_leading_zero(
        fmt::format("{:%S}s", std::chrono::duration_cast<std::chrono::milliseconds>(duration)));
  }

  if (duration < std::chrono::hours(1)) {
    return internal::strip_leading_zero(
        fmt::format("{:%Mm:%S}s", std::chrono::duration_cast<std::chrono::milliseconds>(duration)));
  }

  if (duration < std::chrono::days(1)) {
    return internal::strip_leading_zero(fmt::format(
        "{:%Hh:%Mm:%S}s", std::chrono::duration_cast<std::chrono::milliseconds>(duration)));
  }

  return fmt::format("{:%jd:%Hh:%Mm:%S}s",
                     std::chrono::duration_cast<std::chrono::milliseconds>(duration));
}

}  // namespace modle
