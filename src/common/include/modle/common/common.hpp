// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>  // IWYU pragma: keep
#include <cstdint>
#include <stdexcept>
#include <type_traits>

namespace modle {

// NOLINTBEGIN(cppcoreguidelines-macro-usage)

// See https://clang.llvm.org/docs/ThreadSanitizer.html#has-feature-thread-sanitizer
#define MODLE_NO_TSAN
#if defined(__has_feature)
#if __has_feature(thread_sanitizer)
#undef MODLE_NO_TSAN  // NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MODLE_NO_TSAN __attribute__((no_sanitize("thread")))
#endif
#endif

#if defined(__GNUC__)
#define MODLE_UNREACHABLE_CODE __builtin_unreachable()
#elif defined(_MSC_VER)
#define MODLE_UNREACHABLE_CODE __assume(0)
#else
#define MODLE_UNREACHABLE_CODE
#endif

// NOLINTEND(cppcoreguidelines-macro-usage)

namespace utils {
[[maybe_unused]] [[nodiscard]] constexpr bool ndebug_defined() noexcept {
#ifdef NDEBUG
  return true;
#else
  return false;
#endif
}

[[maybe_unused]] [[nodiscard]] constexpr bool ndebug_not_defined() noexcept {
  return !ndebug_defined();
}

[[noreturn]] inline void unreachable_code() {
  if constexpr (ndebug_not_defined()) {
    throw std::logic_error("Unreachable code");
  }
  MODLE_UNREACHABLE_CODE;
}

// to avoid useless casts (see https://github.com/nlohmann/json/issues/2893#issuecomment-889152324)
template <class T, class U>
[[maybe_unused]] [[nodiscard]] constexpr T conditional_static_cast(U value) {
  if constexpr (std::is_same_v<T, U>) {
    return value;
  } else {
    return static_cast<T>(value);
  }
}

}  // namespace utils

using bp_t = std::uint64_t;
using contacts_t = std::uint32_t;

}  // namespace modle
