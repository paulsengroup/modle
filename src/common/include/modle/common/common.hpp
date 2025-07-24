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

// Define short aliases for common integral types
using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using usize = std::size_t;

using i8 = std::int8_t;
using i16 = std::int16_t;
using i32 = std::int32_t;
using i64 = std::int64_t;
using isize = std::ptrdiff_t;

static_assert(sizeof(u8) == 1);
static_assert(sizeof(u16) == 2);
static_assert(sizeof(u32) == 4);
static_assert(sizeof(u64) == 8);
static_assert(sizeof(i8) == 1);
static_assert(sizeof(i16) == 2);
static_assert(sizeof(i32) == 4);
static_assert(sizeof(i64) == 8);

using u8f = std::uint_fast8_t;
using u16f = std::uint_fast16_t;
using u32f = std::uint_fast32_t;
using u64f = std::uint_fast64_t;

using i8f = std::int_fast8_t;
using i16f = std::int_fast16_t;
using i32f = std::int_fast32_t;
using i64f = std::int_fast64_t;

using bp_t = u64;
using contacts_t = u32;

}  // namespace modle
