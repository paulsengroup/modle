// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/base/macros.h>  // for ABSL_INTERNAL_UNREACHABLE

#include <cstddef>      // IWYU pragma: keep for size_t, ptrdiff_t
#include <cstdint>      // for i8, i16, i32, i64, u8 ...
#include <string_view>  // for string_view

// These variables are defined by CMake based on the value set through the project() command
#if !defined(MODLE_VERSION_MAJOR) || !defined(MODLE_VERSION_MINOR) || !defined(MODLE_VERSION_PATCH)
#error \
    "Project was not configured properly: MODLE_VERSION_* preprocessor variables do not appear to be defined"
#endif

namespace modle {

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

static_assert(MODLE_VERSION_MAJOR >= 0);
static_assert(MODLE_VERSION_MINOR >= 0);
static_assert(MODLE_VERSION_PATCH >= 0);
constexpr u8 modle_version_major{MODLE_VERSION_MAJOR};
constexpr u8 modle_version_minor{MODLE_VERSION_MINOR};
constexpr u8 modle_version_patch{MODLE_VERSION_PATCH};

// String macros
#define MODLE_XSTR(x)               MODLE_STR(x)
#define MODLE_STR(x)                #x
#define MODLE_CONCAT(x, y, z, sep)  x##sep##y##sep##z
#define MODLE_XCONCAT(x, y, z, sep) MODLE_CONCAT(x, y, z, sep)

#define MODLE_VERSION \
  MODLE_XCONCAT(MODLE_VERSION_MAJOR, MODLE_VERSION_MINOR, MODLE_VERSION_PATCH, .)

constexpr std::string_view modle_version_long{"MoDLE-v" MODLE_XSTR(MODLE_VERSION)};
constexpr std::string_view modle_version =
    modle_version_long.substr(std::string_view{"MoDLE-v"}.size());
#undef MODLE_VERSION
#undef MODLE_XSTR
#undef MODLE_STR
#undef MODLE_CONCAT
#undef MODLE_XCONCAT

using bp_t = u64;
using collision_t = std::uint_fast32_t;  // TODO removeme
using contacts_t = u32;

// See https://clang.llvm.org/docs/ThreadSanitizer.html#has-feature-thread-sanitizer
#define MODLE_NO_TSAN
#if defined(__has_feature)
#if __has_feature(thread_sanitizer)
#undef MODLE_NO_TSAN
#define MODLE_NO_TSAN __attribute__((no_sanitize("thread")))
#endif
#endif

#define MODLE_UNREACHABLE_CODE ABSL_INTERNAL_UNREACHABLE

namespace dna {
enum Direction : u8f { none = 0, fwd = 1, rev = 2, both = 3 };
}  // namespace dna

}  // namespace modle
