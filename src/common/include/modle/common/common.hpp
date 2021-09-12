// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for uint32_t, uint_fast8_t
#include <string_view>  // for string_view

// These variables are defined by CMake based on the value set through the project() command
#if !defined(MODLE_VERSION_MAJOR) || !defined(MODLE_VERSION_MINOR) || !defined(MODLE_VERSION_PATCH)
#error \
    "Project was not configured properly: MODLE_VERSION_* preprocessor variables do not appear to be defined"
#endif

namespace modle {

static_assert(MODLE_VERSION_MAJOR >= 0);
static_assert(MODLE_VERSION_MINOR >= 0);
static_assert(MODLE_VERSION_PATCH >= 0);
constexpr uint8_t modle_version_major{MODLE_VERSION_MAJOR};
constexpr uint8_t modle_version_minor{MODLE_VERSION_MINOR};
constexpr uint8_t modle_version_patch{MODLE_VERSION_PATCH};

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

using bp_t = size_t;
using collision_t = uint_fast32_t;
using contacts_t = uint32_t;

namespace dna {
enum Direction : uint_fast8_t { none = 0, fwd = 1, rev = 2, both = 3 };
}  // namespace dna

}  // namespace modle
