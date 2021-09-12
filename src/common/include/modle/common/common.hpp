// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>  // IWYU pragma: keep for size_t
#include <cstdint>  // for uint32_t, uint_fast8_t

// These variables are defined by CMake based on the value set through the project() command
#if !defined(MODLE_VERSION_MAJOR) || !defined(MODLE_VERSION_MINOR) || !defined(MODLE_VERSION_PATCH)
#error \
    "Project was not configured properly: MODLE_VERSION_* preprocessor variables do not appear to be defined"
#endif

namespace modle {

using bp_t = size_t;
using collision_t = uint_fast32_t;
using contacts_t = uint32_t;

namespace dna {
enum Direction : uint_fast8_t { none = 0, fwd = 1, rev = 2, both = 3 };
}  // namespace dna

}  // namespace modle
