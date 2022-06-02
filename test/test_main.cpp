// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#define CATCH_CONFIG_MAIN

// clang-format off
#include "modle/common/suppress_compiler_warnings.hpp"
DISABLE_WARNING_PUSH
DISABLE_WARNING_MAYBE_UNINITIALIZED
#include <catch2/catch.hpp>
DISABLE_WARNING_POP
// clang-format on

#include "modle/test/self_deleting_folder.hpp"  // for SelfDeletingFolder

namespace modle::test {
[[maybe_unused]] inline const SelfDeletingFolder testdir{true};  // NOLINT(cert-err58-cpp)
}  // namespace modle::test
