// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#include "modle/test/self_deleting_folder.hpp"  // for SelfDeletingFolder

namespace modle::test {
[[maybe_unused]] inline const SelfDeletingFolder testdir{true};  // NOLINT(cert-err58-cpp)
}  // namespace modle::test
