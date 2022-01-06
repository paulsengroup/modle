// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#include "modle/common/smartdir.hpp"

namespace modle::test {
[[maybe_unused]] inline const SmartDir testdir{true};  // NOLINT(cert-err58-cpp)
}  // namespace modle::test
