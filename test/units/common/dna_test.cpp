// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/common/dna.hpp"

#include <catch2/catch_test_macros.hpp>

#include "modle/common/common.hpp"

namespace modle::dna::test {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("DNA Ctor", "[common][short]") {
  CHECK(Direction(u8f(0)) == NONE);
  CHECK(Direction(u8f(1)) == REV);
  CHECK(Direction(u8f(2)) == FWD);

  CHECK(Direction::from_char('.') == NONE);
  CHECK(Direction::from_char('-') == REV);
  CHECK(Direction::from_char('+') == FWD);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("DNA complement", "[common][short]") {
  CHECK(NONE.complement() == NONE);
  CHECK(REV.complement() == FWD);
  CHECK(FWD.complement() == REV);
}

}  // namespace modle::dna::test
