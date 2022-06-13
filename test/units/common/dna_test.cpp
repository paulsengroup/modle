// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/common/dna.hpp"

#include <catch2/catch.hpp>

#include "modle/common/common.hpp"

namespace modle::test::dna {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("DNA Ctor", "[common][short]") {
  namespace dna = modle::dna;
  CHECK(dna::Direction(u8f(0)) == dna::NONE);
  CHECK(dna::Direction(u8f(1)) == dna::REV);
  CHECK(dna::Direction(u8f(2)) == dna::FWD);

  CHECK(dna::Direction::from_char('.') == dna::NONE);
  CHECK(dna::Direction::from_char('-') == dna::REV);
  CHECK(dna::Direction::from_char('+') == dna::FWD);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("DNA complement", "[common][short]") {
  namespace dna = modle::dna;
  CHECK(dna::NONE.complement() == dna::NONE);
  CHECK(dna::REV.complement() == dna::FWD);
  CHECK(dna::FWD.complement() == dna::REV);
}

}  // namespace modle::test::dna
