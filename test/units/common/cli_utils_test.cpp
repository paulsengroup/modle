// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch.hpp>  // for operator""_catch_sr, AssertionHandler
// clang-format off
#include "modle/common/suppress_compiler_warnings.hpp"
DISABLE_WARNING_PUSH
DISABLE_WARNING_SHADOW
#include <bitflags/bitflags.hpp>
DISABLE_WARNING_POP
// clang-format on

#include "modle/common/cli_utils.hpp"

namespace modle::test::utils {

BEGIN_BITFLAGS(TestEnum)
FLAG(a)
FLAG(b)
FLAG(c)
FLAG(d)
END_BITFLAGS(TestEnum)

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CliEnumMappings", "[utils][short]") {
  const ::modle::utils::CliEnumMappings<TestEnum> test_enum_map{
      // clang-format off
      std::make_pair("a",               TestEnum::a),
      std::make_pair("b",               TestEnum::b),
      std::make_pair("a-plus-b",        TestEnum::a  | TestEnum::b),
      std::make_pair("a-plus-c",        TestEnum::a  | TestEnum::c),
      std::make_pair("b-plus-c",        TestEnum::b  | TestEnum::c),
      std::make_pair("a-plus-b-plus-c", TestEnum::a  | TestEnum::b  | TestEnum::c)
      // clang-format on
  };

  CHECK(test_enum_map.at("a") == (TestEnum::a));
  CHECK(test_enum_map.at("b") == (TestEnum::b));
  CHECK(test_enum_map.at("a-plus-b") == (TestEnum::a | TestEnum::b));
  CHECK(test_enum_map.at("a-plus-c") == (TestEnum::a | TestEnum::c));
  CHECK(test_enum_map.at("b-plus-c") == (TestEnum::b | TestEnum::c));
  CHECK(test_enum_map.at("a-plus-b-plus-c") == (TestEnum::a | TestEnum::b | TestEnum::c));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("AsGenomicDistance CLI transform operator", "[utils][short]") {
  auto fx = [&](std::string&& s) {
    modle::utils::cli::AsGenomicDistance(s);
    return s;
  };

  CHECK(fx("10000") == "10000");
  CHECK(fx("10kb") == "10000");
  CHECK(fx("10KB") == "10000");
  CHECK(fx("10KBP") == "10000");
  CHECK(fx("10KbP") == "10000");
  CHECK(fx("10k") == "10000");

  CHECK(fx("5.5kb") == "5500");
  CHECK(fx("0.5kb") == "500");
  CHECK(fx("1mbp") == "1000000");

  CHECK(fx("1e6kb") == fx("1gbp"));

  CHECK_THROWS_AS(fx("10kbpp"), CLI::ValidationError);
  CHECK_THROWS_AS(fx("0.0001kb"), CLI::ValidationError);
  CHECK_THROWS_AS(fx("!10gbp"), CLI::ValidationError);
  CHECK_THROWS_AS(fx(" 10gbp"), CLI::ValidationError);
  CHECK_THROWS_AS(fx(""), CLI::ValidationError);
  CHECK_THROWS_AS(fx("kb"), CLI::ValidationError);
}

}  // namespace modle::test::utils
