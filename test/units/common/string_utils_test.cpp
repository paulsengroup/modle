// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/common/string_utils.hpp"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <string_view>
#include <vector>

namespace modle::utils::test {

TEST_CASE("to_lower()", "[utils][short]") {
  CHECK(to_lower("") == "");
  CHECK(to_lower("1") == "1");
  CHECK(to_lower("A") == "a");
  CHECK(to_lower("a") == "a");
  CHECK(to_lower("aBc123!\0") == "abc123!\0");
}

TEST_CASE("strip_trailing_whitespace", "[utils][short]") {
  CHECK(strip_trailing_whitespace("") == "");
  CHECK(strip_trailing_whitespace(" ") == "");
  CHECK(strip_trailing_whitespace("\t") == "");
  CHECK(strip_trailing_whitespace("\n") == "");
  CHECK(strip_trailing_whitespace("  ") == "");
  CHECK(strip_trailing_whitespace(" a") == " a");
  CHECK(strip_trailing_whitespace("a ") == "a");
  CHECK(strip_trailing_whitespace("a  ") == "a");
  CHECK(strip_trailing_whitespace("a b  ") == "a b");
}

TEST_CASE("str_contains", "[utils][short]") {
  CHECK_FALSE(str_contains("", ""));
  CHECK_FALSE(str_contains("", "a"));
  CHECK_FALSE(str_contains("a", ""));

  CHECK(str_contains("abc", "ab"));
  CHECK(str_contains("abc", "bc"));
  CHECK_FALSE(str_contains("abc", "ac"));
  CHECK_FALSE(str_contains("abc", "abcd"));
}

TEST_CASE("str_strip_suffix", "[utils][short]") {
  CHECK(str_strip_suffix("", "") == "");
  CHECK(str_strip_suffix("", "a") == "");
  CHECK(str_strip_suffix("a", "a") == "");
  CHECK(str_strip_suffix("abc", "ab") == "abc");
  CHECK(str_strip_suffix("abc", "bc") == "a");
  CHECK(str_strip_suffix("abc", "ab") == "abc");
  CHECK(str_strip_suffix("abc", "def") == "abc");
}

TEST_CASE("str_split", "[utils][short]") {
  auto test_str_split = [](std::string_view a, const std::vector<std::string_view>& b,
                           std::string_view delim = " ") {
    return std::ranges::equal(str_split(a, delim), b);
  };

  CHECK(str_split("", ' ').empty());
  CHECK(str_split("", "").empty());
  CHECK(std::ranges::equal(str_split("abc", ""), std::vector{"abc"}));
  CHECK(test_str_split("", {}));
  CHECK(test_str_split("abc", {"abc"}));
  CHECK(test_str_split("a bc", {"a", "bc"}));
  CHECK(test_str_split("a  bc", {"a", "", "bc"}));
  CHECK(test_str_split("   ", {"", "", ""}));
  CHECK(test_str_split("a b\tc", {"a", "b", "c"}, " \t"));
}

TEST_CASE("str_replace", "[utils][short]") {
  CHECK(str_replace("", "", "") == "");
  CHECK(str_replace("", "a", "b") == "");
  CHECK(str_replace("", "", "b") == "");
  CHECK(str_replace("", "a", "") == "");
  CHECK(str_replace("abc", "x", "b") == "abc");
  CHECK(str_replace("abc", "a", "b") == "bbc");
  CHECK(str_replace("abc", "ab", "b") == "bc");
  CHECK(str_replace("abc", "ab", "") == "c");
  CHECK(str_replace("abcfooabc", "abc", "bar") == "barfoobar");
}

}  // namespace modle::utils::test
