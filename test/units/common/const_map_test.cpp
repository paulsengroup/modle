// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/common/const_map.hpp"

#include <catch2/catch_test_macros.hpp>
#include <string_view>

namespace modle::utils::test {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ConstMap Ctor", "[utils][short]") {
  using CMap = ConstMap<std::string_view, int, 5>;
  constexpr std::array<CMap::key_type, 5> keys{"1", "2", "3", "4", "5"};
  constexpr std::array<CMap::mapped_type, 5> vals{1, 2, 3, 4, 5};

  CHECK_NOTHROW(CMap{{"1", 1}, {"2", 2}, {"3", 3}, {"4", 4}, {"5", 5}});
  CHECK_NOTHROW(CMap(keys.begin(), keys.end(), vals.begin()));

  CHECK_THROWS(CMap{{"1", 1}, {"2", 2}, {"3", 3}, {"4", 4}});
  CHECK_THROWS(CMap{{"1", 1}, {"2", 2}, {"3", 3}, {"4", 4}, {"5", 5}, {"6", 6}});
  CHECK_THROWS(CMap(keys.begin(), keys.end() - 1, vals.begin()));
  CHECK_THROWS(CMap(keys.begin(), keys.end() + 1, vals.begin()));

  if constexpr (ndebug_not_defined()) {
    CHECK_THROWS(CMap{{"1", 1}, {"1", 2}, {"3", 3}, {"4", 4}, {"5", 5}});
  } else {
    CHECK_NOTHROW(CMap{{"1", 1}, {"1", 2}, {"3", 3}, {"4", 4}, {"5", 5}});
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ConstMap Accessors", "[utils][short]") {
  using CMap = ConstMap<std::string_view, std::size_t, 5>;
  constexpr CMap map{{"1", 1}, {"2", 2}, {"3", 3}, {"4", 4}, {"5", 5}};
  for (std::size_t i = 1; i <= map.size(); ++i) {
    CHECK(map.contains(std::to_string(i)));
    CHECK(map.at(std::to_string(i)) == i);
    const auto it = map.find(std::to_string(i));
    CHECK(*it.first == std::to_string(i));
    CHECK(*it.second == i);
  }

  CHECK(!map.contains("0"));

  const auto& keys = map.keys();
  const auto& values = map.values();
  for (std::size_t i = 0; i < keys.size(); ++i) {
    CHECK(keys[i] == std::to_string(i + 1));
    CHECK(values[i] == i + 1);
  }
}

}  // namespace modle::utils::test
