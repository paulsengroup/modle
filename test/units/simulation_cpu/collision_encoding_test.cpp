// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/collision_encoding.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

namespace modle::libmodle::test {

[[maybe_unused]] static std::string event_to_str(const CollisionEvent<> event) {
  return fmt::format(FMT_STRING("{:08b}"), event());
}

[[maybe_unused]] constexpr usize EVENT_BITS = 5;
[[maybe_unused]] constexpr usize RESERVED_EVENT_BITS = 8;
[[maybe_unused]] constexpr usize INDEX_BITS = (sizeof(u32f) * 8) - RESERVED_EVENT_BITS;
[[maybe_unused]] constexpr u32f EVENT_MASK = (~u32f(0)) << (INDEX_BITS - 1);
[[maybe_unused]] constexpr u32f INDEX_MASK = ~EVENT_MASK;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Collision encoding", "[simulation][short]") {
  usize idx = 5;
  // clang-format off
  const std::array<CollisionEvent<u8f>, 8> events{
        Collision<>::LEF_BAR,
        Collision<>::LEF_LEF_PRIMARY,
        Collision<>::LEF_LEF_SECONDARY,
        Collision<>::CHROM_BOUNDARY,
        Collision<>::COLLISION | Collision<>::LEF_BAR,
        Collision<>::COLLISION | Collision<>::LEF_LEF_PRIMARY,
        Collision<>::COLLISION | Collision<>::LEF_LEF_SECONDARY,
        Collision<>::COLLISION | Collision<>::CHROM_BOUNDARY
  };
  // clang-format on

  for (const auto& event : events) {
    const auto encoded_collision = Collision<>(idx, event);
    CHECK(encoded_collision.decode_index() == idx);
    CHECK(encoded_collision.decode_event() == event);
  }

  if constexpr (utils::ndebug_not_defined()) {
    // clang-format off
    const std::array<CollisionEvent<u8f>, 9> invalid_events{
          Collision<>::COLLISION,
          Collision<>::COLLISION | Collision<>::LEF_BAR | Collision<>::LEF_LEF_PRIMARY,
          Collision<>::LEF_BAR | Collision<>::LEF_LEF_PRIMARY
    };
    // clang-format on

    for (const auto& event : invalid_events) {
      // For some reason Catch2 is not catching the exceptions thrown by Collision<> ctor.
      // This ugly workaround will have to do for the time being
      try {
        [[maybe_unused]] volatile const Collision<> c(idx, event);
      } catch (const std::runtime_error& e) {
        CHECK_THROWS_WITH(throw e, Catch::Matchers::StartsWith("Collision<I>::encode"));
      }
    }

    CollisionEvent<u8f> event{};
    CHECK_THROWS_WITH(
        (Collision<>{idx, event}.collision_occurred()),
        Catch::Matchers::EndsWith("idx must be 0 when no collision has occurred/been avoided"));
    CHECK_THROWS_WITH(
        (Collision<>{idx, event}.collision_avoided()),
        Catch::Matchers::EndsWith("idx must be 0 when no collision has occurred/been avoided"));
  }

  Collision<> collision{};
  collision.set_event(Collision<>::COLLISION | Collision<>::LEF_BAR);
  REQUIRE(collision.collision_occurred(Collision<>::LEF_BAR));

  collision.set_event(Collision<>::LEF_LEF_PRIMARY);
  CHECK(collision.collision_avoided(Collision<>::LEF_LEF_PRIMARY));

  idx = 123;
  collision.set_idx(idx);
  CHECK(collision.collision_avoided(Collision<>::LEF_LEF_PRIMARY));
  CHECK(collision.decode_index() == idx);

  idx = INDEX_MASK;
  collision.set(idx, Collision<>::COLLISION | Collision<>::LEF_LEF_SECONDARY);
  CHECK(collision.collision_occurred(Collision<>::LEF_LEF_SECONDARY));
  CHECK(collision.decode_index() == idx);

  collision.set(idx, Collision<>::COLLISION | Collision<>::LEF_BAR);
  CHECK(collision.collision_occurred(Collision<>::LEF_BAR));
  CHECK(collision.decode_index() == idx);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Collision encoding - Chromosomal boundaries", "[simulation][short]") {
  for (const auto idx : std::array<usize, 2>{3, 5}) {
    auto event = Collision<>::COLLISION | Collision<>::CHROM_BOUNDARY;
    CHECK(Collision<>{idx, event}.collision_occurred());
    CHECK(Collision<>{idx, event}.collision_occurred(Collision<>::CHROM_BOUNDARY));

    CHECK_FALSE(Collision<>{idx, event}.collision_avoided());
    CHECK_FALSE(Collision<>{idx, event}.collision_avoided(Collision<>::CHROM_BOUNDARY));

    event = Collision<>::CHROM_BOUNDARY;
    CHECK_FALSE(Collision<>{idx, event}.collision_occurred());
    CHECK_FALSE(Collision<>{idx, event}.collision_occurred(Collision<>::CHROM_BOUNDARY));

    CHECK(Collision<>{idx, event}.collision_avoided());
    CHECK(Collision<>{idx, event}.collision_avoided(Collision<>::CHROM_BOUNDARY));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Collision encoding - LEF-BAR", "[simulation][short]") {
  const usize idx = 123;

  CollisionEvent<u8f> event = Collision<>::COLLISION() | Collision<>::LEF_BAR();
  CHECK(Collision<>{idx, event}.collision_occurred());
  CHECK(Collision<>{idx, event}.collision_occurred(Collision<>::LEF_BAR));

  CHECK_FALSE(Collision<>{idx, event}.collision_avoided());
  CHECK_FALSE(Collision<>{idx, event}.collision_avoided(Collision<>::LEF_BAR));

  event = Collision<>::LEF_BAR();
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred());
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred(Collision<>::LEF_BAR));

  CHECK(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.collision_avoided(Collision<>::LEF_BAR));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Collision encoding - LEF-LEF primary", "[simulation][short]") {
  const usize idx = 123;

  CollisionEvent<u8f> event = Collision<>::COLLISION() | Collision<>::LEF_LEF_PRIMARY();
  CHECK(Collision<>{idx, event}.collision_occurred());
  CHECK(Collision<>{idx, event}.collision_occurred(Collision<>::LEF_LEF_PRIMARY));

  CHECK_FALSE(Collision<>{idx, event}.collision_avoided());
  CHECK_FALSE(Collision<>{idx, event}.collision_avoided(Collision<>::LEF_LEF_PRIMARY));

  event = Collision<>::LEF_LEF_PRIMARY();
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred());
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred(Collision<>::LEF_LEF_PRIMARY));

  CHECK(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.collision_avoided(Collision<>::LEF_LEF_PRIMARY));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Collision encoding - LEF-LEF secondary", "[simulation][short]") {
  const usize idx = 123;

  CollisionEvent<u8f> event = Collision<>::COLLISION() | Collision<>::LEF_LEF_SECONDARY();
  CHECK(Collision<>{idx, event}.collision_occurred());
  CHECK(Collision<>{idx, event}.collision_occurred(Collision<>::LEF_LEF_SECONDARY));

  CHECK_FALSE(Collision<>{idx, event}.collision_avoided());
  CHECK_FALSE(Collision<>{idx, event}.collision_avoided(Collision<>::LEF_LEF_SECONDARY));

  event = Collision<>::LEF_LEF_SECONDARY();
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred());
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred(Collision<>::LEF_LEF_SECONDARY));

  CHECK(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.collision_avoided(Collision<>::LEF_LEF_SECONDARY));
}

}  // namespace modle::libmodle::test
