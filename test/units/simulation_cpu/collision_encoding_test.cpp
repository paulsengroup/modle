// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/collision_encoding.hpp"

#include <fmt/format.h>

#include <catch2/catch.hpp>  // for SourceLineInfo, StringRef, TEST_CASE
namespace modle::test::libmodle {

[[maybe_unused]] static std::string event_to_str(const CollisionEvent<> event) {
  return fmt::format(FMT_STRING("{:08b}"), event());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Collision encoding", "[simulation][short]") {
  const usize idx = 123;
  // clang-format off
  constexpr std::array<CollisionEvent<u8f>, 8> events{
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
    constexpr std::array<CollisionEvent<u8f>, 9> invalid_events{
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
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Collision encoding - Chromosomal boundaries", "[simulation][short]") {
  usize idx = 3;  // NOLINT
  CollisionEvent<u8f> event = Collision<>::COLLISION() | Collision<>::CHROM_BOUNDARY();
  CHECK(Collision<>{idx, event}.collision_occurred());
  CHECK_FALSE(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.is_chrom_boundary_collision());
  CHECK(Collision<>{idx, event}.is_chrom_boundary_collision_3p());
  CHECK_FALSE(Collision<>{idx, event}.is_chrom_boundary_collision_5p());

  event = Collision<>::CHROM_BOUNDARY();
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred());
  CHECK(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.is_chrom_boundary_collision());
  CHECK(Collision<>{idx, event}.is_chrom_boundary_collision_3p());
  CHECK_FALSE(Collision<>{idx, event}.is_chrom_boundary_collision_5p());

  idx = 5;  // NOLINT
  event = Collision<>::COLLISION() | Collision<>::CHROM_BOUNDARY();
  CHECK(Collision<>{idx, event}.collision_occurred());
  CHECK_FALSE(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.is_chrom_boundary_collision());
  CHECK(Collision<>{idx, event}.is_chrom_boundary_collision_5p());
  CHECK_FALSE(Collision<>{idx, event}.is_chrom_boundary_collision_3p());

  event = Collision<>::CHROM_BOUNDARY();
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred());
  CHECK(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.is_chrom_boundary_collision());
  CHECK(Collision<>{idx, event}.is_chrom_boundary_collision_5p());
  CHECK_FALSE(Collision<>{idx, event}.is_chrom_boundary_collision_3p());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Collision encoding - LEF-BAR", "[simulation][short]") {
  const usize idx = 123;  // NOLINT

  CollisionEvent<u8f> event = Collision<>::COLLISION() | Collision<>::LEF_BAR();
  CHECK(Collision<>{idx, event}.collision_occurred());
  CHECK_FALSE(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.is_lef_bar_collision());

  event = Collision<>::LEF_BAR();
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred());
  CHECK(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.is_lef_bar_collision());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Collision encoding - LEF-LEF primary", "[simulation][short]") {
  const usize idx = 123;  // NOLINT

  CollisionEvent<u8f> event = Collision<>::COLLISION() | Collision<>::LEF_LEF_PRIMARY();
  CHECK(Collision<>{idx, event}.collision_occurred());
  CHECK_FALSE(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.is_lef_lef_primary_collision());

  event = Collision<>::LEF_LEF_PRIMARY();
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred());
  CHECK(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.is_lef_lef_primary_collision());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Collision encoding - LEF-LEF secondary", "[simulation][short]") {
  const usize idx = 123;  // NOLINT

  CollisionEvent<u8f> event = Collision<>::COLLISION() | Collision<>::LEF_LEF_SECONDARY();
  CHECK(Collision<>{idx, event}.collision_occurred());
  CHECK_FALSE(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.is_lef_lef_secondary_collision());

  event = Collision<>::LEF_LEF_SECONDARY();
  CHECK_FALSE(Collision<>{idx, event}.collision_occurred());
  CHECK(Collision<>{idx, event}.collision_avoided());
  CHECK(Collision<>{idx, event}.is_lef_lef_secondary_collision());
}
}  // namespace modle::test::libmodle
