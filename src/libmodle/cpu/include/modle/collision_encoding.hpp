// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>  // for formatter

#include <bitset>  // for bitset
#include <type_traits>

#include "modle/common/common.hpp"  // for u8f, u32f
#include "modle/common/utils.hpp"   // for ndebug_defined, ndebug_not_defined

namespace modle {

template <class I>
class Collision;

/// This class is basically an enum that implements the bitwise operations that we need to represent
/// all kinds of collision events (e.g. an actuall LEF-BAR collision or the avoidance of a primary
/// LEF-LEF collision event)
template <class I = u8f>
class CollisionEvent {
  static_assert(std::is_unsigned_v<I>);

  template <class II>
  friend class Collision;

  I _event{0};

  static constexpr usize EVENT_BITS = 5;
  static constexpr usize RESERVED_EVENT_BITS = 8;

 public:
  // clang and gcc7 complain if the constructors are not explicitly defined here
  constexpr CollisionEvent() noexcept : _event(0) {}
  // NOLINTNEXTLINE(hicpp-explicit-conversions)
  constexpr CollisionEvent(I event) noexcept : _event(event) {}

  constexpr CollisionEvent operator|(CollisionEvent other) const noexcept;
  constexpr CollisionEvent& operator|=(CollisionEvent other) noexcept;
  constexpr CollisionEvent operator&(CollisionEvent other) const noexcept;
  constexpr CollisionEvent& operator&=(CollisionEvent other) noexcept;
  constexpr CollisionEvent operator^(CollisionEvent other) const noexcept;
  constexpr CollisionEvent& operator^=(CollisionEvent other) noexcept;

  constexpr bool operator==(CollisionEvent other) const noexcept;
  constexpr bool operator!=(CollisionEvent other) const noexcept;

  [[nodiscard]] constexpr I operator()() const noexcept;
};

template <class I = u32f>
class Collision {
 private:
  I _collision{0};

  using CollisionEventT = CollisionEvent<u8f>;

  static constexpr usize EVENT_BITS = CollisionEventT::EVENT_BITS;
  static constexpr usize RESERVED_EVENT_BITS = CollisionEventT::RESERVED_EVENT_BITS;
  static constexpr usize INDEX_BITS = (sizeof(I) * 8) - RESERVED_EVENT_BITS;
  static constexpr I EVENT_MASK = (~I(0)) << (INDEX_BITS - 1);
  static constexpr I INDEX_MASK = ~EVENT_MASK;

 public:
  constexpr Collision() noexcept : _collision(0) {}
  constexpr Collision(usize idx, CollisionEventT event) noexcept(utils::ndebug_defined());

  constexpr void set_idx(usize idx) noexcept;
  constexpr void add_event(CollisionEventT event) noexcept;
  constexpr void set_event(CollisionEventT event) noexcept;
  constexpr void set(usize idx, CollisionEventT event) noexcept;
  constexpr void clear() noexcept;

  constexpr bool operator==(Collision other) const noexcept;
  constexpr bool operator!=(Collision other) const noexcept;
  [[nodiscard]] constexpr I operator()() const noexcept;
  [[nodiscard]] constexpr I& operator()() noexcept;

  [[nodiscard]] constexpr usize decode_index() const noexcept;
  [[nodiscard]] constexpr CollisionEventT decode_event() const noexcept;

  [[nodiscard]] constexpr bool collision_occurred() const noexcept;
  [[nodiscard]] constexpr bool collision_avoided() const noexcept;
  [[nodiscard]] constexpr bool collision_occurred(CollisionEventT c) const noexcept;
  [[nodiscard]] constexpr bool collision_avoided(CollisionEventT c) const noexcept;

  // clang-format off
  static constexpr CollisionEventT NO_COLLISION =      0b0000'0000;
  static constexpr CollisionEventT COLLISION =         0b0001'0000;
  static constexpr CollisionEventT CHROM_BOUNDARY =    0b0000'1000;
  static constexpr CollisionEventT LEF_BAR =           0b0000'0100;
  static constexpr CollisionEventT LEF_LEF_PRIMARY =   0b0000'0010;
  static constexpr CollisionEventT LEF_LEF_SECONDARY = 0b0000'0001;
  // clang-format on

 private:
  [[nodiscard]] static constexpr I encode(usize idx,
                                          CollisionEventT event) noexcept(utils::ndebug_defined());
  [[nodiscard]] static constexpr bool is_valid_index(usize idx) noexcept(utils::ndebug_defined());
  [[nodiscard]] static constexpr bool is_valid_event(CollisionEventT event) noexcept(
      utils::ndebug_defined());
  constexpr void validate() const;
  static constexpr void validate(usize idx, CollisionEventT event);
  [[nodiscard]] static constexpr I lshift_event(CollisionEventT event) noexcept;
};

}  // namespace modle

template <>
struct fmt::formatter<modle::Collision<>> {
  char presentation = 's';  // s == short, l == long
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  inline auto format(const modle::Collision<>& c, FormatContext& ctx) -> decltype(ctx.out());
};

#include "../../collision_encoding_impl.hpp"
