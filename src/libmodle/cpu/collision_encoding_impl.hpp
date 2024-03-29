// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once
#include <absl/numeric/bits.h>  // for countl_zero, popcount
#include <fmt/format.h>         // for format, formatter

#include <stdexcept>  // for runtime_error

#include "modle/common/utils.hpp"  // for ndebug_defined, ndebug_not_defined

namespace modle {
template <class I>
constexpr auto CollisionEvent<I>::operator|(CollisionEvent<I> other) const noexcept
    -> CollisionEvent<I> {
  return this->_event | other._event;
}

template <class I>
constexpr auto CollisionEvent<I>::operator|=(CollisionEvent<I> other) noexcept
    -> CollisionEvent<I>& {
  this->_event |= other._event;
  return *this;
}

template <class I>
constexpr auto CollisionEvent<I>::operator&(CollisionEvent<I> other) const noexcept
    -> CollisionEvent<I> {
  return this->_event & other._event;
}

template <class I>
constexpr auto CollisionEvent<I>::operator&=(CollisionEvent<I> other) noexcept
    -> CollisionEvent<I>& {
  this->_event &= other._event;
  return *this;
}

template <class I>
constexpr auto CollisionEvent<I>::operator^(CollisionEvent<I> other) const noexcept
    -> CollisionEvent<I> {
  return this->_event ^ other._event;
}

template <class I>
constexpr auto CollisionEvent<I>::operator^=(CollisionEvent<I> other) noexcept
    -> CollisionEvent<I>& {
  this->_event ^= other._event;
  return *this;
}

template <class I>
constexpr bool CollisionEvent<I>::operator==(CollisionEvent<I> other) const noexcept {
  return (*this)() == other();
}

template <class I>
constexpr bool CollisionEvent<I>::operator!=(CollisionEvent<I> other) const noexcept {
  return !((*this)() == other());
}

template <class I>
constexpr I CollisionEvent<I>::operator()() const noexcept {
  return this->_event;
}

template <class I>
constexpr Collision<I>::Collision(usize idx,
                                  CollisionEventT event) noexcept(utils::ndebug_defined())
    : _collision(encode(idx, event)) {}

template <class I>
constexpr I Collision<I>::encode(const usize idx,
                                 const CollisionEventT event) noexcept(utils::ndebug_defined()) {
  if constexpr (utils::ndebug_not_defined()) {
    if (event == 0 && idx != 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Collision<I>::encode(idx={}, event={:08b}): idx must be 0 when "
                                 "no collision has occurred/been avoided"),
                      idx, event()));
    }
    Collision<I>::validate(idx, event);
  }
  return I(idx) | (I(event()) << INDEX_BITS);
}

template <class I>
constexpr bool Collision<I>::is_valid_index(const usize idx) noexcept(utils::ndebug_defined()) {
  return usize(absl::countl_zero(idx)) >= RESERVED_EVENT_BITS;
}

template <class I>
constexpr bool Collision<I>::is_valid_event(const CollisionEventT event) noexcept(
    utils::ndebug_defined()) {
  // Make sure reserved but unused bits are all 0
  bool ok = usize(absl::countl_zero(event())) >= RESERVED_EVENT_BITS - EVENT_BITS;
  if ((event & COLLISION) != 0) {
    // In case of a collision we require exactly one of the bits recording the reason of collision
    // to be set to 1
    ok &= absl::popcount(event()) == 2;
  } else {
    ok &= absl::popcount(event()) <= 1;
  }
  return ok;
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNUSED_PARAMETER
template <class I>
constexpr void Collision<I>::validate(const usize idx, const CollisionEventT event) {
  if constexpr (utils::ndebug_not_defined()) {
    if (idx == 0 && event == 0) {
      return;
    }

    const auto error_msg = [&]() -> std::string_view {
      const auto idx_ok = is_valid_index(idx);
      const auto event_ok = is_valid_event(event);
      if (!idx_ok || !event_ok) {
        if (!idx_ok && !event_ok) {
          return "invalid collision index and event";
        }
        if (!idx_ok) {
          return "invalid collision index";
        }
        return "invalid collision event";
      }
      if ((event & CHROM_BOUNDARY) != 0 && (idx != 5 && idx != 3)) {
        // 5 == 5'-end, 3 == 3'-end
        return "Index should be either 5 or 3 when event=CHROM_BOUNDARY";
      }
      return "";
    }();

    if (!error_msg.empty()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Collision<I>::encode(idx={}, event={:08b}): {}"), idx, event(), error_msg));
    }
  }
}
DISABLE_WARNING_POP

template <class I>
constexpr void Collision<I>::validate() const {
  Collision<I>::validate(this->decode_index(), this->decode_event());
}

template <class I>
constexpr I Collision<I>::lshift_event(const CollisionEventT event) noexcept {
  return I(event()) << INDEX_BITS;
}

template <class I>
constexpr void Collision<I>::set_idx(const usize idx) noexcept {
  assert(is_valid_index(idx));
  this->_collision |= idx & INDEX_MASK;
}

template <class I>
constexpr void Collision<I>::add_event(const CollisionEventT event) noexcept {
  assert(is_valid_event(event));
  this->_collision |= Collision<I>::lshift_event(event);
  this->validate();
}

template <class I>
constexpr void Collision<I>::set_event(const CollisionEventT event) noexcept {
  assert(is_valid_event(event));
  this->_collision = Collision<I>::lshift_event(event) | (this->_collision & INDEX_MASK);
  assert(is_valid_event(this->decode_event()));
}

template <class I>
constexpr void Collision<I>::set(const usize idx, const CollisionEventT event) noexcept {
  this->clear();
  this->set_idx(idx);
  this->add_event(event);
  this->validate();
}

template <class I>
constexpr void Collision<I>::clear() noexcept {
  this->_collision = 0;
}

template <class I>
constexpr bool Collision<I>::operator==(const Collision other) const noexcept {
  return this->_collision == other._collision;
}

template <class I>
constexpr bool Collision<I>::operator!=(const Collision other) const noexcept {
  return !(*this == other);
}

template <class I>
constexpr I Collision<I>::operator()() const noexcept {
  return this->_collision;
}

template <class I>
constexpr I& Collision<I>::operator()() noexcept {
  return this->_collision;
}

template <class I>
constexpr usize Collision<I>::decode_index() const noexcept {
  return this->_collision & INDEX_MASK;
}

template <class I>
constexpr auto Collision<I>::decode_event() const noexcept -> CollisionEventT {
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_CONVERSION
  return CollisionEventT(this->_collision >> INDEX_BITS);
  DISABLE_WARNING_POP
}

template <class I>
constexpr bool Collision<I>::collision_occurred() const noexcept {
  return this->_collision & Collision<I>::lshift_event(COLLISION);
}

template <class I>
constexpr bool Collision<I>::collision_avoided() const noexcept {
  return !this->collision_occurred() && this->_collision != 0;
}

template <class I>
constexpr bool Collision<I>::collision_occurred(const CollisionEventT c) const noexcept {
  assert((c & COLLISION) == 0);
  return this->decode_event() == (c | COLLISION);
}

template <class I>
constexpr bool Collision<I>::collision_avoided(const CollisionEventT c) const noexcept {
  assert((c & COLLISION) == 0);
  return !this->collision_occurred() && this->decode_event() == c;
}

}  // namespace modle

constexpr auto fmt::formatter<modle::Collision<>>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  auto it = ctx.begin(), end = ctx.end();
  if (it != end && (*it == 's' || *it == 'l')) {
    presentation = *it++;
  }

  // Check if reached the end of the range:
  if (it != end && *it != '}') {
    throw fmt::format_error("invalid format");
  }
  return it;
}

template <typename FormatContext>
auto fmt::formatter<modle::Collision<>>::format(const modle::Collision<>& c,
                                                FormatContext& ctx) const -> decltype(ctx.out()) {
  if (presentation == 's') {
    return fmt::format_to(ctx.out(), FMT_STRING("{}"), c());
  }
  assert(presentation == 'l');
  return fmt::format_to(ctx.out(), FMT_STRING("idx={}; event={:08b};"), c.decode_index(),
                        c.decode_event()());
}
