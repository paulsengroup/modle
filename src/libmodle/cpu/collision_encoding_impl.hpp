// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once
#include <absl/numeric/bits.h>  // for countl_zero, popcount
#include <fmt/format.h>         // for format

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
    const auto idx_ok = is_valid_index(idx);
    const auto event_ok = is_valid_event(event);
    if (!idx_ok || !event_ok) {
      const auto msg = [&]() {
        if (!idx_ok && !event_ok) {
          return "invalid collision index and event";
        }
        if (!idx_ok) {
          return "invalid collision index";
        }
        return "invalid collision event";
      }();
      throw std::runtime_error(fmt::format(
          FMT_STRING("Collision<I>::encode(idx={}, event={:08b}): {}"), idx, event(), msg));
    }
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

template <class I>
constexpr void Collision<I>::set_idx(const usize idx) noexcept {
  assert(is_valid_index(idx));  // NOLINT
  this->_collision |= idx & INDEX_MASK;
}

template <class I>
constexpr void Collision<I>::set_event(const CollisionEventT event) noexcept {
  assert(is_valid_event(event));  // NOLINT
  this->_collision |= (event() << INDEX_BITS);
  assert(is_valid_event(this->decode_event()));  // NOLINT
}

template <class I>
constexpr void Collision<I>::set(const usize idx, const CollisionEventT event) noexcept {
  this->set_idx(idx);
  this->set_event(event);
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
constexpr usize Collision<I>::decode_index() const noexcept {
  return this->_collision & INDEX_MASK;
}

template <class I>
constexpr auto Collision<I>::decode_event() const noexcept -> CollisionEventT {
  using T = decltype(COLLISION());
  return T(this->_collision >> INDEX_BITS);
}

template <class I>
constexpr bool Collision<I>::collision_occurred() const noexcept {
  return this->_collision & (I(COLLISION()) << INDEX_BITS);
}

template <class I>
constexpr bool Collision<I>::collision_avoided() const noexcept {
  return this->_collision != 0 && !this->collision_occurred();
}

template <class I>
constexpr bool Collision<I>::is_chrom_boundary_collision() const noexcept {
  return this->_collision & (I(CHROM_BOUNDARY()) << INDEX_BITS);
}

template <class I>
constexpr bool Collision<I>::is_chrom_boundary_collision_5p() const noexcept {
  return this->_collision & (I(CHROM_BOUNDARY()) << INDEX_BITS) && this->decode_index() == 5;
}

template <class I>
constexpr bool Collision<I>::is_chrom_boundary_collision_3p() const noexcept {
  return this->_collision & (I(CHROM_BOUNDARY()) << INDEX_BITS) && this->decode_index() == 3;
}

template <class I>
constexpr bool Collision<I>::is_lef_bar_collision() const noexcept {
  return this->_collision & (I(LEF_BAR()) << INDEX_BITS);
}

template <class I>
constexpr bool Collision<I>::is_lef_lef_primary_collision() const noexcept {
  return this->_collision & (I(LEF_LEF_PRIMARY()) << INDEX_BITS);
}

template <class I>
constexpr bool Collision<I>::is_lef_lef_secondary_collision() const noexcept {
  return this->_collision & (I(LEF_LEF_SECONDARY()) << INDEX_BITS);
}

}  // namespace modle
