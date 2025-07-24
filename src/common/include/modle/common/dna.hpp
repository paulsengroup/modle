// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <stdexcept>
#include <string>

#include "modle/common/common.hpp"

namespace modle::dna {

/// Class to represent DNA direction or strandness
class Direction {
  std::uint_fast8_t _direction;
  // 0 = none
  // 1 = rev
  // 2 = fwd
  // 3 = both

 public:
  constexpr Direction() = delete;
  template <std::uint_fast8_t direction>
  constexpr Direction() : _direction(direction) {
    static_assert(direction <= 3);
  }
  explicit constexpr Direction(std::uint_fast8_t direction) : _direction(direction) {
    assert(direction <= 3);
  }
  [[nodiscard]] static Direction from_char(char d) {
    switch (d) {
      case '.':
        return Direction(0);
      case '-':
        return Direction(1);
      case '+':
        return Direction(2);
      default:
        throw std::runtime_error("Invalid dna::Direction '" + std::string{d} + "\"");
    }
  }

  [[nodiscard]] constexpr Direction complement() const noexcept {
    switch (_direction) {
      case 0:
        [[fallthrough]];
      case 3:
        return Direction{_direction};
      case 1:
        return Direction{2};
      case 2:
        return Direction{1};
      default:
        utils::unreachable_code();
    }
  }
  friend constexpr bool operator==(Direction a, Direction b) noexcept;
  friend bool operator==(char a, Direction b);
  friend bool operator==(Direction a, char b);
  friend constexpr bool operator!=(Direction a, Direction b) noexcept;
  friend bool operator!=(char a, Direction b);
  friend bool operator!=(Direction a, char b);
};

constexpr bool operator==(Direction a, Direction b) noexcept {
  return a._direction == b._direction;
}
inline bool operator==(char a, Direction b) { return Direction::from_char(a) == b; }

inline bool operator==(Direction a, char b) { return b == a; }
constexpr bool operator!=(Direction a, Direction b) noexcept { return !(a == b); }
inline bool operator!=(char a, Direction b) { return !(a == b); }
inline bool operator!=(Direction a, char b) { return !(a == b); }

[[maybe_unused]] inline constexpr Direction NONE{0};
[[maybe_unused]] inline constexpr Direction REV{1};
[[maybe_unused]] inline constexpr Direction FWD{2};
[[maybe_unused]] inline constexpr Direction BOTH{3};

}  // namespace modle::dna
