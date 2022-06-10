// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>    // for assert
#include <stdexcept>  // for runtime_error
#include <string>     // for string

#include "modle/common/common.hpp"  // for u8f, MODLE_UNREACHABLE_CODE

namespace modle::dna {

/// Class to represent DNA direction or strandness
class Direction {
  u8f _direction;
  // 0 = none
  // 1 = rev
  // 2 = fwd
  // 3 = both

  explicit Direction(char d) {
    switch (d) {
      case '.':
        this->_direction = 0;
        break;
      case '-':
        this->_direction = 1;
        break;
      case '+':
        this->_direction = 2;
        break;
      default:
        throw std::runtime_error("Invalid dna::Direction '" + std::string{d} + "\"");
    }
  }

 public:
  constexpr Direction() = delete;
  template <u8f direction>
  constexpr Direction() : _direction(direction) {
    static_assert(direction <= 3);
  }
  explicit constexpr Direction(u8f direction) : _direction(direction) { assert(direction <= 3); }

  [[nodiscard]] constexpr Direction complement() const noexcept {
    switch (_direction) {
      case 0:
        [[fallthrough]];
      case 3:
        return Direction{this->_direction};
      case 1:
        return Direction{u8f(2)};
      case 2:
        return Direction{u8f(3)};
      default:
        MODLE_UNREACHABLE_CODE;
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
inline bool operator==(char a, Direction b) { return Direction(a) == b; }

inline bool operator==(Direction a, char b) { return b == a; }
constexpr bool operator!=(Direction a, Direction b) noexcept { return !(a == b); }
inline bool operator!=(char a, Direction b) { return !(a == b); }
inline bool operator!=(Direction a, char b) { return !(a == b); }

[[maybe_unused]] inline constexpr Direction NONE{u8f(0)};
[[maybe_unused]] inline constexpr Direction REV{u8f(1)};
[[maybe_unused]] inline constexpr Direction FWD{u8f(2)};
[[maybe_unused]] inline constexpr Direction BOTH{u8f(3)};

}  // namespace modle::dna
