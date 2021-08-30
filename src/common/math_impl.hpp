#pragma once

#include <cmath>

namespace modle::math {
template <class T, class>
constexpr T abs_diff(const T a, const T b) noexcept {
  return a > b ? a - b : b - a;
}

}  // namespace modle::math
