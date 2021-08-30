#pragma once

#include <type_traits>

namespace modle::math {

template <class T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr T abs_diff(T a, T b) noexcept;

}  // namespace modle::math
#include "../../../math_impl.hpp"  // IWYU pragma: export
