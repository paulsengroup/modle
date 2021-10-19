// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>  // for declval, enable_if_t

#include "modle/common/utils.hpp"  // for identity

namespace modle::stats {
enum binomial_test_alternative : std::uint_fast8_t { TWO_SIDED, GREATER, LESS };

template <binomial_test_alternative alternative = TWO_SIDED, class FP = double, class I,
          class = std::enable_if<std::is_integral_v<I> && std::is_floating_point_v<FP>>>
[[nodiscard]] inline FP binomial_test(I k, I n, FP p = 0.5);

}  // namespace modle::stats
#include "../../../test_impl.hpp"  // IWYU pragma: export
