// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>      // for uint_fast8_t
#include <type_traits>  // for enable_if_t

#include "modle/common/utils.hpp"  // for identity

namespace modle::stats {
enum binomial_test_alternative : u8f { TWO_SIDED, GREATER, LESS };

template <binomial_test_alternative alternative = TWO_SIDED, class FP = double, class I,
          class = std::enable_if<std::is_integral_v<I> && std::is_floating_point_v<FP>>>
[[nodiscard]] inline FP binomial_test(I x_, I x0_, FP x1_ = 0.5);

}  // namespace modle::stats
#include "../../../tests_impl.hpp"  // IWYU pragma: export
