// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

#pragma once

#include <algorithm>  // IWYU pragma: keep for min
#include <iterator>   // for distance

#include "modle/common/random.hpp"  // for uniform_int_distribution

namespace modle {
// Adapted from libcxx
// https://github.com/llvm/llvm-project/blob/339c7340422425755cec4d28a6ff1d1f6ea4a528/libcxx/include/__algorithm/sample.h#L48
template <class PopulationIterator, class SampleIterator, class Distance,
          class UniformRandomNumberGenerator>
inline SampleIterator random_sample(PopulationIterator first, PopulationIterator last,
                                    SampleIterator output_iter, Distance n,
                                    UniformRandomNumberGenerator& g) {
  auto unsampled_sz = static_cast<Distance>(std::distance(first, last));
  for (n = std::min(n, unsampled_sz); n != 0; ++first) {
    auto r = random::uniform_int_distribution<Distance>(0, --unsampled_sz)(g);
    if (r < n) {
      *output_iter++ = *first;
      --n;
    }
  }
  return output_iter;
}
}  // namespace modle
