//===----------------------------------------------------------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
#pragma once

#include <algorithm>
#include <boost/random/uniform_int_distribution.hpp>
#include <iterator>

namespace modle {
template <class PopulationIterator, class SampleIterator, class Distance,
          class UniformRandomNumberGenerator>
inline SampleIterator random_sample(PopulationIterator first, PopulationIterator last,
                                    SampleIterator output_iter, Distance n,
                                    UniformRandomNumberGenerator& g) {
  auto unsampled_sz = static_cast<Distance>(std::distance(first, last));
  for (n = std::min(n, unsampled_sz); n != 0; ++first) {
    auto r = boost::random::uniform_int_distribution<Distance>(0, --unsampled_sz)(g);
    if (r < n) {
      *output_iter++ = *first;
      --n;
    }
  }
  return output_iter;
}
}  // namespace modle
