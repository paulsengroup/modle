#pragma once

#include <absl/types/span.h>

#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <cstdint>
#include <vector>

#include "modle/contacts.hpp"

namespace modle::tools {

template <typename I>
std::size_t compute_number_of_contacts_after_depletion(const ContactMatrix<I>& cmatrix,
                                                       absl::Span<const uint64_t> hist,
                                                       std::size_t effective_nbins,
                                                       double depletion_multiplier) {
  assert(hist.size() == cmatrix.nrows());
  // This histogram contains the average contact number (instead of the total)
  std::vector<uint64_t> row_wise_avg_contacts(hist.size());
  std::transform(hist.begin(), hist.end(), row_wise_avg_contacts.begin(), [&](const auto n) {
    return static_cast<uint64_t>(std::round((depletion_multiplier * static_cast<double>(n)) /
                                            static_cast<double>(effective_nbins)));
  });

  std::size_t depl_contacts{0};
  for (auto i = 0UL; i < cmatrix.ncols(); ++i) {
    for (auto j = i; j < i + cmatrix.nrows() && j < cmatrix.ncols(); ++j) {
      // j - i corresponds to the distance from the diagonal
      // Only process bins that more contacts than average
      if (const auto n = cmatrix.get(j, i); n > row_wise_avg_contacts[j - i]) {
        depl_contacts += n - row_wise_avg_contacts[j - i];
      }
    }
  }

  return depl_contacts;
}
}  // namespace modle::tools
