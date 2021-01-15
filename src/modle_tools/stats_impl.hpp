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
void compute_row_wise_contact_histogram(const ContactMatrix<I>& cmatrix,
                                        std::vector<uint64_t>& buff) {
  buff.resize(cmatrix.nrows());
  std::fill(buff.begin(), buff.end(), 0);

  for (auto i = 0UL; i < cmatrix.ncols(); ++i) {
    for (auto j = i; j < i + cmatrix.nrows() && j < cmatrix.ncols(); ++j) {
      // j - i corresponds to the distance from the diagonal
      buff[j - i] += cmatrix.get(j, i);
    }
  }
}

template <typename I>
std::vector<uint64_t> compute_row_wise_contact_histogram(const ContactMatrix<I>& cmatrix) {
  std::vector<uint64_t> buff;
  compute_row_wise_contact_histogram(cmatrix, buff);
  return buff;
}

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

template <typename I>
ContactMatrix<I> compute_depl_cmatrix(const ContactMatrix<I>& cmatrix1,
                                      absl::Span<const uint64_t> hist,
                                      const boost::dynamic_bitset<>& mask,
                                      double depletion_multiplier) {
  assert(hist.size() == cmatrix1.nrows());

  // This histogram contains the average contact number (instead of the total)
  std::vector<uint64_t> row_wise_avg_contacts(hist.size());
  std::transform(hist.begin(), hist.end(), row_wise_avg_contacts.begin(), [&](const auto n) {
    return static_cast<uint64_t>(std::round(depletion_multiplier * (static_cast<double>(n)) /
                                            static_cast<double>(mask.count())));
  });

  auto cmatrix2 = cmatrix1;
  for (auto i = 0UL; i < cmatrix2.ncols(); ++i) {
    for (auto j = i; j < i + cmatrix2.nrows() && j < cmatrix2.ncols(); ++j) {
      // If none of the bins are graylisted, remove the average contacts from the observed contacts.
      // If avg. contacts > obs. contacts, clamp the value to 0
      if (mask[i] && mask[j]) {
        const auto n = cmatrix2.get(i, j);
        cmatrix2.set(j, i, n > row_wise_avg_contacts[j - i] ? n - row_wise_avg_contacts[j - i] : 0);
      } else {
        assert(cmatrix2.get(i, j) == 0);
      }
    }
  }
  return cmatrix2;
}
}  // namespace modle::tools
