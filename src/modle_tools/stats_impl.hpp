#pragma once

#include <absl/types/span.h>

#include <cstdint>
#include <vector>

#include "modle/contacts.hpp"

namespace modle::tools {

template <typename I>
void compute_expected_contacts(const ContactMatrix<I>& cmatrix, std::vector<uint64_t>& buff) {
  buff.resize(cmatrix.nrows());
  std::fill(buff.begin(), buff.end(), 0);

  for (auto i = 0UL; i < cmatrix.ncols(); ++i) {
    for (auto j = i; j < i + cmatrix.nrows() && j < cmatrix.ncols(); ++j) {
      buff[i] += cmatrix.get(j - i, j);
    }
  }
}

template <typename I>
std::vector<uint64_t> compute_expected_contacts(const ContactMatrix<I>& cmatrix) {
  std::vector<uint64_t> buff;
  compute_expected_contacts(cmatrix, buff);
  return buff;
}

template <typename I>
ContactMatrix<I> subtract_expected_contacts(const ContactMatrix<I>& cmatrix1,
                                            absl::Span<const uint64_t> hist) {
  ContactMatrix<I> cmatrix2(cmatrix1.nrows(), cmatrix1.ncols());
  const auto nbins = cmatrix1.generate_mask_for_empty_rows().count();
  std::vector<double> expected_counts(hist.size());
  for (auto i = 0UL; i < cmatrix1.ncols(); ++i) {
    for (auto j = i; j < i + cmatrix1.nrows() && j < cmatrix1.ncols(); ++j) {
      const auto n1 = cmatrix1.get(i, j);
      const auto n_expected =
          static_cast<double>(hist[j - i]) / static_cast<double>(nbins - (j - i));
      const auto n2 = n1 > n_expected ? static_cast<I>(std::round(n1 - n_expected)) : 0;
      cmatrix2.set(i, j, n2);
    }
  }
  return cmatrix2;
}

}  // namespace modle::tools
