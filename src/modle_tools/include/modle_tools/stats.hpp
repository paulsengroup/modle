#pragma once

#include <absl/types/span.h>

#include <boost/dynamic_bitset.hpp>
#include <cstdint>
#include <vector>

#include "modle/contacts.hpp"

namespace modle::tools {
template <typename I>
[[nodiscard]] inline std::vector<uint64_t> compute_row_wise_contact_histogram(
    const ContactMatrix<I>& cmatrix);

template <typename I>
inline void compute_row_wise_contact_histogram(const ContactMatrix<I>& cmatrix,
                                               std::vector<uint64_t>& buff);

template <typename I>
[[nodiscard]] inline std::size_t compute_number_of_contacts_after_depletion(
    const ContactMatrix<I>& cmatrix, absl::Span<const uint64_t> hist, std::size_t effective_nbins,
    double depletion_multiplier);

template <typename I>
[[nodiscard]] inline ContactMatrix<I> compute_depl_cmatrix(const ContactMatrix<I>& cmatrix1,
                                                           absl::Span<const uint64_t> hist,
                                                           const boost::dynamic_bitset<>& mask,
                                                           double depletion_multiplier);
}  // namespace modle::tools

#include "../../stats_impl.hpp"