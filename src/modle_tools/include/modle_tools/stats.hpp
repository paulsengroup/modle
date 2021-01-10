#pragma once

#include <absl/types/span.h>

#include <cstdint>
#include <vector>

#include "modle/contacts.hpp"

namespace modle::tools {
template <typename I>
[[nodiscard]] inline std::vector<uint64_t> compute_expected_contacts(
    const ContactMatrix<I>& cmatrix);

template <typename I>
inline void compute_expected_contacts(const ContactMatrix<I>& cmatrix, std::vector<uint64_t>& buff);

template <typename I>
[[nodiscard]] inline ContactMatrix<I> subtract_expected_contacts(const ContactMatrix<I>& cmatrix1,
                                                                 absl::Span<const uint64_t> hist);
}  // namespace modle::tools

#include "../../stats_impl.hpp"
