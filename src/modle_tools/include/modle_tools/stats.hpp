#pragma once

#include <absl/types/span.h>  // for Span

#include <cstddef>  // IWYU pragma: keep for size_t
#include <cstdint>  // for uint64_t

namespace modle {
template <typename I>
class ContactMatrix;
}

namespace modle::tools {
template <typename I>
[[nodiscard]] inline size_t compute_number_of_contacts_after_depletion(
    const ContactMatrix<I>& cmatrix, absl::Span<const uint64_t> hist, size_t effective_nbins,
    double depletion_multiplier);
}  // namespace modle::tools

#include "../../stats_impl.hpp"  // IWYU pragma: keep
