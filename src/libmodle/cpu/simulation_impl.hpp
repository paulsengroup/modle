#pragma once

#include <absl/types/span.h>  // for Span, MakeSpan, MakeConstSpan

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <type_traits>

#include "modle/common.hpp"
#include "modle/dna.hpp"

namespace modle {
template <typename MaskT>
void Simulation::bind_lefs(const Chromosome* chrom, const absl::Span<Lef> lefs,
                           const absl::Span<size_t> rev_lef_ranks,
                           const absl::Span<size_t> fwd_lef_ranks, const MaskT& mask,
                           modle::PRNG_t& rand_eng,
                           bool first_epoch) noexcept(utils::ndebug_defined()) {
  static_assert(
      std::is_same_v<boost::dynamic_bitset<>, MaskT> ||
          std::is_integral_v<
              std::decay_t<decltype(std::declval<MaskT&>().operator[](std::declval<size_t>()))>>,
      "mask should be a vector of integral numbers or a boost::dynamic_bitset.");
  {
    assert(lefs.size() <= mask.size() || mask.empty());             // NOLINT
    assert(chrom);                                                  // NOLINT
    assert(std::all_of(rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
                       [&](const auto i) { return i < lefs.size(); }));
    assert(std::all_of(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),  // NOLINT
                       [&](const auto i) { return i < lefs.size(); }));
  }

  chrom_pos_generator_t pos_generator{chrom->start_pos(), chrom->end_pos() - 1};
  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (mask.empty() || mask[i]) {  // Bind all LEFs when mask is empty
      lefs[i].bind_at_pos(pos_generator(rand_eng));
    }
  }

  {
    assert(std::all_of(lefs.begin(), lefs.end(), [&](const auto& lef) {  // NOLINT
      return lef.rev_unit >= chrom->start_pos() && lef.rev_unit < chrom->end_pos();
    }));
    assert(std::all_of(lefs.begin(), lefs.end(), [&](const auto& lef) {  // NOLINT
      return lef.fwd_unit >= chrom->start_pos() && lef.fwd_unit < chrom->end_pos();
    }));
  }

  Simulation::rank_lefs(lefs, rev_lef_ranks, fwd_lef_ranks, !first_epoch);
  {
    assert(std::all_of(lefs.begin(), lefs.end(), [&](const auto& lef) {  // NOLINT
      return lef.rev_unit >= chrom->start_pos() && lef.rev_unit < chrom->end_pos();
    }));
    assert(std::all_of(lefs.begin(), lefs.end(), [&](const auto& lef) {  // NOLINT
      return lef.fwd_unit >= chrom->start_pos() && lef.fwd_unit < chrom->end_pos();
    }));
    assert(std::all_of(rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
                       [&](const auto i) { return i < lefs.size(); }));
    assert(std::all_of(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),  // NOLINT
                       [&](const auto i) { return i < lefs.size(); }));
  }
}
}  // namespace modle
