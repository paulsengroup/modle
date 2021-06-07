#pragma once

// IWYU pragma: private, include "modle/simulation.hpp"

#include <absl/types/span.h>  // for Span

#include <algorithm>                   // for min
#include <boost/asio/thread_pool.hpp>  // for thread_pool
#include <cassert>                     // for assert
#include <cstddef>                     // for size_t
#include <thread>                      // for thread
#include <type_traits>                 // for declval, decay_t

#include "modle/common.hpp"                      // for PRNG_t
#include "modle/dna.hpp"                         // for Chromosome
#include "modle/extrusion_factors.hpp"           // for Lef
#include "modle/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNING_...
#include "modle/utils.hpp"                       // for ndebug_defined

namespace modle {
template <typename MaskT>
void Simulation::bind_lefs(const Chromosome* chrom, const absl::Span<Lef> lefs,
                           const absl::Span<size_t> rev_lef_ranks,
                           const absl::Span<size_t> fwd_lef_ranks, const MaskT& mask,
                           modle::PRNG_t& rand_eng,
                           bool first_epoch) noexcept(utils::ndebug_defined()) {
  using T = std::decay_t<decltype(std::declval<MaskT&>().operator[](std::declval<size_t>()))>;
  static_assert(std::is_integral_v<T> || std::is_same_v<MaskT, boost::dynamic_bitset<>>,
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
    for (auto i = 0UL; i < lefs.size(); ++i) {
      if (mask.empty() || mask[i]) {
        assert(lefs[i].rev_unit >= chrom->start_pos() &&
               lefs[i].rev_unit < chrom->end_pos());  // NOLINT
        assert(lefs[i].fwd_unit >= chrom->start_pos() &&
               lefs[i].fwd_unit < chrom->end_pos());  // NOLINT
      }
    }
  }

  Simulation::rank_lefs(lefs, rev_lef_ranks, fwd_lef_ranks, !first_epoch);
  {
    using IT = std::decay_t<decltype(rev_lef_ranks.front())>;
    (void)static_cast<IT*>(nullptr);
    assert(std::all_of(
        rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
        [&](const auto i) { return i < lefs.size() || i == std::numeric_limits<IT>::max(); }));
  }
}

template <typename MaskT>
void Simulation::select_lefs_to_bind(const absl::Span<const Lef> lefs,
                                     MaskT& mask) noexcept(utils::ndebug_defined()) {
  using T = std::decay_t<decltype(std::declval<MaskT&>().operator[](std::declval<size_t>()))>;
  static_assert(std::is_integral_v<T> || std::is_same_v<MaskT, boost::dynamic_bitset<>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");
  assert(lefs.size() == mask.size());  // NOLINT
  for (auto i = 0UL; i < lefs.size(); ++i) {
    mask[i] = !lefs[i].is_bound();
  }
}

template <typename I>
boost::asio::thread_pool Simulation::instantiate_thread_pool(I nthreads, bool clamp_nthreads) {
  static_assert(std::is_integral_v<I>, "nthreads should have an integral type.");
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  if (clamp_nthreads) {
    return boost::asio::thread_pool(
        std::min(std::thread::hardware_concurrency(), static_cast<unsigned int>(nthreads)));
  }
  assert(nthreads > 0);
  return boost::asio::thread_pool(static_cast<unsigned int>(nthreads));
  DISABLE_WARNING_POP
}
}  // namespace modle
