// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "modle/simulation.hpp"

#include <absl/strings/str_join.h>              // for StrJoin
#include <absl/types/span.h>                    // for Span
#include <cpp-sort/sorters/counting_sorter.h>   // for counting_sort, counting_sorter
#include <cpp-sort/sorters/insertion_sorter.h>  // for insertion_sorter
#include <cpp-sort/sorters/ska_sorter.h>        // for ska_sort, ska_sorter
#include <cpp-sort/sorters/split_sorter.h>      // for split_sort, split_sorter
#include <fmt/compile.h>
#include <spdlog/spdlog.h>  // for warn

#include <algorithm>                         // for min
#include <boost/range/adaptor/reversed.hpp>  // for reversed_range, reverse
#include <cassert>                           // for assert
#include <cstddef>                           // for size_t
#include <thread>                            // for thread
#include <thread_pool/thread_pool.hpp>
#include <type_traits>  // for declval, decay_t

#include "modle/common/common.hpp"  // for random::PRNG_t
#include "modle/common/math.hpp"
#include "modle/common/random_sampling.hpp"             // for random_sampe
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNING_...
#include "modle/common/utils.hpp"                       // for ndebug_defined
#include "modle/extrusion_factors.hpp"                  // for Lef
#include "modle/genome.hpp"                             // for Chromosome

namespace modle {

template <typename MaskT>
void Simulation::bind_lefs(const bp_t start_pos, const bp_t end_pos, const absl::Span<Lef> lefs,
                           const absl::Span<size_t> rev_lef_ranks,
                           const absl::Span<size_t> fwd_lef_ranks, const MaskT& mask,
                           random::PRNG_t& rand_eng,
                           size_t current_epoch) noexcept(utils::ndebug_defined()) {
  using T = std::decay_t<decltype(std::declval<MaskT&>().operator[](std::declval<size_t>()))>;
  static_assert(std::is_integral_v<T> || std::is_same_v<MaskT, boost::dynamic_bitset<>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");
  {
    assert(lefs.size() <= mask.size() || mask.empty());             // NOLINT
    assert(std::all_of(rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
                       [&](const auto i) { return i < lefs.size(); }));
    assert(std::all_of(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),  // NOLINT
                       [&](const auto i) { return i < lefs.size(); }));
  }

  chrom_pos_generator_t pos_generator{start_pos, end_pos - 1};
  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (mask.empty() || mask[i]) {  // Bind all LEFs when mask is empty
      lefs[i].bind_at_pos(current_epoch, pos_generator(rand_eng));
    }
  }

  if constexpr (utils::ndebug_not_defined()) {
    for (auto i = 0UL; i < lefs.size(); ++i) {
      if (mask.empty() || mask[i]) {
        assert(lefs[i].rev_unit >= start_pos && lefs[i].rev_unit < end_pos);  // NOLINT
        assert(lefs[i].fwd_unit >= start_pos && lefs[i].fwd_unit < end_pos);  // NOLINT
      }
    }
  }

  Simulation::rank_lefs(lefs, rev_lef_ranks, fwd_lef_ranks, current_epoch != 0);
  {
    using IT = std::decay_t<decltype(rev_lef_ranks.front())>;
    (void)static_cast<IT*>(nullptr);
    assert(std::all_of(
        rev_lef_ranks.begin(), rev_lef_ranks.end(),  // NOLINT
        [&](const auto i) { return i < lefs.size() || i == (std::numeric_limits<IT>::max)(); }));
  }
}

template <typename MaskT>
void Simulation::bind_lefs(const Chromosome& chrom, const absl::Span<Lef> lefs,
                           const absl::Span<size_t> rev_lef_ranks,
                           const absl::Span<size_t> fwd_lef_ranks, const MaskT& mask,
                           random::PRNG_t& rand_eng,
                           size_t current_epoch) noexcept(utils::ndebug_defined()) {
  Simulation::bind_lefs(chrom.start_pos(), chrom.end_pos(), lefs, rev_lef_ranks, fwd_lef_ranks,
                        mask, rand_eng, current_epoch);
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
thread_pool Simulation::instantiate_thread_pool(I nthreads_, bool clamp_nthreads) {
  static_assert(std::is_integral_v<I>, "nthreads should have an integral type.");
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  if (clamp_nthreads) {
    return thread_pool(
        std::min(std::thread::hardware_concurrency(), static_cast<uint32_t>(nthreads_)));
  }
  assert(nthreads_ > 0);
  return thread_pool(static_cast<uint32_t>(nthreads_));
  DISABLE_WARNING_POP
}

}  // namespace modle

constexpr auto fmt::formatter<modle::Simulation::Task>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<modle::Simulation::Task>::format(const modle::Simulation::Task& t,
                                                     FormatContext& ctx) -> decltype(ctx.out()) {
  assert(t.chrom);  // NOLINT
  return fmt::format_to(ctx.out(), FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\t{}"), t.id, t.chrom->name(),
                        t.cell_id, t.num_target_epochs, t.num_target_contacts, t.num_lefs,
                        t.barriers.size());
}

constexpr auto fmt::formatter<modle::Simulation::TaskPW>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<modle::Simulation::TaskPW>::format(const modle::Simulation::TaskPW& t,
                                                       FormatContext& ctx) -> decltype(ctx.out()) {
  return fmt::format_to(
      ctx.out(), FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"), t.id,
      t.chrom ? t.chrom->name() : "null", t.cell_id, t.num_target_epochs, t.num_target_contacts,
      t.num_lefs, t.barriers.size(), t.deletion_begin, t.deletion_size, t.window_start,
      t.window_end, t.active_window_start, t.active_window_end, t.feats1.size(), t.feats2.size());
}

constexpr auto fmt::formatter<modle::Simulation::State>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<modle::Simulation::State>::format(const modle::Simulation::State& s,
                                                      FormatContext& ctx) -> decltype(ctx.out()) {
  return fmt::format_to(ctx.out(),
                        FMT_STRING("State:\n"
                                   " - TaskID: {:d}\n"
                                   " - CellID: {:d}\n"
                                   " - Chrom: {}:{:d}-{:d}\n"
                                   " - Deletion: {:d}-{:d}\n"
                                   " - Window: {:d}-{:d}\n"
                                   " - Active window: {:d}-{:d}\n"
                                   " - Current epoch: {:d}\n"
                                   " - Burn-in completed: {}\n"
                                   " - Target epochs: {:d}\n"
                                   " - Target contacts: {:d}\n"
                                   " - # of LEFs: {:d}\n"
                                   " - # of active LEFs: {:d}\n"
                                   " - # Extrusion barriers: {:d}\n"
                                   " - # Type I features: {:d}\n"
                                   " - # Type II features: {:d}\n"
                                   " - # of contacts registered: {:d}\n"
                                   " - seed: {:d}"),
                        s.id, s.cell_id, s.chrom ? s.chrom->name() : "null",
                        s.chrom ? static_cast<int64_t>(s.chrom->start_pos()) : -1,
                        s.chrom ? static_cast<int64_t>(s.chrom->end_pos()) : -1, s.deletion_begin,
                        s.deletion_begin + s.deletion_size, s.window_start, s.window_end,
                        s.active_window_start, s.active_window_end, s.epoch,
                        s.burnin_completed ? "True" : "False", s.num_target_epochs,
                        s.num_target_contacts, s.num_lefs, s.num_active_lefs, s.barriers.size(),
                        s.feats1.size(), s.feats2.size(), s.num_contacts, s.seed);
}
