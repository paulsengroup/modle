// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "modle/simulation.hpp"

#include <absl/types/span.h>  // for Span
#include <fmt/format.h>       // for format_parse_context, format_error

#include <BS_thread_pool.hpp>  // for BS::thread_pool
#include <algorithm>           // for min
#include <cassert>             // for assert
#include <limits>              // for numeric_limits
#include <thread>              // for thread
#include <type_traits>         // for declval, decay_t

#include "modle/common/common.hpp"                      // for usize, bp_t, i64, u32
#include "modle/common/random.hpp"                      // for PRNG_t
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARN...
#include "modle/common/utils.hpp"                       // for ndebug_defined, ndebug_not_defined
#include "modle/extrusion_factors.hpp"                  // for Lef, ExtrusionUnit
#include "modle/genome.hpp"                             // for Chromosome

namespace modle {

template <typename MaskT>
void Simulation::bind_lefs(const bp_t start_pos, const bp_t end_pos, const absl::Span<Lef> lefs,
                           const absl::Span<usize> rev_lef_ranks,
                           const absl::Span<usize> fwd_lef_ranks, const MaskT& mask,
                           random::PRNG_t& rand_eng,
                           usize current_epoch) noexcept(utils::ndebug_defined()) {
  using T = std::decay_t<decltype(std::declval<MaskT&>().operator[](std::declval<usize>()))>;
  static_assert(std::is_integral_v<T> || std::is_same_v<MaskT, boost::dynamic_bitset<>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");
  {
    assert(lefs.size() <= mask.size() || mask.empty());
    assert(std::all_of(rev_lef_ranks.begin(), rev_lef_ranks.end(),
                       [&](const auto i) { return i < lefs.size(); }));
    assert(std::all_of(fwd_lef_ranks.begin(), fwd_lef_ranks.end(),
                       [&](const auto i) { return i < lefs.size(); }));
  }

  chrom_pos_generator_t pos_generator{start_pos, end_pos - 1};
  for (usize i = 0; i < lefs.size(); ++i) {
    if (mask.empty() || mask[i]) {  // Bind all LEFs when mask is empty
      lefs[i].bind_at_pos(current_epoch, pos_generator(rand_eng));
    }
  }

  if constexpr (utils::ndebug_not_defined()) {
    for (usize i = 0; i < lefs.size(); ++i) {
      if (mask.empty() || mask[i]) {
        assert(lefs[i].rev_unit >= start_pos && lefs[i].rev_unit < end_pos);
        assert(lefs[i].fwd_unit >= start_pos && lefs[i].fwd_unit < end_pos);
      }
    }
  }

  Simulation::rank_lefs(lefs, rev_lef_ranks, fwd_lef_ranks, current_epoch != 0);

  assert(std::all_of(rev_lef_ranks.begin(), rev_lef_ranks.end(), [&](const auto i) {
    using IT = std::decay_t<decltype(rev_lef_ranks.front())>;
    return i < lefs.size() || i == (std::numeric_limits<IT>::max)();
  }));
}

template <typename MaskT>
void Simulation::bind_lefs(const Chromosome& chrom, const absl::Span<Lef> lefs,
                           const absl::Span<usize> rev_lef_ranks,
                           const absl::Span<usize> fwd_lef_ranks, const MaskT& mask,
                           random::PRNG_t& rand_eng,
                           usize current_epoch) noexcept(utils::ndebug_defined()) {
  Simulation::bind_lefs(chrom.start_pos(), chrom.end_pos(), lefs, rev_lef_ranks, fwd_lef_ranks,
                        mask, rand_eng, current_epoch);
}

template <typename MaskT>
void Simulation::select_lefs_to_bind(const absl::Span<const Lef> lefs,
                                     MaskT& mask) noexcept(utils::ndebug_defined()) {
  using T = std::decay_t<decltype(std::declval<MaskT&>().operator[](std::declval<usize>()))>;
  static_assert(std::is_integral_v<T> || std::is_same_v<MaskT, boost::dynamic_bitset<>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");
  assert(lefs.size() == mask.size());
  for (usize i = 0; i < lefs.size(); ++i) {
    mask[i] = !lefs[i].is_bound();
  }
}

template <typename I>
BS::thread_pool Simulation::instantiate_thread_pool(I nthreads_, bool clamp_nthreads) {
  static_assert(std::is_integral_v<I>, "nthreads should have an integral type.");
  if (clamp_nthreads) {
    return BS::thread_pool(std::min(std::thread::hardware_concurrency(),
                                    utils::conditional_static_cast<u32>(nthreads_)));
  }
  assert(nthreads_ > 0);
  return BS::thread_pool(static_cast<u32>(nthreads_));
}

template <class TaskT>
usize Simulation::consume_tasks_blocking(moodycamel::BlockingConcurrentQueue<TaskT>& task_queue,
                                         moodycamel::ConsumerToken& ctok,
                                         absl::FixedArray<TaskT>& task_buff) {
  while (this->ok()) {
    const auto avail_tasks = task_queue.wait_dequeue_bulk_timed(
        ctok, task_buff.begin(), task_buff.size(), std::chrono::milliseconds(10));
    // Check whether dequeue operation timed-out before any task became available
    if (avail_tasks == 0) {
      // Reached end of simulation (i.e. all tasks have been processed)
      if (this->_end_of_simulation) {
        return 0;
      }
      // Keep waiting until one or more tasks become available
      continue;
    }
    return avail_tasks;
  }
  return 0;
}

constexpr bool Simulation::run_lef_lef_collision_trial(random::PRNG_t& rand_eng) const noexcept {
  return this->probability_of_extrusion_unit_bypass == 0.0 ||
         random::bernoulli_trial{1.0 - this->probability_of_extrusion_unit_bypass}(rand_eng);
}

constexpr bool Simulation::run_lef_bar_collision_trial(const double pblock,
                                                       random::PRNG_t& rand_eng) const noexcept {
  return pblock == 1.0 || random::bernoulli_trial{pblock}(rand_eng);
}

}  // namespace modle

constexpr auto fmt::formatter<modle::Simulation::Task>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.end();
}

template <typename FormatContext>
auto fmt::formatter<modle::Simulation::Task>::format(const modle::Simulation::Task& t,
                                                     FormatContext& ctx) const
    -> decltype(ctx.out()) {
  assert(t.chrom);
  return fmt::format_to(ctx.out(), FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\t{}"), t.id, t.chrom->name(),
                        t.cell_id, t.num_target_epochs, t.num_target_contacts, t.num_lefs,
                        t.barriers.size());
}

constexpr auto fmt::formatter<modle::Simulation::TaskPW>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.end();
}

template <typename FormatContext>
auto fmt::formatter<modle::Simulation::TaskPW>::format(const modle::Simulation::TaskPW& t,
                                                       FormatContext& ctx) const
    -> decltype(ctx.out()) {
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
  return ctx.end();
}

template <typename FormatContext>
auto fmt::formatter<modle::Simulation::State>::format(const modle::Simulation::State& s,
                                                      FormatContext& ctx) const
    -> decltype(ctx.out()) {
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
                        s.chrom ? static_cast<modle::i64>(s.chrom->start_pos()) : -1,
                        s.chrom ? static_cast<modle::i64>(s.chrom->end_pos()) : -1,
                        s.deletion_begin, s.deletion_begin + s.deletion_size, s.window_start,
                        s.window_end, s.active_window_start, s.active_window_end, s.epoch,
                        s.burnin_completed ? "True" : "False", s.num_target_epochs,
                        s.num_target_contacts, s.num_lefs, s.num_active_lefs, s.barriers.size(),
                        s.feats1.size(), s.feats2.size(), s.num_contacts, s.seed);
}
