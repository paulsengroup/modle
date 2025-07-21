// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "modle/simulation.hpp"

#include <absl/types/span.h>
#include <fmt/format.h>

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <cassert>
#include <limits>
#include <thread>
#include <type_traits>

#include "modle/common/common.hpp"
#include "modle/common/random.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"
#include "modle/common/utils.hpp"
#include "modle/extrusion_factors.hpp"
#include "modle/genome.hpp"

namespace modle {
constexpr const Config& Simulation::config() const noexcept { return this->_config; }
constexpr const Config& Simulation::c() const noexcept { return this->config(); }

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
void Simulation::bind_lefs(const GenomicInterval& interval, const absl::Span<Lef> lefs,
                           const absl::Span<usize> rev_lef_ranks,
                           const absl::Span<usize> fwd_lef_ranks, const MaskT& mask,
                           random::PRNG_t& rand_eng,
                           usize current_epoch) noexcept(utils::ndebug_defined()) {
  Simulation::bind_lefs(interval.start(), interval.end(), lefs, rev_lef_ranks, fwd_lef_ranks, mask,
                        rand_eng, current_epoch);
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

constexpr bool Simulation::run_lef_lef_collision_trial(random::PRNG_t& rand_eng) const noexcept {
  return c().probability_of_extrusion_unit_bypass == 0.0 ||
         random::bernoulli_trial{1.0 - c().probability_of_extrusion_unit_bypass}(rand_eng);
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
  assert(t.interval);
  return fmt::format_to(ctx.out(), "{}\t{}\t{}\t{}\t{}\t{}\t{}", t.id, t.interval->chrom(),
                        t.cell_id, t.num_target_epochs, t.num_target_contacts, t.num_lefs,
                        !!t.interval ? t.interval->barriers().size() : 0);
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
  assert(s.interval);
  // clang-format off
  return fmt::format_to(
      ctx.out(),
      "State:\n"
      " - TaskID: {:d}\n"
      " - CellID: {:d}\n"
      " - Interval: {}\n"
      " - Current epoch: {:d}\n"
      " - Burn-in completed: {}\n"
      " - Target epochs: {:d}\n"
      " - Target contacts: {:d}\n"
      " - # of LEFs: {:d}\n"
      " - # of active LEFs: {:d}\n"
      " - # Extrusion barriers: {:d}\n"
      " - # of contacts registered: {:d}\n"
      " - seed: {:d}",
      s.id,
      s.cell_id,
      *s.interval,
      s.epoch,
      s.burnin_completed ? "True" : "False",
      s.num_target_epochs,
      s.num_target_contacts,
      s.num_lefs,
      s.num_active_lefs,
      s.barriers.size(),
      s.num_contacts,
      s.seed);
  // clang-format on
}
