#pragma once

// IWYU pragma: private, include "modle/simulation.hpp"

#include <absl/strings/str_join.h>              // for StrJoin
#include <absl/types/span.h>                    // for Span
#include <cpp-sort/sorters/counting_sorter.h>   // for counting_sort, counting_sorter
#include <cpp-sort/sorters/insertion_sorter.h>  // for insertion_sorter
#include <cpp-sort/sorters/ska_sorter.h>        // for ska_sort, ska_sorter
#include <cpp-sort/sorters/split_sorter.h>      // for split_sort, split_sorter
#include <fmt/compile.h>

#include <algorithm>                         // for min
#include <boost/range/adaptor/reversed.hpp>  // for reversed_range, reverse
#include <cassert>                           // for assert
#include <cstddef>                           // for size_t
#include <thread>                            // for thread
#include <thread_pool/thread_pool.hpp>
#include <type_traits>  // for declval, decay_t

#include "modle/common/common.hpp"                      // for random::PRNG_t
#include "modle/common/random_sampling.hpp"             // for random_sampe
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNING_...
#include "modle/common/utils.hpp"                       // for ndebug_defined
#include "modle/extrusion_factors.hpp"                  // for Lef
#include "modle/genome.hpp"                             // for Chromosome

namespace modle {

template <typename StateT>
inline void print_avg_loop_size(const size_t epoch, const size_t num_active_lefs,
                                const bool burnin_completed, const StateT& state) {
  static_assert(std::is_same<StateT, Simulation::State>() ||
                std::is_same<StateT, Simulation::StatePW>());
  std::vector<bp_t> loop_sizes(num_active_lefs, 0);
  std::transform(state.lef_buff.begin(),
                 state.lef_buff.begin() + static_cast<int32_t>(num_active_lefs), loop_sizes.begin(),
                 [](const auto& lef) { return lef.fwd_unit.pos() - lef.rev_unit.pos(); });
  // chrom_name, cell_id, epoch, burnin_completed, num_lefs, num_active_lefs, loop_sizes
  fmt::print(stdout, FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{}\t{}\n"), state.chrom->name(),
             state.cell_id, epoch, burnin_completed ? "True" : "False", state.num_lefs,
             num_active_lefs, fmt::join(loop_sizes, "\t"));
}

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
  return Simulation::bind_lefs(chrom.start_pos(), chrom.end_pos(), lefs, rev_lef_ranks,
                               fwd_lef_ranks, mask, rand_eng, current_epoch);
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
thread_pool Simulation::instantiate_thread_pool(I nthreads, bool clamp_nthreads) {
  static_assert(std::is_integral_v<I>, "nthreads should have an integral type.");
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  if (clamp_nthreads) {
    return thread_pool(
        std::min(std::thread::hardware_concurrency(), static_cast<uint32_t>(nthreads)));
  }
  assert(nthreads > 0);
  return thread_pool(static_cast<uint32_t>(nthreads));
  DISABLE_WARNING_POP
}

template <typename StateT, typename>
void Simulation::simulate_one_cell(StateT& s) const {
  constexpr auto normal_simulation = std::is_same_v<std::remove_reference_t<StateT>, State>;
  {
    assert(s.num_lefs == s.lef_buff.size());         // NOLINT
    assert(s.num_lefs == s.rank_buff1.size());       // NOLINT
    assert(s.num_lefs == s.rank_buff2.size());       // NOLINT
    assert(s.num_lefs == s.moves_buff1.size());      // NOLINT
    assert(s.num_lefs == s.moves_buff2.size());      // NOLINT
    assert(s.num_lefs == s.idx_buff.size());         // NOLINT
    assert(s.num_lefs == s.collision_buff1.size());  // NOLINT
    assert(s.num_lefs == s.collision_buff2.size());  // NOLINT
  }
  auto epoch = 0UL;  // The epoch needs to be declared outside of the try-catch body so that we can
  // use the current epoch when generating error messages
  try {
    // Seed is computed based on chrom. name, size and cellid
    s.seed = s.chrom->hash(s.xxh_state.get(), this->seed, s.cell_id);
    s.rand_eng = random::PRNG(s.seed);

    const auto lef_binding_rate_burnin =
        static_cast<double>(s.num_lefs) / static_cast<double>(this->burnin_lef_binding_epochs);

    // Compute the avg. # of LEFs to use to sample contact every iterations
    const auto avg_nlefs_to_sample =
        static_cast<double>(s.num_lefs) * this->lef_fraction_contact_sampling;

    // Compute the avg. # of LEFs to be released every iterations
    const auto avg_nlefs_to_release =
        static_cast<double>((this->rev_extrusion_speed + this->fwd_extrusion_speed) * s.num_lefs) /
        static_cast<double>(this->average_lef_lifetime);
    const auto avg_nlefs_to_release_burnin =
        static_cast<double>((this->rev_extrusion_speed_burnin + this->fwd_extrusion_speed_burnin) *
                            s.num_lefs) /
        static_cast<double>(this->average_lef_lifetime);

    // Local counter for the number of contacts generated by the current kernel instance
    auto n_contacts = 0UL;

    assert(avg_nlefs_to_sample <= static_cast<double>(s.num_lefs));   // NOLINT
    assert(avg_nlefs_to_release <= static_cast<double>(s.num_lefs));  // NOLINT
    assert((s.num_target_contacts != 0) ==                            // NOLINT
           (s.num_target_epochs == (std::numeric_limits<size_t>::max)()));

    // Declare the spans used by the burnin phase as well as the simulation itself
    auto lefs = absl::MakeSpan(s.lef_buff);
    auto lef_unloader_affinity = absl::MakeSpan(s.lef_unloader_affinity);
    const auto barriers = absl::MakeConstSpan(s.barriers);
    auto rev_lef_ranks = absl::MakeSpan(s.rank_buff1);
    auto fwd_lef_ranks = absl::MakeSpan(s.rank_buff2);
    auto rev_moves = absl::MakeSpan(s.moves_buff1);
    auto fwd_moves = absl::MakeSpan(s.moves_buff2);
    auto rev_collision_mask = absl::MakeSpan(s.collision_buff1);
    auto fwd_collision_mask = absl::MakeSpan(s.collision_buff2);

    // Generate initial extr. barrier states, so that they are already at or close to equilibrium
    for (auto i = 0UL; i < s.barrier_mask.size(); ++i) {
      s.barrier_mask[i] = random::bernoulli_trial{barriers[i].occupancy()}(s.rand_eng);
    }

    size_t nlefs_to_release;  // NOLINT
    size_t num_active_lefs = this->skip_burnin ? s.num_lefs : 0;
    auto target_burnin_epochs = this->burnin_epochs;
    auto burnin_completed = this->skip_burnin;
    // Start the burnin phase (followed by the actual simulation)
    for (; epoch < s.num_target_epochs; ++epoch) {
      print_avg_loop_size(epoch, num_active_lefs, burnin_completed, s);
      if (num_active_lefs != s.num_lefs) {
        s.num_target_epochs = std::max(s.num_target_epochs + 1, s.num_target_epochs);
        ++target_burnin_epochs;
        const auto num_lefs_to_bind =
            modle::random::poisson_distribution<size_t>{lef_binding_rate_burnin}(s.rand_eng);
        num_active_lefs = std::min(num_active_lefs + num_lefs_to_bind, s.num_lefs);
        if (num_active_lefs == 0) {  // Guard against the rare occasion where the poisson prng
                                     // samples 0 in the first epoch
          continue;
        }

        lefs = absl::MakeSpan(s.lef_buff.data(), num_active_lefs);
        lef_unloader_affinity = absl::MakeSpan(s.lef_unloader_affinity.data(), num_active_lefs);
        rev_lef_ranks = absl::MakeSpan(s.rank_buff1.data(), num_active_lefs);
        fwd_lef_ranks = absl::MakeSpan(s.rank_buff2.data(), num_active_lefs);
        rev_moves = absl::MakeSpan(s.moves_buff1.data(), num_active_lefs);
        fwd_moves = absl::MakeSpan(s.moves_buff2.data(), num_active_lefs);
        rev_collision_mask = absl::MakeSpan(s.collision_buff1.data(), num_active_lefs);
        fwd_collision_mask = absl::MakeSpan(s.collision_buff2.data(), num_active_lefs);

        // Sample nlefs to be released while in burn-in phase
        nlefs_to_release =
            std::min(lefs.size(), random::poisson_distribution<size_t>{
                                      static_cast<double>((this->rev_extrusion_speed_burnin +
                                                           this->fwd_extrusion_speed_burnin) *
                                                          num_active_lefs) /
                                      static_cast<double>(this->average_lef_lifetime)}(s.rand_eng));
      } else {  // Sample nlefs to be released after the burn-in phase has been completed
        burnin_completed = epoch > target_burnin_epochs;
        nlefs_to_release = std::min(
            lefs.size(),
            random::poisson_distribution<size_t>{
                burnin_completed ? avg_nlefs_to_release : avg_nlefs_to_release_burnin}(s.rand_eng));
      }

      ////////////////////////
      // Simulate one epoch //
      ////////////////////////

      {  // Select inactive LEFs and bind them
        auto lef_mask = absl::MakeSpan(s.idx_buff.data(), lefs.size());
        Simulation::select_lefs_to_bind(lefs, lef_mask);
        if constexpr (normal_simulation) {
          Simulation::bind_lefs(*s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, lef_mask, s.rand_eng,
                                epoch);
        } else {
          Simulation::bind_lefs(s.window_start, s.window_end, lefs, rev_lef_ranks, fwd_lef_ranks,
                                lef_mask, s.rand_eng, epoch);
        }
      }

      if (burnin_completed) {                        // Register contacts
        assert(fwd_lef_ranks.size() == s.num_lefs);  // NOLINT

        auto nlefs_to_sample = std::min(
            lefs.size(), random::poisson_distribution<size_t>{avg_nlefs_to_sample}(s.rand_eng));

        if (s.num_target_contacts != 0) {  // When using the target contact density as stopping
          // criterion, don't overshoot the target number of contacts
          nlefs_to_sample = std::min(nlefs_to_sample, s.num_target_contacts - n_contacts);
        }

        // Select LEFs to be used for contact registration.
        // We are sampling from fwd ranks to avoid having to allocate a vector of indices just to
        // do this sampling
        const auto lef_idx = absl::MakeSpan(s.idx_buff.data(), nlefs_to_sample);
        random_sample(fwd_lef_ranks.begin(), fwd_lef_ranks.end(), lef_idx.begin(), lef_idx.size(),
                      s.rand_eng);
        if (!std::all_of(lef_idx.begin(), lef_idx.end(),
                         [&](const auto i) { return i < lefs.size(); })) {
          throw std::runtime_error(
              fmt::format(FMT_STRING("lef_idx.size()={}; num_lefs={};\nlef_idx=[{}]\n"),
                          lef_idx.size(), lefs.size(), fmt::join(lef_idx, ", ")));
        }
        if constexpr (normal_simulation) {
          if (this->randomize_contacts) {
            n_contacts +=
                this->register_contacts_w_randomization(*s.chrom, lefs, lef_idx, s.rand_eng);
          } else {
            n_contacts += this->register_contacts(*s.chrom, lefs, lef_idx);
          }
        } else {
          if (this->randomize_contacts) {
            n_contacts += this->register_contacts_w_randomization(
                s.window_start + 1, s.window_end - 1, s.contacts, lefs, lef_idx, s.rand_eng);
          } else {
            n_contacts += this->register_contacts(s.window_start + 1, s.window_end - 1, s.contacts,
                                                  lefs, lef_idx);
          }
        }

        if (s.num_target_contacts != 0 && n_contacts >= s.num_target_contacts) {
          return;  // Enough contact have been generated. Yay!
        }
      }

      this->generate_moves(*s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                           burnin_completed, s.rand_eng);

      CTCF::update_states(barriers, s.barrier_mask, s.rand_eng);

      // Reset collision masks
      std::fill(rev_collision_mask.begin(), rev_collision_mask.end(), NO_COLLISION);
      std::fill(fwd_collision_mask.begin(), fwd_collision_mask.end(), NO_COLLISION);

      // Detect collision and correct moves
      const auto& [num_rev_units_at_5prime, num_fwd_units_at_3prime] =
          Simulation::process_collisions(*s.chrom, lefs, barriers, s.barrier_mask, rev_lef_ranks,
                                         fwd_lef_ranks, rev_moves, fwd_moves, rev_collision_mask,
                                         fwd_collision_mask, s.rand_eng);
      // Advance LEFs
      Simulation::extrude(*s.chrom, lefs, rev_moves, fwd_moves, num_rev_units_at_5prime,
                          num_fwd_units_at_3prime);

      // The vector of affinities is used to bias LEF release towards LEFs that are not in a hard
      // stall condition
      this->generate_lef_unloader_affinities(lefs, barriers, rev_collision_mask, fwd_collision_mask,
                                             lef_unloader_affinity);

      // Reusing this buffer is ok, as at this point we don't need access to collision information
      const auto lef_idx = absl::MakeSpan(s.idx_buff.data(), nlefs_to_release);
      Simulation::select_lefs_to_release(lef_idx, lef_unloader_affinity, s.rand_eng);
      Simulation::release_lefs(lefs, lef_idx);
    }
  } catch (const std::exception& err) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Caught an exception while simulating epoch #{} in cell #{} of {}:\n{}"), epoch,
        s.cell_id, s.chrom->name(), err.what()));
  }
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
                        FMT_STRING("State:\n - TaskID {}\n"
                                   " - Chrom: {}[{}-{}]\n"
                                   " - CellID: {}\n"
                                   " - Target epochs: {}\n"
                                   " - Target contacts: {}\n"
                                   " - # of LEFs: {}\n"
                                   " - # Extrusion barriers: {}\n"
                                   " - seed: {}\n"),
                        s.id, s.chrom ? s.chrom->name() : "null",
                        s.chrom ? static_cast<int64_t>(s.chrom->start_pos()) : -1,
                        s.chrom ? static_cast<int64_t>(s.chrom->end_pos()) : -1, s.cell_id,
                        s.num_target_epochs, s.num_target_contacts, s.num_lefs, s.barriers.size(),
                        s.seed);
}
