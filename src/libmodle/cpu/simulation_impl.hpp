#pragma once

// IWYU pragma: private, include "modle/simulation.hpp"

#include <absl/strings/str_join.h>              // for StrJoin
#include <absl/types/span.h>                    // for Span
#include <cpp-sort/sorters/counting_sorter.h>   // for counting_sort, counting_sorter
#include <cpp-sort/sorters/insertion_sorter.h>  // for insertion_sorter
#include <cpp-sort/sorters/ska_sorter.h>        // for ska_sort, ska_sorter
#include <cpp-sort/sorters/split_sorter.h>      // for split_sort, split_sorter

#include <algorithm>                         // for min
#include <boost/asio/thread_pool.hpp>        // for thread_pool
#include <boost/range/adaptor/reversed.hpp>  // for reversed_range, reverse
#include <cassert>                           // for assert
#include <cstddef>                           // for size_t
#include <thread>                            // for thread
#include <type_traits>                       // for declval, decay_t

#include "modle/common/common.hpp"                      // for random::PRNG_t
#include "modle/common/random_sampling.hpp"             // for random_sampe
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNING_...
#include "modle/common/utils.hpp"                       // for ndebug_defined
#include "modle/extrusion_factors.hpp"                  // for Lef
#include "modle/genome.hpp"                             // for Chromosome

namespace modle {
template <typename MaskT>
void Simulation::bind_lefs(const Chromosome& chrom, const absl::Span<Lef> lefs,
                           const absl::Span<size_t> rev_lef_ranks,
                           const absl::Span<size_t> fwd_lef_ranks, const MaskT& mask,
                           random::PRNG_t& rand_eng, size_t current_epoch,
                           const bp_t deletion_begin,
                           const bp_t deletion_size) noexcept(utils::ndebug_defined()) {
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

  chrom_pos_generator_t pos_generator{chrom.start_pos(), chrom.end_pos() - 1};
  for (auto i = 0UL; i < lefs.size(); ++i) {
    if (mask.empty() || mask[i]) {  // Bind all LEFs when mask is empty
      auto pos = pos_generator(rand_eng);
      if (deletion_size > 0) {
        const auto deletion_end = deletion_begin + deletion_size;
        while (pos >= deletion_begin && pos < deletion_end) {
          pos = pos_generator(rand_eng);
        }
      }
      lefs[i].bind_at_pos(current_epoch, pos);
    }
  }

  {
    for (auto i = 0UL; i < lefs.size(); ++i) {
      if (mask.empty() || mask[i]) {
        assert(lefs[i].rev_unit >= chrom.start_pos() &&
               lefs[i].rev_unit < chrom.end_pos());  // NOLINT
        assert(lefs[i].fwd_unit >= chrom.start_pos() &&
               lefs[i].fwd_unit < chrom.end_pos());  // NOLINT
      }
    }
  }

  Simulation::rank_lefs(lefs, rev_lef_ranks, fwd_lef_ranks, current_epoch != 0);
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
    assert(s.num_lefs == s.epoch_buff.size());       // NOLINT
  }
  auto epoch = 0UL;  // The epoch needs to be declared outside of the try-catch body so that we can
  // use the current epoch when generating error messages
  try {
    // Seed is computed based on chrom. name, size and cellid
    s.seed = s.chrom->hash(s.xxh_state.get(), this->seed, s.cell_id);
    s.rand_eng = random::PRNG(s.seed);

    // Generate the epoch at which each LEF is supposed to be initially loaded
    auto lef_initial_loading_epoch = absl::MakeSpan(s.epoch_buff);
    // lef_initial_loading_epoch.resize(this->skip_burnin ? 0 : s.nlefs);

    if (!this->skip_burnin) {
      // TODO Consider using a Poisson process instead of sampling from an uniform distribution
      random::uniform_int_distribution<size_t> round_gen{
          0, (4 * this->average_lef_lifetime) / this->bin_size};
      std::generate(lef_initial_loading_epoch.begin(), lef_initial_loading_epoch.end(),
                    [&]() { return round_gen(s.rand_eng); });

      // Sort epochs in descending order
      if (round_gen.max() > 2048) {
        // Counting sort uses n + r space in memory, where r is the number of unique values in the
        // range to be sorted. For this reason it is not a good idea to use it when the sampling
        // interval is relatively large. Whether 2048 is a reasonable threshold has yet to be tested
        cppsort::ska_sort(lef_initial_loading_epoch.rbegin(), lef_initial_loading_epoch.rend());
      } else {
        cppsort::counting_sort(lef_initial_loading_epoch.rbegin(),
                               lef_initial_loading_epoch.rend());
      }
    }

    // Shift epochs so that the first epoch == 0
    if (const auto offset = lef_initial_loading_epoch.back(); offset != 0) {
      std::for_each(lef_initial_loading_epoch.begin(), lef_initial_loading_epoch.end(),
                    [&](auto& n) { n -= offset; });
    }

    // The highest loading epoch equals to the number of burnin epochs
    const auto n_burnin_epochs = this->skip_burnin ? 0 : lef_initial_loading_epoch.front();
    s.num_target_epochs += s.num_target_contacts != 0 ? 0 : n_burnin_epochs;

    // Compute the avg. # of LEFs to use to sample contact every iterations
    const auto avg_nlefs_to_sample =
        static_cast<double>(s.num_lefs) * this->lef_fraction_contact_sampling;

    // Compute the avg. # of LEFs to be released every iterations
    const auto avg_nlefs_to_release =
        static_cast<double>((this->rev_extrusion_speed + this->fwd_extrusion_speed) * s.num_lefs) /
        static_cast<double>(this->average_lef_lifetime);

    // Local counter for the number of contacts generated by the current kernel instance
    auto n_contacts = 0UL;

    assert(avg_nlefs_to_sample <= static_cast<double>(s.num_lefs));   // NOLINT
    assert(avg_nlefs_to_release <= static_cast<double>(s.num_lefs));  // NOLINT
    assert((s.num_target_contacts != 0) ==                            // NOLINT
           (s.num_target_epochs == std::numeric_limits<size_t>::max()));

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
    // Start the burnin phase (followed by the actual simulation)
    for (; epoch < s.num_target_epochs; ++epoch) {
      if (!lef_initial_loading_epoch.empty()) {  // Execute this branch only during the burnin phase
        if (this->skip_burnin) {
          // The following will cause every LEF to be loaded in the current iteration
          lef_initial_loading_epoch.subspan(0, 0);
        }

        // Consume epochs for LEFs that are supposed to be loaded in the current epoch
        auto nlefs_to_bind = 0UL;
        for (const auto n : boost::adaptors::reverse(lef_initial_loading_epoch)) {
          if (n != epoch) {
            break;
          }
          ++nlefs_to_bind;
        }
        lef_initial_loading_epoch.remove_suffix(nlefs_to_bind);
        // Don't remove this assertion! It is very useful to prevent empty Spans and other
        // issues that can cause weird issues later on during the simulation NOLINTNEXTLINE
        assert(lef_initial_loading_epoch.size() < s.num_lefs);  // NOLINT

        // Compute the current number of active LEFs (i.e. LEFs that are either bound to DNA, or
        // that are valid candidates for re-binding). Grow spans accordingly
        const auto nlefs = s.num_lefs - lef_initial_loading_epoch.size();

        lefs = absl::MakeSpan(s.lef_buff.data(), nlefs);
        lef_unloader_affinity = absl::MakeSpan(s.lef_unloader_affinity.data(), nlefs);
        rev_lef_ranks = absl::MakeSpan(s.rank_buff1.data(), nlefs);
        fwd_lef_ranks = absl::MakeSpan(s.rank_buff2.data(), nlefs);
        rev_moves = absl::MakeSpan(s.moves_buff1.data(), nlefs);
        fwd_moves = absl::MakeSpan(s.moves_buff2.data(), nlefs);
        rev_collision_mask = absl::MakeSpan(s.collision_buff1.data(), nlefs);
        fwd_collision_mask = absl::MakeSpan(s.collision_buff2.data(), nlefs);

        // Sample nlefs to be released while in burn-in phase
        nlefs_to_release =
            std::min(lefs.size(),
                     random::poisson_distribution<size_t>{
                         static_cast<double>(
                             (this->rev_extrusion_speed + this->fwd_extrusion_speed) * s.num_lefs) /
                         static_cast<double>(this->average_lef_lifetime)}(s.rand_eng));
      } else {  // Sample nlefs to be released after the burn-in phase has been completed
        nlefs_to_release = std::min(
            lefs.size(), random::poisson_distribution<size_t>{avg_nlefs_to_release}(s.rand_eng));
      }

      ////////////////////////
      // Simulate one epoch //
      ////////////////////////

      {  // Select inactive LEFs and bind them
        auto lef_mask = absl::MakeSpan(s.idx_buff.data(), lefs.size());
        Simulation::select_lefs_to_bind(lefs, lef_mask);
        // if constexpr (normal_simulation) {
        Simulation::bind_lefs(*s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, lef_mask, s.rand_eng,
                              epoch);
        // } else {
        //  Simulation::bind_lefs(*s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, lef_mask,
        //  s.rand_eng, epoch, s.deletion_begin, s.deletion_size);
        //}
      }

      if (epoch > n_burnin_epochs) {                 // Register contacts
        assert(fwd_lef_ranks.size() == s.num_lefs);  // NOLINT

        auto nlefs_to_sample = std::min(
            lefs.size(), random::poisson_distribution<size_t>{avg_nlefs_to_sample}(s.rand_eng));

        if (s.num_target_contacts != 0) {  // When using the target contact density as stopping
          // criterion, don't overshoot the target number of contacts
          nlefs_to_sample = std::min(nlefs_to_sample, s.num_target_contacts - n_contacts);
        }

        // Select LEFs to be used for contact registration.
        // We are sampling from fwd ranks to avoid having to allocate a vector of indices just to do
        // this sampling
        const auto lef_idx = absl::MakeSpan(s.idx_buff.data(), nlefs_to_sample);
        random_sample(fwd_lef_ranks.begin(), fwd_lef_ranks.end(), lef_idx.begin(), lef_idx.size(),
                      s.rand_eng);
        if (!std::all_of(lef_idx.begin(), lef_idx.end(),
                         [&](const auto i) { return i < lefs.size(); })) {
          throw std::runtime_error(fmt::format("lef_idx.size()={}; num_lefs={};\nlef_idx=[{}]\n",
                                               lef_idx.size(), lefs.size(),
                                               absl::StrJoin(lef_idx, ", ")));
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
                s.window_start, s.window_end, s.contacts, lefs, lef_idx, s.rand_eng,
                s.deletion_begin, s.deletion_size);
          } else {
            n_contacts +=
                this->register_contacts(s.window_start, s.window_end, s.contacts, lefs, lef_idx);
          }
        }

        if (s.num_target_contacts != 0 && n_contacts >= s.num_target_contacts) {
          return;  // Enough contact have been generated. Yay!
        }
      }

      this->generate_moves(*s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                           s.rand_eng);

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
      // if constexpr (normal_simulation) {
      Simulation::extrude(*s.chrom, lefs, rev_moves, fwd_moves, num_rev_units_at_5prime,
                          num_fwd_units_at_3prime);
      //} else {
      //  Simulation::extrude(*s.chrom, lefs, rev_moves, fwd_moves, num_rev_units_at_5prime,
      //                     num_fwd_units_at_3prime, s.deletion_begin, s.deletion_size);
      //}

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
