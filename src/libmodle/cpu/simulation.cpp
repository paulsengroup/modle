// IWYU pragma: private, include "modle/simulation.hpp"

#include "modle/simulation.hpp"

#include <H5Cpp.h>                              // IWYU pragma: keep
#include <absl/container/btree_set.h>           // for btree_iterator
#include <absl/strings/str_join.h>              // for StrJoin
#include <absl/types/span.h>                    // for Span, MakeSpan, MakeConstSpan
#include <cpp-sort/sorter_facade.h>             // for sorter_facade
#include <cpp-sort/sorters/counting_sorter.h>   // for counting_sort, counting_sorter
#include <cpp-sort/sorters/insertion_sorter.h>  // for insertion_sorter
#include <cpp-sort/sorters/pdq_sorter.h>        // for pdq_sort, pdq_sorter
#include <cpp-sort/sorters/ska_sorter.h>        // for ska_sort, ska_sorter
#include <cpp-sort/sorters/split_sorter.h>      // for split_sort, split_sorter
#include <fmt/format.h>                         // for format, print, FMT_STRING
#include <fmt/ostream.h>                        // for formatbuf<>::int_type

#include <algorithm>                                // for fill, min, max, clamp, for_each, gene...
#include <atomic>                                   // for atomic
#include <boost/asio/thread_pool.hpp>               // for thread_pool
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset, dynamic_bitset<>::ref...
#include <boost/range/adaptor/reversed.hpp>         // for reversed_range, reverse
#include <cassert>                                  // for assert
#include <chrono>                                   // for microseconds
#include <cmath>                                    // for round
#include <cstddef>                                  // for size_t
#include <cstdio>                                   // for stderr
#include <cstdlib>                                  // for abs
#include <deque>                                    // for deque
#include <filesystem>                               // for operator<<, path
#include <ios>                                      // for streamsize
#include <limits>                                   // for numeric_limits
#include <memory>                                   // for unique_ptr, make_unique, _MakeUniq<>:...
#include <mutex>                                    // for mutex
#include <numeric>                                  // for iota
#include <stdexcept>                                // for runtime_error
#include <string>                                   // for basic_string, string
#include <string_view>                              // for string_view
#include <utility>                                  // for make_pair, pair
#include <vector>                                   // for vector, vector<>::iterator

#include "modle/common/common.hpp"           // for BOOST_LIKELY, BOOST_UNLIKELY bp_t...
#include "modle/common/config.hpp"           // for Config
#include "modle/common/random_sampling.hpp"  // for random_sampling
#include "modle/common/utils.hpp"            // for ndebug_defined
#include "modle/cooler.hpp"                  // for Cooler, Cooler::WRITE_ONLY
#include "modle/extrusion_barriers.hpp"      // for update_states, ExtrusionBarrier
#include "modle/extrusion_factors.hpp"       // for Lef, ExtrusionUnit
#include "modle/genome.hpp"                  // for Chromosome, Genome

#ifndef BOOST_STACKTRACE_USE_NOOP
#include <boost/exception/get_error_info.hpp>  // for get_error_info
#include <boost/stacktrace/stacktrace.hpp>     // for operator<<
#include <ios>                                 // IWYU pragma: keep for streamsize
#endif

namespace modle {

Simulation::Simulation(const Config& c, bool import_chroms)
    : Config(c),
      _config(&c),
      _genome(import_chroms
                  ? Genome(path_to_chrom_sizes, path_to_extr_barriers, path_to_chrom_subranges,
                           path_to_feature_bed_files, ctcf_occupied_self_prob,
                           ctcf_not_occupied_self_prob, write_contacts_for_ko_chroms)
                  : Genome{}) {}

size_t Simulation::size() const { return this->_genome.size(); }

size_t Simulation::simulated_size() const { return this->_genome.simulated_size(); }

void Simulation::write_contacts_to_disk(std::deque<std::pair<Chromosome*, size_t>>& progress_queue,
                                        std::mutex& progress_queue_mutex,
                                        std::atomic<bool>& end_of_simulation)
    const {  // This thread is in charge of writing contacts to disk
  Chromosome* chrom_to_be_written = nullptr;
  const auto max_str_length =
      std::max_element(  // Find chrom with the longest name
          this->_genome.begin(), this->_genome.end(),
          [](const auto& c1, const auto& c2) { return c1.name().size() < c2.name().size(); })
          ->name()
          .size();

  auto c = this->skip_output ? nullptr
                             : std::make_unique<cooler::Cooler>(this->path_to_output_file_cool,
                                                                cooler::Cooler::WRITE_ONLY,
                                                                this->bin_size, max_str_length);

  auto sleep_us = 100;
  while (true) {  // Structuring the loop in this way allows us to sleep without holding the mutex
    sleep_us = std::min(500000, sleep_us * 2);  // NOLINT
    std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
    {
      std::scoped_lock l(progress_queue_mutex);
      if (progress_queue.empty()) {
        // There are no contacts to write to disk at the moment. Go back to sleep
        continue;
      }

      // chrom == nullptr is the end-of-queue signal
      if (auto& [chrom, count] = progress_queue.front(); chrom == nullptr) {
        end_of_simulation = true;
        return;
      }
      // count == ncells signals that we are done simulating the current chromosome
      else if (count == ncells) {  // NOLINT
        chrom_to_be_written = chrom;
        progress_queue.pop_front();
      } else {
        assert(count < ncells);  // NOLINT
        continue;
      }
    }
    sleep_us = 100;
    try {
      if (c) {  // c == nullptr only when --skip-output is used
        fmt::print(stderr, "Writing contacts for '{}' to file {}...\n", chrom_to_be_written->name(),
                   c->get_path());
        c->write_or_append_cmatrix_to_file(
            chrom_to_be_written->contacts(), chrom_to_be_written->name(),
            chrom_to_be_written->start_pos(), chrom_to_be_written->end_pos(),
            chrom_to_be_written->size(), true);
        fmt::print(stderr, "Written {} contacts for '{}' in {:.2f}M pixels to file {}.\n",
                   chrom_to_be_written->contacts().get_tot_contacts(), chrom_to_be_written->name(),
                   static_cast<double>(chrom_to_be_written->contacts().npixels()) / 1.0e6,
                   c->get_path());
      }
      // Deallocate the contact matrix to free up unused memory
      chrom_to_be_written->deallocate_contacts();
    } catch (const std::runtime_error& err) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("The following error occurred while writing contacts for '{}' to file {}: {}"),
          chrom_to_be_written->name(), c->get_path(), err.what()));
    }
  }
}

void Simulation::simulate_extrusion_kernel(Simulation::State& s) const {
  {
    assert(s.nlefs == s.lef_buff.size());     // NOLINT
    assert(s.nlefs == s.rank_buff1.size());   // NOLINT
    assert(s.nlefs == s.rank_buff2.size());   // NOLINT
    assert(s.nlefs == s.moves_buff1.size());  // NOLINT
    assert(s.nlefs == s.moves_buff2.size());  // NOLINT
    assert(s.nlefs == s.idx_buff1.size());    // NOLINT
    assert(s.nlefs == s.idx_buff2.size());    // NOLINT
    assert(s.nlefs == s.epoch_buff.size());   // NOLINT
  }
  auto epoch = 0UL;  // The epoch needs to be declared outside of the try-catch body so that we can
                     // use the current epoch when generating error messages
  try {
    // Seed is computed based on chrom. name, size and cellid
    s.seed = s.chrom->hash(this->seed, s.cell_id);
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
    s.n_target_epochs += s.n_target_contacts != 0 ? 0 : n_burnin_epochs;

    // Compute the avg. # of LEFs to use to sample contact every iterations
    const auto avg_nlefs_to_sample =
        static_cast<double>(s.nlefs) * this->lef_fraction_contact_sampling;

    // Compute the avg. # of LEFs to be released every iterations
    const auto avg_nlefs_to_release =
        static_cast<double>((this->rev_extrusion_speed + this->fwd_extrusion_speed) * s.nlefs) /
        static_cast<double>(this->average_lef_lifetime);

    // Local counter for the number of contacts generated by the current kernel instance
    auto n_contacts = 0UL;

    assert(avg_nlefs_to_sample <= static_cast<double>(s.nlefs));   // NOLINT
    assert(avg_nlefs_to_release <= static_cast<double>(s.nlefs));  // NOLINT
    assert((s.n_target_contacts != 0) ==                           // NOLINT
           (s.n_target_epochs == std::numeric_limits<size_t>::max()));

    // Declare the spans used by the burnin phase as well as the simulation itself
    auto lefs = absl::MakeSpan(s.lef_buff);
    auto lef_unloader_affinity = absl::MakeSpan(s.lef_unloader_affinity);
    const auto barriers = absl::MakeConstSpan(s.barriers);
    auto rev_lef_ranks = absl::MakeSpan(s.rank_buff1);
    auto fwd_lef_ranks = absl::MakeSpan(s.rank_buff2);
    auto rev_moves = absl::MakeSpan(s.moves_buff1);
    auto fwd_moves = absl::MakeSpan(s.moves_buff2);
    auto rev_collision_mask = absl::MakeSpan(s.idx_buff1);
    auto fwd_collision_mask = absl::MakeSpan(s.idx_buff2);

    // Generate initial extr. barrier states, so that they are already at or close to equilibrium
    for (auto i = 0UL; i < s.barrier_mask.size(); ++i) {
      s.barrier_mask[i] =
          random::bernoulli_trial{this->probability_of_extrusion_barrier_block}(s.rand_eng);
    }

    size_t nlefs_to_release;  // NOLINT
    // Start the burnin phase (followed by the actual simulation)
    for (; epoch < s.n_target_epochs; ++epoch) {
      if (!lef_initial_loading_epoch.empty()) {  // Execute this branch only during the burnin phase
        if (this->skip_burnin) {
          // The following will cause every LEF to be loaded in the current iteration
          lef_initial_loading_epoch.subspan(0, 0);
          // lef_initial_loading_epoch.clear();
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
        assert(lef_initial_loading_epoch.size() < s.nlefs);  // NOLINT

        // Compute the current number of active LEFs (i.e. LEFs that are either bound to DNA, or
        // that are valid candidates for re-binding). Grow spans accordingly
        const auto nlefs = s.nlefs - lef_initial_loading_epoch.size();

        lefs = absl::MakeSpan(s.lef_buff.data(), nlefs);
        lef_unloader_affinity = absl::MakeSpan(s.lef_unloader_affinity.data(), nlefs);
        rev_lef_ranks = absl::MakeSpan(s.rank_buff1.data(), nlefs);
        fwd_lef_ranks = absl::MakeSpan(s.rank_buff2.data(), nlefs);
        rev_moves = absl::MakeSpan(s.moves_buff1.data(), nlefs);
        fwd_moves = absl::MakeSpan(s.moves_buff2.data(), nlefs);
        rev_collision_mask = absl::MakeSpan(s.idx_buff1.data(), nlefs);
        fwd_collision_mask = absl::MakeSpan(s.idx_buff2.data(), nlefs);

        // Sample nlefs to be released while in burn-in phase
        nlefs_to_release =
            std::min(lefs.size(),
                     random::poisson_distribution<size_t>{
                         static_cast<double>(
                             (this->rev_extrusion_speed + this->fwd_extrusion_speed) * s.nlefs) /
                         static_cast<double>(this->average_lef_lifetime)}(s.rand_eng));
      } else {  // Sample nlefs to be released after the burn-in phase has been completed
        nlefs_to_release = std::min(
            lefs.size(), random::poisson_distribution<size_t>{avg_nlefs_to_release}(s.rand_eng));
      }

      ////////////////////////
      // Simulate one epoch //
      ////////////////////////

      {  // Select inactive LEFs and bind them
        auto lef_mask = absl::MakeSpan(s.idx_buff1.data(), lefs.size());
        Simulation::select_lefs_to_bind(lefs, lef_mask);
        Simulation::bind_lefs(s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, lef_mask, s.rand_eng,
                              epoch);
      }

      if (epoch > n_burnin_epochs) {              // Register contacts
        assert(fwd_lef_ranks.size() == s.nlefs);  // NOLINT

        auto nlefs_to_sample = std::min(
            lefs.size(), random::poisson_distribution<size_t>{avg_nlefs_to_sample}(s.rand_eng));

        if (s.n_target_contacts != 0) {  // When using the target contact density as stopping
                                         // criterion, don't overshoot the target number of contacts
          nlefs_to_sample = std::min(nlefs_to_sample, s.n_target_contacts - n_contacts);
        }

        // Select LEFs to be used for contact registration.
        // We are sampling from fwd ranks to avoid having to allocate a vector of indices just to do
        // this sampling
        const auto lef_idx = absl::MakeSpan(s.idx_buff1.data(), nlefs_to_sample);
        random_sample(fwd_lef_ranks.begin(), fwd_lef_ranks.end(), lef_idx.begin(), lef_idx.size(),
                      s.rand_eng);
        if (!std::all_of(lef_idx.begin(), lef_idx.end(),
                         [&](const auto i) { return i < lefs.size(); })) {
          throw std::runtime_error(fmt::format("lef_idx.size()={}; num_lefs={};\nlef_idx=[{}]\n",
                                               lef_idx.size(), lefs.size(),
                                               absl::StrJoin(lef_idx, ", ")));
        }
        n_contacts += this->register_contacts(s.chrom, lefs, lef_idx);

        if (s.n_target_contacts != 0 && n_contacts >= s.n_target_contacts) {
          return;  // Enough contact have been generated. Yay!
        }
      }

      this->generate_moves(s.chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                           s.rand_eng);

      CTCF::update_states(barriers, s.barrier_mask, s.rand_eng);

      // Reset collision masks
      std::fill(rev_collision_mask.begin(), rev_collision_mask.end(), NO_COLLISION);
      std::fill(fwd_collision_mask.begin(), fwd_collision_mask.end(), NO_COLLISION);

      // Detect collision and correct moves
      const auto& [num_rev_units_at_5prime, num_fwd_units_at_3prime] =
          Simulation::process_collisions(s.chrom, lefs, barriers, s.barrier_mask, rev_lef_ranks,
                                         fwd_lef_ranks, rev_moves, fwd_moves, rev_collision_mask,
                                         fwd_collision_mask, s.rand_eng);

      // Advance LEFs
      Simulation::extrude(s.chrom, lefs, rev_moves, fwd_moves, num_rev_units_at_5prime,
                          num_fwd_units_at_3prime);

      // The vector of affinities is used to bias LEF release towards LEFs that are not in a hard
      // stall condition
      this->generate_lef_unloader_affinities(lefs, barriers, rev_collision_mask, fwd_collision_mask,
                                             lef_unloader_affinity);

      // Reusing this buffer is ok, as at this point we don't need access to collision information
      const auto lef_idx = absl::MakeSpan(s.idx_buff1.data(), nlefs_to_release);
      Simulation::select_lefs_to_release(lef_idx, lef_unloader_affinity, s.rand_eng);
      Simulation::release_lefs(lefs, lef_idx);
    }
  } catch (const std::exception& err) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Caught an exception while simulating epoch #{} in cell #{} of {}:\n{}"), epoch,
        s.cell_id, s.chrom->name(), err.what()));
  }
}

bp_t Simulation::generate_rev_move(const Chromosome* const chrom, const ExtrusionUnit& unit,
                                   random::PRNG_t& rand_eng) const {
  assert(unit.pos() >= chrom->start_pos());  // NOLINT
  if (this->rev_extrusion_speed_std == 0) {  // When std == 0 always return the avg. extrusion speed
    // (except when unit is close to chrom start pos.)
    return std::min(this->rev_extrusion_speed, unit.pos() - chrom->start_pos());
  }
  // Generate the move distance (and make sure it is not a negative distance)
  // NOTE: on my laptop generating doubles from a normal distribution, rounding them, then
  // casting double to uint is a lot faster (~4x) than drawing uints directly from a Poisson
  // distr.
  return std::clamp(static_cast<bp_t>(std::round(
                        lef_move_generator_t{static_cast<double>(this->rev_extrusion_speed),
                                             this->rev_extrusion_speed_std}(rand_eng))),
                    0UL, unit.pos() - chrom->start_pos());
}

bp_t Simulation::generate_fwd_move(const Chromosome* const chrom, const ExtrusionUnit& unit,
                                   random::PRNG_t& rand_eng) const {
  // See Simulation::generate_rev_move for comments
  assert(unit.pos() < chrom->end_pos());  // NOLINT
  if (this->fwd_extrusion_speed_std == 0) {
    return std::min(this->fwd_extrusion_speed, (chrom->end_pos() - 1) - unit.pos());
  }
  return std::clamp(static_cast<bp_t>(std::round(
                        lef_move_generator_t{static_cast<double>(this->fwd_extrusion_speed),
                                             this->fwd_extrusion_speed_std}(rand_eng))),
                    0UL, (chrom->end_pos() - 1) - unit.pos());
}

void Simulation::generate_moves(const Chromosome* const chrom, const absl::Span<const Lef> lefs,
                                const absl::Span<const size_t> rev_lef_ranks,
                                const absl::Span<const size_t> fwd_lef_ranks,
                                const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
                                random::PRNG_t& rand_eng, bool adjust_moves_) const
    noexcept(utils::ndebug_defined()) {
  {
    assert(lefs.size() == fwd_lef_ranks.size());  // NOLINT
    assert(lefs.size() == rev_lef_ranks.size());  // NOLINT
    assert(lefs.size() == fwd_moves.size());      // NOLINT
    assert(lefs.size() == rev_moves.size());      // NOLINT
  }

  // As long as a LEF is bound to DNA, always generate a move
  for (auto i = 0UL; i < lefs.size(); ++i) {
    rev_moves[i] = lefs[i].is_bound() ? generate_rev_move(chrom, lefs[i].rev_unit, rand_eng) : 0UL;
    fwd_moves[i] = lefs[i].is_bound() ? generate_fwd_move(chrom, lefs[i].fwd_unit, rand_eng) : 0UL;
  }

  if (adjust_moves_) {  // Adjust moves of consecutive extr. units to make LEF behavior more
                        // realistic See comments in adjust_moves_of_consecutive_extr_units for more
                        // details on what this entails
    this->adjust_moves_of_consecutive_extr_units(chrom, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                 rev_moves, fwd_moves);
  }
}

void Simulation::adjust_moves_of_consecutive_extr_units(
    const Chromosome* chrom, absl::Span<const Lef> lefs, absl::Span<const size_t> rev_lef_ranks,
    absl::Span<const size_t> fwd_lef_ranks, absl::Span<bp_t> rev_moves,
    absl::Span<bp_t> fwd_moves) const noexcept(utils::ndebug_defined()) {
  (void)chrom;

  // Loop over pairs of consecutive extr. units.
  // Extr. units moving in rev direction are processed in 3'-5' order, while units moving in fwd
  // direction are processed in 5'-3' direction
  const auto rev_offset = lefs.size() - 1;
  for (auto i = 0UL; i < lefs.size() - 1; ++i) {
    const auto& idx1 = rev_lef_ranks[rev_offset - 1 - i];
    const auto& idx2 = rev_lef_ranks[rev_offset - i];

    if (lefs[idx1].is_bound() && lefs[idx2].is_bound()) {
      assert(lefs[idx1].rev_unit.pos() >= chrom->start_pos() + rev_moves[idx1]);  // NOLINT
      assert(lefs[idx2].rev_unit.pos() >= chrom->start_pos() + rev_moves[idx2]);  // NOLINT

      const auto pos1 = lefs[idx1].rev_unit.pos() - rev_moves[idx1];
      const auto pos2 = lefs[idx2].rev_unit.pos() - rev_moves[idx2];

      // If moving extr. units 1 and 2 by their respective moves would cause unit 1 (which comes
      // first in 3'-5' order) to be surpassed by unit 2 (which comes second in 3'-5' order),
      // increment the move for unit 1 so that after moving both units, unit 1 will be located one
      // bp upstream of unit 2 in 5'-3' direction.
      // This mimics what would probably happen in a real system, where extr. unit 2 would most
      // likely push extr. unit 1, temporarily increasing extr. speed of unit 1.
      if (pos2 < pos1) {
        rev_moves[idx1] += pos1 - pos2;
      }
    }

    const auto& idx3 = fwd_lef_ranks[i];
    const auto& idx4 = fwd_lef_ranks[i + 1];

    // See above for detailed comments. The logic is the same used on rev units (but mirrored!)
    if (lefs[idx3].is_bound() && lefs[idx4].is_bound()) {
      assert(lefs[idx3].fwd_unit.pos() + fwd_moves[idx3] < chrom->end_pos());  // NOLINT
      assert(lefs[idx4].fwd_unit.pos() + fwd_moves[idx4] < chrom->end_pos());  // NOLINT

      const auto pos3 = lefs[idx3].fwd_unit.pos() + fwd_moves[idx3];
      const auto pos4 = lefs[idx4].fwd_unit.pos() + fwd_moves[idx4];

      if (pos3 > pos4) {
        fwd_moves[idx4] += pos3 - pos4;
      }
    }
  }
}

void Simulation::rank_lefs(const absl::Span<const Lef> lefs,
                           const absl::Span<size_t> rev_lef_rank_buff,
                           const absl::Span<size_t> fwd_lef_rank_buff,
                           bool ranks_are_partially_sorted,
                           bool init_buffers) noexcept(utils::ndebug_defined()) {
  assert(lefs.size() == fwd_lef_rank_buff.size());  // NOLINT
  assert(lefs.size() == rev_lef_rank_buff.size());  // NOLINT

  auto rev_comparator = [&](const auto r1, const auto r2) constexpr noexcept {
    assert(r1 < lefs.size());  // NOLINT
    assert(r2 < lefs.size());  // NOLINT
    return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
  };

  // See comments for rev_comparator.
  auto fwd_comparator = [&](const auto r1, const auto r2) constexpr noexcept {
    assert(r1 < lefs.size());  // NOLINT
    assert(r2 < lefs.size());  // NOLINT
    return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
  };

  if (BOOST_UNLIKELY(init_buffers)) {  // Init rank buffers
    std::iota(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), 0);
    std::iota(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), 0);
  }

  if (BOOST_LIKELY(ranks_are_partially_sorted)) {
    cppsort::split_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), rev_comparator);
    cppsort::split_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), fwd_comparator);
  } else {
    // Fallback to pattern-defeating quicksort we have no information regarding the level of
    // pre-sortedness of LEFs
    cppsort::pdq_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), rev_comparator);
    cppsort::pdq_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), fwd_comparator);
  }

  // TODO Figure out a better way to deal with ties
  auto begin = 0UL;
  auto end = 0UL;
  for (auto i = 1UL; i < rev_lef_rank_buff.size(); ++i) {
    const auto& r1 = rev_lef_rank_buff[i - 1];
    const auto& r2 = rev_lef_rank_buff[i];
    if (BOOST_UNLIKELY(lefs[r1].rev_unit.pos() == lefs[r2].rev_unit.pos())) {
      begin = i - 1;
      for (; i < rev_lef_rank_buff.size(); ++i) {
        const auto& r11 = rev_lef_rank_buff[i - 1];
        const auto& r22 = rev_lef_rank_buff[i];
        if (lefs[r11].rev_unit.pos() != lefs[r22].rev_unit.pos()) {
          break;
        }
      }

      end = i;
      cppsort::insertion_sort(rev_lef_rank_buff.begin() + begin, rev_lef_rank_buff.begin() + end,
                              [&lefs](const auto r11, const auto r22) {
                                assert(r11 < lefs.size());  // NOLINT
                                assert(r22 < lefs.size());  // NOLINT
                                return lefs[r11].binding_epoch < lefs[r22].binding_epoch;
                              });
      begin = end;
    }
  }

  for (auto i = 1UL; i < fwd_lef_rank_buff.size(); ++i) {
    const auto& r1 = fwd_lef_rank_buff[i - 1];
    const auto& r2 = fwd_lef_rank_buff[i];
    if (BOOST_UNLIKELY(lefs[r1].fwd_unit.pos() == lefs[r2].fwd_unit.pos())) {
      begin = i - 1;
      for (; i < fwd_lef_rank_buff.size(); ++i) {
        const auto& r11 = fwd_lef_rank_buff[i - 1];
        const auto& r22 = fwd_lef_rank_buff[i];
        if (lefs[r11].fwd_unit.pos() != lefs[r22].fwd_unit.pos()) {
          break;
        }
      }

      end = i;
      cppsort::insertion_sort(fwd_lef_rank_buff.begin() + begin, fwd_lef_rank_buff.begin() + end,
                              [&lefs](const auto r11, const auto r22) {
                                assert(r11 < lefs.size());  // NOLINT
                                assert(r22 < lefs.size());  // NOLINT
                                return lefs[r22].binding_epoch < lefs[r11].binding_epoch;
                              });
      begin = end;
    }
  }
}

void Simulation::extrude(const Chromosome* chrom, const absl::Span<Lef> lefs,
                         const absl::Span<const bp_t> rev_moves,
                         const absl::Span<const bp_t> fwd_moves, size_t num_rev_units_at_5prime,
                         size_t num_fwd_units_at_3prime) noexcept(utils::ndebug_defined()) {
  {
    assert(lefs.size() == rev_moves.size());         // NOLINT
    assert(lefs.size() == fwd_moves.size());         // NOLINT
    assert(lefs.size() >= num_rev_units_at_5prime);  // NOLINT
    assert(lefs.size() >= num_fwd_units_at_3prime);  // NOLINT
    (void)chrom;
  }

  auto i1 = num_rev_units_at_5prime == 0 ? 0UL : num_rev_units_at_5prime - 1;
  const auto i2 = lefs.size() - num_fwd_units_at_3prime;
  for (; i1 < i2; ++i1) {
    auto& lef = lefs[i1];
    if (BOOST_UNLIKELY(!lef.is_bound())) {  // Do not process inactive LEFs
      continue;
    }
    assert(lef.rev_unit.pos() <= lef.fwd_unit.pos());  // NOLINT

    // Extrude rev unit
    assert(lef.rev_unit.pos() >= chrom->start_pos() + rev_moves[i1]);  // NOLINT
    lef.rev_unit._pos -= rev_moves[i1];  // Advance extr. unit in 3'-5' direction

    // Extrude fwd unit
    assert(lef.fwd_unit.pos() + fwd_moves[i1] <= chrom->end_pos() - 1);  // NOLINT
    lef.fwd_unit._pos += fwd_moves[i1];  // Advance extr. unit in 5'-3' direction

    assert(lef.rev_unit.pos() <= lef.fwd_unit.pos());  // NOLINT
  }
}

std::pair<bp_t, bp_t> Simulation::compute_lef_lef_collision_pos(const ExtrusionUnit& rev_unit,
                                                                const ExtrusionUnit& fwd_unit,
                                                                bp_t rev_move, bp_t fwd_move) {
  const auto& rev_speed = rev_move;
  const auto& fwd_speed = fwd_move;
  const auto& rev_pos = rev_unit.pos();
  const auto& fwd_pos = fwd_unit.pos();

  const auto relative_speed = rev_speed + fwd_speed;
  const auto time_to_collision =
      static_cast<double>(rev_pos - fwd_pos) / static_cast<double>(relative_speed);

  const auto collision_pos =
      fwd_pos + static_cast<bp_t>(std::round((static_cast<double>(fwd_speed) * time_to_collision)));
  assert(collision_pos <= rev_pos);  // NOLINT
#ifndef NDEBUG
  const auto collision_pos_ =
      static_cast<double>(rev_pos) - static_cast<double>(rev_speed) * time_to_collision;
  assert(abs(static_cast<double>(collision_pos) - collision_pos_) < 1.0);  // NOLINT
#endif
  if (BOOST_UNLIKELY(collision_pos == fwd_pos)) {
    assert(collision_pos >= fwd_pos);      // NOLINT
    assert(collision_pos + 1 <= rev_pos);  // NOLINT
    return std::make_pair(collision_pos + 1, collision_pos);
  }
  assert(collision_pos > 0);             // NOLINT
  assert(collision_pos - 1 >= fwd_pos);  // NOLINT
  return std::make_pair(collision_pos, collision_pos - 1);
}

size_t Simulation::register_contacts(Chromosome* chrom, const absl::Span<const Lef> lefs,
                                     const absl::Span<const size_t> selected_lef_idx) const
    noexcept(utils::ndebug_defined()) {
  // Register contacts for the selected LEFs (excluding LEFs that have one of their units at the
  // beginning/end of a chromosome)
  size_t new_contacts = 0;
  for (const auto i : selected_lef_idx) {
    assert(i < lefs.size());  // NOLINT
    const auto& lef = lefs[i];
    if (BOOST_LIKELY(lef.is_bound() && lef.rev_unit.pos() > chrom->start_pos() &&
                     lef.rev_unit.pos() < chrom->end_pos() - 1 &&
                     lef.fwd_unit.pos() > chrom->start_pos() &&
                     lef.fwd_unit.pos() < chrom->end_pos() - 1)) {
      chrom->increment_contacts(lef.rev_unit.pos(), lef.fwd_unit.pos(), this->bin_size);
      ++new_contacts;
    }
  }
  return new_contacts;
}

void Simulation::generate_lef_unloader_affinities(
    const absl::Span<const Lef> lefs, const absl::Span<const ExtrusionBarrier> barriers,
    const absl::Span<const collision_t> rev_collisions,
    const absl::Span<const collision_t> fwd_collisions,
    const absl::Span<double> lef_unloader_affinity) const noexcept(utils::ndebug_defined()) {
  assert(lefs.size() == rev_collisions.size());         // NOLINT
  assert(lefs.size() == fwd_collisions.size());         // NOLINT
  assert(lefs.size() == lef_unloader_affinity.size());  // NOLINT

  auto is_lef_bar_collision = [&](const auto i) { return i < barriers.size(); };

  for (auto i = 0UL; i < lefs.size(); ++i) {
    const auto& lef = lefs[i];
    if (!lef.is_bound()) {
      lef_unloader_affinity[i] = 0.0;
    } else if (BOOST_LIKELY(!is_lef_bar_collision(rev_collisions[i]) ||
                            !is_lef_bar_collision(fwd_collisions[i]))) {
      lef_unloader_affinity[i] = 1.0;
    } else {
      const auto& rev_barrier = barriers[rev_collisions[i]];
      const auto& fwd_barrier = barriers[fwd_collisions[i]];

      if (BOOST_UNLIKELY(rev_barrier.blocking_direction_major() == dna::rev &&
                         fwd_barrier.blocking_direction_major() == dna::fwd)) {
        lef_unloader_affinity[i] = 1.0 / this->hard_stall_multiplier;
      } else {
        lef_unloader_affinity[i] = 1.0;
      }
    }
  }
}

void Simulation::select_lefs_to_release(
    const absl::Span<size_t> lef_idx, const absl::Span<const double> lef_unloader_affinity,
    random::PRNG_t& rand_eng) noexcept(utils::ndebug_defined()) {
  random::discrete_distribution<size_t> idx_gen(lef_unloader_affinity.begin(),
                                                lef_unloader_affinity.end());
  std::generate(lef_idx.begin(), lef_idx.end(), [&]() { return idx_gen(rand_eng); });
}

void Simulation::release_lefs(const absl::Span<Lef> lefs,
                              const absl::Span<const size_t> lef_idx) noexcept {
  for (const auto i : lef_idx) {
    assert(i < lefs.size());  // NOLINT
    lefs[i].release();
  }
}

boost::asio::thread_pool Simulation::instantiate_thread_pool() const {
  return boost::asio::thread_pool(this->nthreads);
}

Simulation::State& Simulation::State::operator=(const Task& task) {
  this->id = task.id;
  this->chrom = task.chrom;
  this->cell_id = task.cell_id;
  this->n_target_epochs = task.n_target_epochs;
  this->n_target_contacts = task.n_target_contacts;
  this->nlefs = task.nlefs;
  this->barriers = task.barriers;
  return *this;
}

void Simulation::State::resize(size_t size) {
  if (size == std::numeric_limits<size_t>::max()) {
    size = this->nlefs;
  }
  lef_buff.resize(size);
  lef_unloader_affinity.resize(size);
  rank_buff1.resize(size);
  rank_buff2.resize(size);
  barrier_mask.resize(this->barriers.size());
  moves_buff1.resize(size);
  moves_buff2.resize(size);
  idx_buff1.resize(size);
  idx_buff2.resize(size);
  epoch_buff.resize(size);
}

void Simulation::State::reset() {  // TODO figure out which resets are redundant
  std::for_each(lef_buff.begin(), lef_buff.end(), [](auto& lef) { lef.reset(); });
  std::fill(lef_unloader_affinity.begin(), lef_unloader_affinity.end(), 0.0);
  std::iota(rank_buff1.begin(), rank_buff1.end(), 0);
  std::copy(rank_buff1.begin(), rank_buff1.end(), rank_buff2.begin());
  barrier_mask.reset();
  std::fill(moves_buff1.begin(), moves_buff1.end(), 0);
  std::fill(moves_buff2.begin(), moves_buff2.end(), 0);
  std::fill(idx_buff1.begin(), idx_buff1.end(), 0);
  std::fill(idx_buff2.begin(), idx_buff2.end(), 0);
  std::fill(epoch_buff.begin(), epoch_buff.end(), 0);
}

std::string Simulation::State::to_string() const noexcept {
  return fmt::format(FMT_STRING("State:\n - TaskID {}\n"
                                " - Chrom: {}[{}-{}]\n"
                                " - CellID: {}\n"
                                " - Target epochs: {}\n"
                                " - Target contacts: {}\n"
                                " - # of LEFs: {}\n"
                                " - # Extrusion barriers: {}\n"
                                " - seed: {}\n"),
                     id, chrom->name(), chrom->start_pos(), chrom->end_pos(), cell_id,
                     n_target_epochs, n_target_contacts, nlefs, barriers.size(), seed);
}

std::pair<size_t, size_t> Simulation::process_collisions(
    const Chromosome* chrom, const absl::Span<const Lef> lefs,
    const absl::Span<const ExtrusionBarrier> barriers, const boost::dynamic_bitset<>& barrier_mask,
    const absl::Span<const size_t> rev_lef_ranks, const absl::Span<const size_t> fwd_lef_ranks,
    const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
    const absl::Span<collision_t> rev_collisions, const absl::Span<collision_t> fwd_collisions,
    random::PRNG_t& rand_eng) const noexcept(utils::ndebug_defined()) {
  const auto& [num_rev_units_at_5prime, num_fwd_units_at_3prime] =
      Simulation::detect_units_at_chrom_boundaries(chrom, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                   rev_moves, fwd_moves, rev_collisions,
                                                   fwd_collisions);

  this->detect_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                                  barriers, barrier_mask, rev_collisions, fwd_collisions, rand_eng,
                                  num_rev_units_at_5prime, num_fwd_units_at_3prime);

  this->detect_primary_lef_lef_collisions(lefs, barriers, rev_lef_ranks, fwd_lef_ranks, rev_moves,
                                          fwd_moves, rev_collisions, fwd_collisions, rand_eng,
                                          num_rev_units_at_5prime, num_fwd_units_at_3prime);
  Simulation::correct_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                   rev_collisions, fwd_collisions);

  Simulation::correct_moves_for_primary_lef_lef_collisions(lefs, barriers, rev_lef_ranks,
                                                           fwd_lef_ranks, rev_moves, fwd_moves,
                                                           rev_collisions, fwd_collisions);
  this->process_secondary_lef_lef_collisions(
      chrom, lefs, barriers.size(), rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
      rev_collisions, fwd_collisions, rand_eng, num_rev_units_at_5prime, num_fwd_units_at_3prime);
  return std::make_pair(num_rev_units_at_5prime, num_fwd_units_at_3prime);
}

}  // namespace modle
