// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// IWYU pragma: private, include "modle/simulation.hpp"

#include "modle/simulation.hpp"

#include <absl/container/btree_set.h>           // for btree_iterator
#include <absl/strings/str_split.h>             // for StrSplit, Splitter
#include <absl/types/span.h>                    // for Span, MakeConstSpan, MakeSpan
#include <cpp-sort/sorter_facade.h>             // for sorter_facade
#include <cpp-sort/sorters/insertion_sorter.h>  // for insertion_sort, insertion_so...
#include <cpp-sort/sorters/pdq_sorter.h>        // for pdq_sort, pdq_sorter
#include <cpp-sort/sorters/split_sorter.h>      // for split_sort, split_sorter
#include <fmt/ostream.h>                        // for formatbuf<>::int_type
#include <spdlog/spdlog.h>                      // for info, warn

#include <algorithm>         // for max, fill, min, copy, clamp
#include <atomic>            // for atomic
#include <boost/config.hpp>  // IWYU pragma: keep for BOOST_LIKELY, BOOST_UNLIKELY
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset, dynamic_bits...
#include <boost/filesystem/path.hpp>                // for operator<<, path
#include <cassert>                                  // for assert
#include <chrono>                                   // for microseconds
#include <cmath>                                    // for log, round, exp, floor, sqrt
#include <cstddef>                                  // for size_t, ptrdiff_t
#include <cstdlib>                                  // for abs
#include <deque>                                    // for _Deque_iterator<>::_Self
#include <iosfwd>                                   // for streamsize
#include <limits>                                   // for numeric_limits
#include <memory>                                   // for shared_ptr, unique_ptr, make...
#include <mutex>                                    // for mutex
#include <numeric>                                  // for iota
#include <stdexcept>                                // for runtime_error
#include <string>                                   // for string
#include <string_view>                              // for string_view
#include <thread>                                   // for sleep_for
#include <thread_pool/thread_pool.hpp>              // for thread_pool
#include <utility>                                  // for make_pair, pair
#include <vector>                                   // for vector, vector<>::iterator

#include "modle/common/common.hpp"                         // for bp_t, collision_t, contacts_t
#include "modle/common/config.hpp"                         // for Config
#include "modle/common/genextreme_value_distribution.hpp"  // for genextreme_value_distribution
#include "modle/common/random.hpp"           // for bernoulli_trial, poisson_distribution
#include "modle/common/random_sampling.hpp"  // for random_sample
#include "modle/common/utils.hpp"            // for parse_numeric_or_throw, ndeb...
#include "modle/cooler.hpp"                  // for Cooler, Cooler::WRITE_ONLY
#include "modle/extrusion_barriers.hpp"      // for ExtrusionBarrier, update_states
#include "modle/extrusion_factors.hpp"       // for Lef, ExtrusionUnit
#include "modle/genome.hpp"                  // for Genome::iterator, Chromosome
#include "modle/interval_tree.hpp"           // for IITree, IITree::data

namespace modle {

Simulation::Simulation(const Config& c, bool import_chroms)
    : Config(c),
      _genome(import_chroms
                  ? Genome(path_to_chrom_sizes, path_to_extr_barriers, path_to_chrom_subranges,
                           path_to_feature_bed_files, ctcf_occupied_self_prob,
                           ctcf_not_occupied_self_prob, write_contacts_for_ko_chroms)
                  : Genome{}),
      _tpool(c.nthreads + 1) {}

bool Simulation::ok() const noexcept { return !this->_exception_thrown; }

size_t Simulation::size() const { return this->_genome.size(); }

size_t Simulation::simulated_size() const { return this->_genome.simulated_size(); }

void Simulation::write_contacts_to_disk(std::deque<std::pair<Chromosome*, size_t>>& progress_queue,
                                        std::mutex& progress_queue_mutex) {
  // This thread is in charge of writing contacts to disk
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

  try {
    if (c && !this->argv_json.empty()) {
      c->write_metadata_attribute(this->argv_json);
    }
    // NOLINTNEXTLINE(readability-magic-numbers), cppcoreguidelines-avoid-magic-numbers)
    auto sleep_us = 100;
    while (this->ok()) {  // Structuring the loop in this way allows us to sleep without
                          // holding the mutex
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
          break;
        }
        // count == ncells signals that we are done simulating the current chromosome
        else if (count == num_cells) {  // NOLINT
          chrom_to_be_written = chrom;
          progress_queue.pop_front();
        } else {
          assert(count < num_cells);  // NOLINT
          continue;
        }
      }
      sleep_us = 100;  // NOLINT(readability-magic-numbers), cppcoreguidelines-avoid-magic-numbers)
      if (c) {         // c == nullptr only when --skip-output is used
        // NOTE here we have to use pointers instead of references because
        // chrom_to_be_written.contacts() == nullptr is used to signal an empty matrix.
        // In this case, c->write_or_append_cmatrix_to_file() will create an entry in the chroms and
        // bins datasets, as well as update the appropriate index
        if (chrom_to_be_written->contacts_ptr()) {
          spdlog::info(FMT_STRING("Writing contacts for '{}' to file {}..."),
                       chrom_to_be_written->name(), c->get_path());
        } else {
          spdlog::info(FMT_STRING("Creating an empty entry for '{}' in file {}..."),
                       chrom_to_be_written->name(), c->get_path());
        }

        c->write_or_append_cmatrix_to_file(
            chrom_to_be_written->contacts_ptr().get(), chrom_to_be_written->name(),
            chrom_to_be_written->start_pos(), chrom_to_be_written->end_pos(),
            chrom_to_be_written->size());

        if (chrom_to_be_written->contacts_ptr()) {
          spdlog::info(
              FMT_STRING("Written {} contacts for '{}' in {:.2f}M pixels to file {}."),
              chrom_to_be_written->contacts().get_tot_contacts(), chrom_to_be_written->name(),
              static_cast<double>(chrom_to_be_written->contacts().npixels()) / 1.0e6,  // NOLINT
              c->get_path());
        } else {
          spdlog::info(FMT_STRING("Created an entry for '{}' in file {}."),
                       chrom_to_be_written->name(), c->get_path());
        }
      }
      // Deallocate the contact matrix to free up unused memory
      chrom_to_be_written->deallocate_contacts();
    }
  } catch (const std::exception& err) {
    std::scoped_lock l(this->_exceptions_mutex);
    if (chrom_to_be_written) {
      this->_exceptions.emplace_back(std::make_exception_ptr(std::runtime_error(fmt::format(
          FMT_STRING("The following error occurred while writing contacts for '{}' to file {}: {}"),
          chrom_to_be_written->name(), c->get_path(), err.what()))));
    } else {
      this->_exceptions.emplace_back(std::make_exception_ptr(std::runtime_error(fmt::format(
          FMT_STRING("The following error occurred while writing contacts to file {}: {}"),
          c->get_path(), err.what()))));
    }
    this->_exception_thrown = true;
  } catch (...) {
    std::scoped_lock l(this->_exceptions_mutex);
    this->_exceptions.emplace_back(std::make_exception_ptr(
        std::runtime_error("An unhandled exception was caught! This should never happen! "
                           "If you see this message, please file an issue on GitHub.")));
    this->_exception_thrown = true;
  }
  this->_end_of_simulation = true;
}

bp_t Simulation::generate_rev_move(const Chromosome& chrom, const ExtrusionUnit& unit,
                                   const double avg_extr_speed, const double extr_speed_std,
                                   random::PRNG_t& rand_eng) {
  assert(unit.pos() >= chrom.start_pos());  // NOLINT
  if (extr_speed_std == 0.0) {              // When std == 0 always return the avg. extrusion speed
    // (except when unit is close to chrom start pos.)
    return std::min(static_cast<bp_t>(avg_extr_speed), unit.pos() - chrom.start_pos());
  }
  // Generate the move distance (and make sure it is not a negative distance)
  // NOTE: on my laptop generating doubles from a normal distribution, rounding them, then
  // casting double to uint is a lot faster (~4x) than drawing uints directly from a Poisson
  // distr.
  return std::clamp(
      static_cast<bp_t>(std::round(lef_move_generator_t{avg_extr_speed, extr_speed_std}(rand_eng))),
      bp_t(0), unit.pos() - chrom.start_pos());
}

bp_t Simulation::generate_fwd_move(const Chromosome& chrom, const ExtrusionUnit& unit,
                                   const double avg_extr_speed, const double extr_speed_std,
                                   random::PRNG_t& rand_eng) {
  // See Simulation::generate_rev_move for comments
  assert(unit.pos() < chrom.end_pos());  // NOLINT
  if (extr_speed_std == 0.0) {
    return std::min(static_cast<bp_t>(avg_extr_speed), (chrom.end_pos() - 1) - unit.pos());
  }
  return std::clamp(
      static_cast<bp_t>(std::round(lef_move_generator_t{avg_extr_speed, extr_speed_std}(rand_eng))),
      bp_t(0), (chrom.end_pos() - 1) - unit.pos());
}

void Simulation::generate_moves(const Chromosome& chrom, const absl::Span<const Lef> lefs,
                                const absl::Span<const size_t> rev_lef_ranks,
                                const absl::Span<const size_t> fwd_lef_ranks,
                                const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
                                const bool burnin_completed, random::PRNG_t& rand_eng,
                                bool adjust_moves_) const noexcept(utils::ndebug_defined()) {
  {
    assert(lefs.size() == fwd_lef_ranks.size());  // NOLINT
    assert(lefs.size() == rev_lef_ranks.size());  // NOLINT
    assert(lefs.size() == fwd_moves.size());      // NOLINT
    assert(lefs.size() == rev_moves.size());      // NOLINT
  }

  // As long as a LEF is bound to DNA, always generate a move
  const auto rev_extr_speed = static_cast<double>(
      burnin_completed ? this->rev_extrusion_speed : this->rev_extrusion_speed_burnin);
  const auto fwd_extr_speed = static_cast<double>(
      burnin_completed ? this->fwd_extrusion_speed : this->fwd_extrusion_speed_burnin);
  for (auto i = 0UL; i < lefs.size(); ++i) {
    rev_moves[i] = lefs[i].is_bound() ? generate_rev_move(chrom, lefs[i].rev_unit, rev_extr_speed,
                                                          this->rev_extrusion_speed_std, rand_eng)
                                      : 0UL;
    fwd_moves[i] = lefs[i].is_bound() ? generate_fwd_move(chrom, lefs[i].fwd_unit, fwd_extr_speed,
                                                          this->fwd_extrusion_speed_std, rand_eng)
                                      : 0UL;
  }

  if (adjust_moves_) {  // Adjust moves of consecutive extr. units to make LEF behavior more
    // realistic See comments in adjust_moves_of_consecutive_extr_units for more
    // details on what this entails
    Simulation::adjust_moves_of_consecutive_extr_units(chrom, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                       rev_moves, fwd_moves);
  }
}

void Simulation::adjust_moves_of_consecutive_extr_units(
    [[maybe_unused]] const Chromosome& chrom, absl::Span<const Lef> lefs,
    absl::Span<const size_t> rev_lef_ranks, absl::Span<const size_t> fwd_lef_ranks,
    absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves) noexcept(utils::ndebug_defined()) {
  assert(!lefs.empty());  // NOLINT

  // Loop over pairs of consecutive extr. units.
  // Extr. units moving in rev direction are processed in 3'-5' order, while units moving in fwd
  // direction are processed in 5'-3' direction
  const auto rev_offset = lefs.size() - 1;
  for (auto i = 0UL; i < lefs.size() - 1; ++i) {
    {
      const auto idx1 = rev_lef_ranks[rev_offset - 1 - i];
      const auto idx2 = rev_lef_ranks[rev_offset - i];

      if (lefs[idx1].is_bound() && lefs[idx2].is_bound()) {
        assert(lefs[idx1].rev_unit.pos() >= chrom.start_pos() + rev_moves[idx1]);  // NOLINT
        assert(lefs[idx2].rev_unit.pos() >= chrom.start_pos() + rev_moves[idx2]);  // NOLINT

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
    }

    const auto idx1 = fwd_lef_ranks[i];
    const auto idx2 = fwd_lef_ranks[i + 1];

    // See above for detailed comments. The logic is the same used on rev units (but mirrored!)
    if (lefs[idx1].is_bound() && lefs[idx2].is_bound()) {
      assert(lefs[idx1].fwd_unit.pos() + fwd_moves[idx1] < chrom.end_pos());  // NOLINT
      assert(lefs[idx1].fwd_unit.pos() + fwd_moves[idx2] < chrom.end_pos());  // NOLINT

      const auto pos1 = lefs[idx1].fwd_unit.pos() + fwd_moves[idx1];
      const auto pos2 = lefs[idx2].fwd_unit.pos() + fwd_moves[idx2];

      if (pos1 > pos2) {
        fwd_moves[idx2] += pos1 - pos2;
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

void Simulation::extrude([[maybe_unused]] const Chromosome& chrom, const absl::Span<Lef> lefs,
                         const absl::Span<const bp_t> rev_moves,
                         const absl::Span<const bp_t> fwd_moves,
                         const size_t num_rev_units_at_5prime,
                         const size_t num_fwd_units_at_3prime) noexcept(utils::ndebug_defined()) {
  {
    assert(lefs.size() == rev_moves.size());         // NOLINT
    assert(lefs.size() == fwd_moves.size());         // NOLINT
    assert(lefs.size() >= num_rev_units_at_5prime);  // NOLINT
    assert(lefs.size() >= num_fwd_units_at_3prime);  // NOLINT
  }

  auto i1 = num_rev_units_at_5prime == 0 ? 0UL : num_rev_units_at_5prime - 1;
  const auto i2 = lefs.size() - num_fwd_units_at_3prime;
  for (; i1 < i2; ++i1) {
    auto& lef = lefs[i1];
    if (BOOST_UNLIKELY(!lef.is_bound())) {  // Do not process inactive LEFs
      continue;
    }
    assert(lef.rev_unit.pos() <= lef.fwd_unit.pos());                   // NOLINT
    assert(lef.rev_unit.pos() >= chrom.start_pos() + rev_moves[i1]);    // NOLINT
    assert(lef.fwd_unit.pos() + fwd_moves[i1] <= chrom.end_pos() - 1);  // NOLINT

    // Extrude rev unit
    lef.rev_unit._pos -= rev_moves[i1];  // Advance extr. unit in 3'-5' direction

    // Extrude fwd unit
    lef.fwd_unit._pos += fwd_moves[i1];                // Advance extr. unit in 5'-3' direction
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

size_t Simulation::register_contacts(Chromosome& chrom, const absl::Span<const Lef> lefs,
                                     const absl::Span<const size_t> selected_lef_idx) const {
  return this->register_contacts(chrom.start_pos() + 1, chrom.end_pos() - 1, chrom.contacts(), lefs,
                                 selected_lef_idx);
}

size_t Simulation::register_contacts(const bp_t start_pos, const bp_t end_pos,
                                     ContactMatrix<contacts_t>& contacts,
                                     const absl::Span<const Lef> lefs,
                                     const absl::Span<const size_t> selected_lef_idx) const {
  // Register contacts for the selected LEFs (excluding LEFs that have one of their units at the
  // beginning/end of a chromosome)
  size_t new_contacts = 0;
  for (const auto i : selected_lef_idx) {
    assert(i < lefs.size());  // NOLINT
    const auto& lef = lefs[i];
    if (BOOST_LIKELY(lef.is_bound() && lef.rev_unit.pos() > start_pos &&
                     lef.rev_unit.pos() < end_pos && lef.fwd_unit.pos() > start_pos &&
                     lef.fwd_unit.pos() < end_pos)) {
      const auto pos1 = lef.rev_unit.pos() - start_pos;
      const auto pos2 = lef.fwd_unit.pos() - start_pos;
      contacts.increment(pos1 / this->bin_size, pos2 / this->bin_size);
      ++new_contacts;
    }
  }
  return new_contacts;
}

size_t Simulation::register_contacts_w_randomization(Chromosome& chrom, absl::Span<const Lef> lefs,
                                                     absl::Span<const size_t> selected_lef_idx,
                                                     random::PRNG_t& rand_eng) const {
  return this->register_contacts_w_randomization(chrom.start_pos() + 1, chrom.end_pos() - 1,
                                                 chrom.contacts(), lefs, selected_lef_idx,
                                                 rand_eng);
}

size_t Simulation::register_contacts_w_randomization(bp_t start_pos, bp_t end_pos,
                                                     ContactMatrix<contacts_t>& contacts,
                                                     absl::Span<const Lef> lefs,
                                                     absl::Span<const size_t> selected_lef_idx,
                                                     random::PRNG_t& rand_eng) const {
  auto noise_gen = genextreme_value_distribution<double>{
      this->genextreme_mu, this->genextreme_sigma, this->genextreme_xi};
  const auto start_pos_dbl = static_cast<double>(start_pos);
  const auto end_pos_dbl = static_cast<double>(end_pos);

  size_t new_contacts = 0;
  for (const auto i : selected_lef_idx) {
    assert(i < lefs.size());  // NOLINT
    const auto& lef = lefs[i];
    if (BOOST_LIKELY(lef.is_bound() && lef.rev_unit.pos() > start_pos &&
                     lef.rev_unit.pos() < end_pos && lef.fwd_unit.pos() > start_pos &&
                     lef.fwd_unit.pos() < end_pos)) {
      // We are performing most operations using double to deal with the possibility that the noise
      // generated to compute p1 is larger than the pos of the rev unit
      const auto p1 = static_cast<double>(lef.rev_unit.pos()) - noise_gen(rand_eng);
      const auto p2 = static_cast<double>(lef.fwd_unit.pos()) + noise_gen(rand_eng);

      if (BOOST_UNLIKELY(p1 < start_pos_dbl || p2 < start_pos_dbl || p1 >= end_pos_dbl ||
                         p2 >= end_pos_dbl)) {
        continue;
      }

      const auto pos1 = static_cast<bp_t>(p1) - start_pos;
      const auto pos2 = static_cast<bp_t>(p2) - start_pos;
      contacts.increment(pos1 / this->bin_size, pos2 / this->bin_size);
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

  // Changing ii -> i raises a -Werror=shadow on GCC 7.5
  auto is_lef_bar_collision = [&](const auto ii) { return ii < barriers.size(); };

  for (auto i = 0UL; i < lefs.size(); ++i) {
    const auto& lef = lefs[i];
    if (BOOST_UNLIKELY(!lef.is_bound())) {
      // Inactive LEF
      lef_unloader_affinity[i] = 0.0;
    } else if (BOOST_UNLIKELY(is_lef_bar_collision(rev_collisions[i]) &&
                              is_lef_bar_collision(fwd_collisions[i]))) {
      // LEF is stalled on both sides by two extrusion barriers
      const auto& rev_barrier = barriers[rev_collisions[i]];
      const auto& fwd_barrier = barriers[fwd_collisions[i]];

      if (BOOST_UNLIKELY(rev_barrier.blocking_direction_major() == dna::rev &&
                         fwd_barrier.blocking_direction_major() == dna::fwd)) {
        // LEF is blocked by a pair of convergent extrusion barriers
        lef_unloader_affinity[i] = 1.0 / this->hard_stall_multiplier;
      } else {
        lef_unloader_affinity[i] = 1.0 / this->soft_stall_multiplier;
      }
    } else {
      // LEF is stalled on at least one side but not by an extrusion barrier
      lef_unloader_affinity[i] = 1.0;
    }
  }
}

size_t Simulation::release_lefs(const absl::Span<Lef> lefs,
                                const absl::Span<const ExtrusionBarrier> barriers,
                                const absl::Span<const collision_t> rev_collisions,
                                const absl::Span<const collision_t> fwd_collisions,
                                random::PRNG_t& rand_eng) const noexcept {
  auto compute_lef_unloader_affinity = [&](const auto j) {
    assert(lefs[j].is_bound());  // NOLINT

    auto is_lef_bar_collision = [&](const auto k) { return k < barriers.size(); };

    const auto& rev_collision = rev_collisions[j];
    const auto& fwd_collision = fwd_collisions[j];
    if (BOOST_LIKELY(!is_lef_bar_collision(rev_collision) ||
                     !is_lef_bar_collision(fwd_collision))) {
      return 1.0 / this->soft_stall_multiplier;
    }

    assert(is_lef_bar_collision(rev_collision));  // NOLINT
    assert(is_lef_bar_collision(fwd_collision));  // NOLINT
    const auto& rev_barrier = barriers[rev_collision];
    const auto& fwd_barrier = barriers[fwd_collision];

    if (BOOST_UNLIKELY(rev_barrier.blocking_direction_major() == dna::rev &&
                       fwd_barrier.blocking_direction_major() == dna::fwd)) {
      return 1.0 / this->hard_stall_multiplier;
    }
    return 1.0;
  };

  size_t lefs_released = 0;
  for (size_t i = 0; i < lefs.size(); ++i) {
    if (BOOST_LIKELY(lefs[i].is_bound())) {
      const auto p = compute_lef_unloader_affinity(i) * this->prob_of_lef_release;
      assert(p >= 0 && p <= 1);  // NOLINT
      if (BOOST_UNLIKELY(random::bernoulli_trial{p}(rand_eng))) {
        ++lefs_released;
        lefs[i].release();
      }
    }
  }
  return lefs_released;
}

thread_pool Simulation::instantiate_thread_pool() const { return thread_pool(this->nthreads); }

Simulation::Task Simulation::Task::from_string(std::string_view serialized_task, Genome& genome) {
  [[maybe_unused]] const size_t ntoks = 7;
  const auto sep = '\t';
  const std::vector<std::string_view> toks = absl::StrSplit(serialized_task, sep);
  assert(toks.size() == ntoks);  // NOLINT

  Task t;
  try {
    size_t i = 0;
    utils::parse_numeric_or_throw(toks[i++], t.id);
    const auto& chrom_name = toks[i++];
    utils::parse_numeric_or_throw(toks[i++], t.cell_id);
    utils::parse_numeric_or_throw(toks[i++], t.num_target_epochs);
    utils::parse_numeric_or_throw(toks[i++], t.num_target_contacts);
    utils::parse_numeric_or_throw(toks[i++], t.num_lefs);
    size_t num_barriers_expected;  // NOLINT
    utils::parse_numeric_or_throw(toks[i++], num_barriers_expected);

    auto chrom_it = genome.find(chrom_name);
    if (chrom_it == genome.end()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Unable to find a chromosome named \"{}\""), chrom_name));
    }
    t.chrom = &(*chrom_it);

    if (t.chrom->num_barriers() != num_barriers_expected) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected {} extrusion barriers for chromosome \"{}\". Found {}"),
                      num_barriers_expected, t.chrom->name(), t.chrom->num_barriers()));
    }
    t.barriers = absl::MakeConstSpan(t.chrom->barriers().data());

  } catch (const std::runtime_error& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occourred while parsing the following task definition:\n"
                               "Task string: \"{}\"\n"
                               "Error reason: {}"),
                    serialized_task, e.what()));
  }
  return t;
}

Simulation::TaskPW Simulation::TaskPW::from_string(std::string_view serialized_task,
                                                   Genome& genome) {
  [[maybe_unused]] const size_t ntoks = 15;
  const auto sep = '\t';
  const std::vector<std::string_view> toks = absl::StrSplit(serialized_task, sep);
  assert(toks.size() == ntoks);  // NOLINT

  TaskPW t;
  try {
    size_t i = 0;
    utils::parse_numeric_or_throw(toks[i++], t.id);
    const auto& chrom_name = toks[i++];
    utils::parse_numeric_or_throw(toks[i++], t.cell_id);
    utils::parse_numeric_or_throw(toks[i++], t.num_target_epochs);
    utils::parse_numeric_or_throw(toks[i++], t.num_target_contacts);
    utils::parse_numeric_or_throw(toks[i++], t.num_lefs);
    size_t num_barriers_expected;  // NOLINT
    utils::parse_numeric_or_throw(toks[i++], num_barriers_expected);
    utils::parse_numeric_or_throw(toks[i++], t.deletion_begin);
    utils::parse_numeric_or_throw(toks[i++], t.deletion_size);
    utils::parse_numeric_or_throw(toks[i++], t.window_start);
    utils::parse_numeric_or_throw(toks[i++], t.window_end);
    utils::parse_numeric_or_throw(toks[i++], t.active_window_start);
    utils::parse_numeric_or_throw(toks[i++], t.active_window_end);
    size_t num_feats1_expected;  // NOLINT
    size_t num_feats2_expected;  // NOLINT
    utils::parse_numeric_or_throw(toks[i++], num_feats1_expected);
    utils::parse_numeric_or_throw(toks[i++], num_feats2_expected);

    auto chrom_it = genome.find(chrom_name);
    if (chrom_it == genome.end()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Unable to find a chromosome named \"{}\""), chrom_name));
    }
    t.chrom = &(*chrom_it);

    if (t.chrom->num_barriers() < num_barriers_expected) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected {} extrusion barriers for chromosome \"{}\". Found {}"),
                      num_barriers_expected, t.chrom->name(), t.chrom->num_barriers()));
    }
    t.barriers = absl::MakeConstSpan(t.chrom->barriers().data());

    assert(t.chrom->get_features().size() == 2);  // NOLINT
    Simulation::map_barriers_to_window(t, *t.chrom);
    Simulation::map_features_to_window(t, *t.chrom);

    if (t.feats1.size() != num_feats1_expected || t.feats2.size() != num_feats2_expected) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Expected {} and {} features for chromosome \"{}\". Found {} and {}"),
          num_feats1_expected, num_feats2_expected, t.chrom->name(), t.feats1.size(),
          t.feats2.size()));
    }

  } catch (const std::runtime_error& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occourred while parsing task definition #{}.\n"
                               "Task definition: \"{}\"\n"
                               "Reason: {}"),
                    toks.front(), serialized_task, e.what()));
  }
  return t;
}

void Simulation::State::resize_buffers(size_t new_size) {
  if (new_size == (std::numeric_limits<size_t>::max)()) {
    new_size = this->num_lefs;
  }
  lef_buff.resize(new_size);
  rank_buff1.resize(new_size);
  rank_buff2.resize(new_size);
  moves_buff1.resize(new_size);
  moves_buff2.resize(new_size);
  idx_buff.resize(new_size);
  collision_buff1.resize(new_size);
  collision_buff2.resize(new_size);

  this->get_barrier_mask().resize(this->barriers.size());
}

void Simulation::State::reset_buffers() {  // TODO figure out which resets are redundant
  std::for_each(lef_buff.begin(), lef_buff.end(), [](auto& lef) { lef.reset(); });
  std::iota(rank_buff1.begin(), rank_buff1.end(), 0);
  std::copy(rank_buff1.begin(), rank_buff1.end(), rank_buff2.begin());
  barrier_mask.reset();
  std::fill(moves_buff1.begin(), moves_buff1.end(), 0);
  std::fill(moves_buff2.begin(), moves_buff2.end(), 0);
  std::fill(collision_buff1.begin(), collision_buff1.end(), 0);
  std::fill(collision_buff2.begin(), collision_buff2.end(), 0);
  barrier_tmp_buff.clear();
  cfx_of_variation_buff.clear();
  avg_loop_size_buff.clear();
}

absl::Span<Lef> Simulation::State::get_lefs(size_t size) noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  assert(size <= this->lef_buff.size());  // NOLINT
  return absl::MakeSpan(this->lef_buff.data(), size);
}
absl::Span<size_t> Simulation::State::get_rev_ranks(size_t size) noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->rank_buff1.data(), size);
}
absl::Span<size_t> Simulation::State::get_fwd_ranks(size_t size) noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->rank_buff2.data(), size);
}
boost::dynamic_bitset<>& Simulation::State::get_barrier_mask() noexcept {
  return this->barrier_mask;
}
absl::Span<bp_t> Simulation::State::get_rev_moves(size_t size) noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->moves_buff1.data(), size);
}
absl::Span<bp_t> Simulation::State::get_fwd_moves(size_t size) noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->moves_buff2.data(), size);
}
absl::Span<size_t> Simulation::State::get_idx_buff(size_t size) noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->idx_buff.data(), size);
}
absl::Span<collision_t> Simulation::State::get_rev_collisions(size_t size) noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->collision_buff1.data(), size);
}
absl::Span<collision_t> Simulation::State::get_fwd_collisions(size_t size) noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->collision_buff2.data(), size);
}
std::deque<double>& Simulation::State::get_cfx_of_variation() noexcept {
  return this->cfx_of_variation_buff;
}
std::deque<double>& Simulation::State::get_avg_loop_sizes() noexcept {
  return this->avg_loop_size_buff;
}

absl::Span<const Lef> Simulation::State::get_lefs(size_t size) const noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->lef_buff.data(), size);
}

absl::Span<const size_t> Simulation::State::get_rev_ranks(size_t size) const noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->rank_buff1.data(), size);
}
absl::Span<const size_t> Simulation::State::get_fwd_ranks(size_t size) const noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->rank_buff2.data(), size);
}
const boost::dynamic_bitset<>& Simulation::State::get_barrier_mask() const noexcept {
  return this->barrier_mask;
}
absl::Span<const bp_t> Simulation::State::get_rev_moves(size_t size) const noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->moves_buff1.data(), size);
}
absl::Span<const bp_t> Simulation::State::get_fwd_moves(size_t size) const noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->moves_buff2.data(), size);
}
absl::Span<const size_t> Simulation::State::get_idx_buff(size_t size) const noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->idx_buff.data(), size);
}
absl::Span<const collision_t> Simulation::State::get_rev_collisions(size_t size) const noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->collision_buff1.data(), size);
}
absl::Span<const collision_t> Simulation::State::get_fwd_collisions(size_t size) const noexcept {
  if (size == 0) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->collision_buff2.data(), size);
}
const std::deque<double>& Simulation::State::get_cfx_of_variation() const noexcept {
  return this->cfx_of_variation_buff;
}
const std::deque<double>& Simulation::State::get_avg_loop_sizes() const noexcept {
  return this->avg_loop_size_buff;
}

bool Simulation::State::is_modle_pert_state() const noexcept { return !this->is_modle_sim_state(); }

bool Simulation::State::is_modle_sim_state() const noexcept {
  return feats1.empty() && feats2.empty();
}

Simulation::State& Simulation::State::operator=(const Task& task) {
  this->epoch = 0;
  this->num_burnin_epochs = 0;
  this->burnin_completed = false;
  this->num_active_lefs = 0;
  this->num_contacts = 0;

  this->id = task.id;
  this->chrom = task.chrom;
  this->cell_id = task.cell_id;
  this->num_target_epochs = task.num_target_epochs;
  this->num_target_contacts = task.num_target_contacts;
  this->num_lefs = task.num_lefs;
  this->barriers = task.barriers;

  this->deletion_begin = 0;
  this->deletion_size = 0;
  this->window_start = 0;
  this->window_end = 0;
  this->active_window_start = 0;
  this->active_window_end = 0;

  this->feats1 = absl::Span<const bed::BED>{};
  this->feats2 = absl::Span<const bed::BED>{};

  if (this->chrom->contacts_ptr()) {
    this->contacts = this->chrom->contacts_ptr();
  }
  return *this;
}

Simulation::State& Simulation::State::operator=(const TaskPW& task) {
  this->epoch = 0;
  this->num_burnin_epochs = 0;
  this->burnin_completed = false;
  this->num_active_lefs = 0;
  this->num_contacts = 0;

  this->id = task.id;
  this->chrom = task.chrom;
  this->cell_id = task.cell_id;
  this->num_target_epochs = task.num_target_epochs;
  this->num_target_contacts = task.num_target_contacts;
  this->num_lefs = task.num_lefs;
  this->barriers = task.barriers;

  this->deletion_begin = task.deletion_begin;
  this->deletion_size = task.deletion_size;
  this->window_start = task.window_start;
  this->window_end = task.window_end;
  this->active_window_start = task.active_window_start;
  this->active_window_end = task.active_window_end;

  this->feats1 = task.feats1;
  this->feats2 = task.feats2;
  return *this;
}

std::string Simulation::State::to_string() const noexcept {
  return fmt::format(FMT_STRING("StatePW:\n - TaskID {}\n"
                                " - Chrom: {}[{}-{}]\n"
                                " - Range start: {}\n"
                                " - Range end: {}\n"
                                " - CellID: {}\n"
                                " - Target epochs: {}\n"
                                " - Target contacts: {}\n"
                                " - # of LEFs: {}\n"
                                " - # Extrusion barriers: {}\n"
                                " - seed: {}\n"),
                     id, chrom->name(), chrom->start_pos(), chrom->end_pos(), window_start,
                     window_end, cell_id, num_target_epochs, num_target_contacts, num_lefs,
                     barriers.size(), seed);
}

std::pair<size_t, size_t> Simulation::process_collisions(
    const Chromosome& chrom, const absl::Span<const Lef> lefs,
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

void Simulation::compute_loop_size_stats(const absl::Span<const Lef> lefs,
                                         std::deque<double>& cfx_of_variation_buff,
                                         std::deque<double>& avg_loop_size_buff,
                                         const size_t buff_capacity) noexcept {
  assert(cfx_of_variation_buff.size() == avg_loop_size_buff.size());  // NOLINT
  if (lefs.empty()) {
    cfx_of_variation_buff.clear();
    avg_loop_size_buff.clear();
    return;
  }

  const auto avg_loop_size =
      math::mean(lefs.begin(), lefs.end(), [](const auto& lef) { return lef.loop_size(); });

  const auto avg_loop_size_std = math::standard_dev(
      lefs.begin(), lefs.end(), avg_loop_size, [](const auto& lef) { return lef.loop_size(); });

  if (avg_loop_size_buff.size() == buff_capacity) {
    avg_loop_size_buff.pop_front();
    cfx_of_variation_buff.pop_front();
  }

  avg_loop_size_buff.push_back(avg_loop_size);
  cfx_of_variation_buff.push_back(avg_loop_size_std / avg_loop_size_buff.back());
}

bool Simulation::evaluate_burnin(const std::deque<double>& cfx_of_variation_buff,
                                 [[maybe_unused]] const std::deque<double>& avg_loop_size_buff,
                                 const size_t buff_capacity, const size_t window_size) noexcept {
  assert(cfx_of_variation_buff.size() == avg_loop_size_buff.size());  // NOLINT
  assert(cfx_of_variation_buff.size() <= buff_capacity);              // NOLINT

  if (cfx_of_variation_buff.size() != buff_capacity) {
    return false;
  }

  // Count the number of adjacent pairs of values where the first value is larger than the
  // second This visually corresponds to a local dip in the avg loop size plot When we are in a
  // stable state the above line should be relatively flat. For this reason, we expect n to be
  // roughly equal to 1/2 of all the measurements
  size_t n = 0;
  assert(window_size < buff_capacity);  // NOLINT
  for (auto it1 = cfx_of_variation_buff.begin() + 1,
            it2 = cfx_of_variation_buff.begin() + static_cast<std::ptrdiff_t>(window_size + 1);
       it2 != cfx_of_variation_buff.end(); ++it1, ++it2) {
    const auto n1 = math::mean(it1 - 1, it2 - 1);
    const auto n2 = math::mean(it1, it2);
    n += static_cast<size_t>(n1 > n2);
  }
  const auto r1 = static_cast<double>(n) / static_cast<double>(buff_capacity - window_size - n);

  const auto cfx_of_variation_is_stable = r1 >= 0.95 && r1 <= 1.05;
  if (!cfx_of_variation_is_stable) {
    return false;
  }

  n = 0;
  for (auto it1 = avg_loop_size_buff.begin() + 1,
            it2 = avg_loop_size_buff.begin() + static_cast<std::ptrdiff_t>(window_size + 1);
       it2 != avg_loop_size_buff.end(); ++it1, ++it2) {
    const auto n1 = math::mean(it1 - 1, it2 - 1);
    const auto n2 = math::mean(it1, it2);
    n += static_cast<size_t>(n1 > n2);
  }
  const auto r2 = static_cast<double>(n) / static_cast<double>(buff_capacity - window_size - n);

  const auto avg_loop_size_is_stable = r2 >= 0.95 && r2 <= 1.05;
  return avg_loop_size_is_stable;
}

void Simulation::run_burnin(State& s, const double lef_binding_rate_burnin) const {
  do {
    ++s.num_burnin_epochs;
    if (s.num_active_lefs != s.num_lefs) {
      const auto num_lefs_to_bind =
          modle::random::poisson_distribution<size_t>{lef_binding_rate_burnin}(s.rand_eng);
      s.num_active_lefs = std::min(s.num_active_lefs + num_lefs_to_bind, s.num_lefs);
    } else {
      Simulation::compute_loop_size_stats(s.get_lefs(), s.get_cfx_of_variation(),
                                          s.get_avg_loop_sizes(), this->burnin_history_length);

      s.burnin_completed = Simulation::evaluate_burnin(
          s.get_cfx_of_variation(), s.get_avg_loop_sizes(), this->burnin_history_length,
          this->burnin_smoothing_window_size);

      if (!s.burnin_completed && s.epoch >= this->max_burnin_epochs) {
        s.burnin_completed = true;
        spdlog::warn(
            FMT_STRING("The burn-in phase for '{}' from cell #{} was terminated earlier than "
                       "it should've as its simulation instance failed to stabilize "
                       "within --max-burnin-epochs={} epochs."),
            s.chrom->name(), s.cell_id, this->max_burnin_epochs);
      }
    }
    // Guard against the rare occasion where the poisson prng samples 0 in the first epoch
  } while (s.num_active_lefs == 0);
}

void Simulation::simulate_one_cell(State& s) const {
  assert(s.epoch == 0);              // NOLINT
  assert(s.num_burnin_epochs == 0);  // NOLINT
  assert(!s.burnin_completed);       // NOLINT
  assert(s.num_active_lefs == 0);    // NOLINT

  // use the current epoch when generating error messages
  try {
    // Seed is computed based on chrom. name, size and cellid
    s.seed = s.chrom->hash(s.xxh_state.get(), this->seed, s.cell_id);
    s.rand_eng = random::PRNG(s.seed);

    const auto lef_binding_rate_burnin =
        static_cast<double>(s.num_lefs) /
        static_cast<double>(this->burnin_target_epochs_for_lef_activation);

    // Compute the avg. # of LEFs to use to sample contact every iterations
    const auto avg_nlefs_to_sample =
        static_cast<double>(s.num_lefs) * this->lef_fraction_contact_sampling;

    // Local counter for the number of contacts generated by the current kernel instance
    s.num_contacts = 0;
    assert((s.num_target_contacts != 0) ==  // NOLINT
           (s.num_target_epochs == (std::numeric_limits<size_t>::max)()));

    // Generate initial extr. barrier states, so that they are already at or close to
    // equilibrium
    for (auto i = 0UL; i < s.barriers.size(); ++i) {
      s.get_barrier_mask()[i] = random::bernoulli_trial{s.barriers[i].occupancy()}(s.rand_eng);
    }

    if (this->skip_burnin) {
      s.num_active_lefs = s.num_lefs;
      s.burnin_completed = true;
    }
    for (; s.epoch < std::max(s.num_burnin_epochs + s.num_target_epochs, s.num_target_epochs);
         ++s.epoch) {
      if (!s.burnin_completed) {
        this->Simulation::run_burnin(s, lef_binding_rate_burnin);
      }

      ////////////////////////
      // Simulate one epoch //
      ////////////////////////

      // Select inactive LEFs and bind them
      Simulation::select_and_bind_lefs(s);
      if (s.burnin_completed) {  // Register contacts
        this->sample_and_register_contacts(s, avg_nlefs_to_sample);
        if (s.num_target_contacts != 0 && s.num_contacts >= s.num_target_contacts) {
          return;  // Enough contact have been generated. Yay!
        }
      }

      this->generate_moves(*s.chrom, s.get_lefs(), s.get_rev_ranks(), s.get_fwd_ranks(),
                           s.get_rev_moves(), s.get_fwd_moves(), s.burnin_completed, s.rand_eng);

      CTCF::update_states(s.barriers, s.get_barrier_mask(), s.rand_eng);

      // Reset collision masks
      std::fill(s.get_rev_collisions().begin(), s.get_rev_collisions().end(), NO_COLLISION);
      std::fill(s.get_fwd_collisions().begin(), s.get_fwd_collisions().end(), NO_COLLISION);

      // Detect collision and correct moves
      const auto [num_rev_units_at_5prime, num_fwd_units_at_3prime] =
          Simulation::process_collisions(*s.chrom, s.get_lefs(), s.barriers, s.get_barrier_mask(),
                                         s.get_rev_ranks(), s.get_fwd_ranks(), s.get_rev_moves(),
                                         s.get_fwd_moves(), s.get_rev_collisions(),
                                         s.get_fwd_collisions(), s.rand_eng);
      // Advance LEFs
      Simulation::extrude(*s.chrom, s.get_lefs(), s.get_rev_moves(), s.get_fwd_moves(),
                          num_rev_units_at_5prime, num_fwd_units_at_3prime);

      // Select LEFs to be released in the current epoch and release them
      this->release_lefs(s.get_lefs(), s.barriers, s.get_rev_collisions(), s.get_fwd_collisions(),
                         s.rand_eng);
    }
  } catch (const std::exception& err) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Caught an exception while simulating epoch #{} in cell #{} of {}:\n{}"),
        s.epoch, s.cell_id, s.chrom->name(), err.what()));
  }
}

void Simulation::select_and_bind_lefs(State& s) noexcept(utils::ndebug_defined()) {
  const auto lef_mask = s.get_idx_buff();
  Simulation::select_lefs_to_bind(s.get_lefs(), lef_mask);
  if (s.is_modle_sim_state()) {
    Simulation::bind_lefs(*s.chrom, s.get_lefs(), s.get_rev_ranks(), s.get_fwd_ranks(), lef_mask,
                          s.rand_eng, s.epoch);
    return;
  }
  Simulation::bind_lefs(s.window_start, s.window_end, s.get_lefs(), s.get_rev_ranks(),
                        s.get_fwd_ranks(), lef_mask, s.rand_eng, s.epoch);
}

void Simulation::sample_and_register_contacts(State& s, const double avg_nlefs_to_sample) const {
  assert(s.num_active_lefs == s.num_lefs);  // NOLINT
  auto nlefs_to_sample =
      std::min(s.num_lefs, random::poisson_distribution<size_t>{avg_nlefs_to_sample}(s.rand_eng));

  if (s.num_target_contacts != 0) {  // When using the target contact density as stopping
    // criterion, don't overshoot the target number of contacts
    nlefs_to_sample = std::min(nlefs_to_sample, s.num_target_contacts - s.num_contacts);
  }

  // Select LEFs to be used for contact registration.
  // We are sampling from fwd ranks to avoid having to allocate a vector of indices just to
  // do this sampling
  const auto lef_idx = s.get_idx_buff(nlefs_to_sample);
  random_sample(s.get_fwd_ranks().begin(), s.get_fwd_ranks().end(), lef_idx.begin(), lef_idx.size(),
                s.rand_eng);
#ifndef NDEBUG  // GCC 9.5.0 chokes if we use if constexpr (utils::ndebug_not_defined()) here
  if (!std::all_of(lef_idx.begin(), lef_idx.end(), [&](const auto i) { return i < s.num_lefs; })) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("lef_idx.size()={}; num_lefs={};\nlef_idx=[{}]\n"), lef_idx.size(),
                    s.num_lefs, fmt::join(lef_idx, ", ")));
  }
#endif
  const auto start_pos = s.is_modle_sim_state() ? s.chrom->start_pos() : s.window_start;
  const auto end_pos = s.is_modle_sim_state() ? s.chrom->end_pos() : s.window_end;

  assert(s.contacts);  // NOLINT
  if (this->randomize_contacts) {
    s.num_contacts += this->register_contacts_w_randomization(
        start_pos + 1, end_pos - 1, *s.contacts, s.get_lefs(), lef_idx, s.rand_eng);
  } else {
    s.num_contacts +=
        this->register_contacts(start_pos + 1, end_pos - 1, *s.contacts, s.get_lefs(), lef_idx);
  }
}
}  // namespace modle
