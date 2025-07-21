// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// IWYU pragma: private, include "modle/simulation.hpp"

#include "modle/simulation.hpp"

#include <absl/time/clock.h>
#include <absl/types/span.h>                    // for Span, MakeConstSpan, MakeSpan
#include <cpp-sort/sorters/insertion_sorter.h>  // for insertion_sort, insertion_so...
#include <cpp-sort/sorters/pdq_sorter.h>        // for pdq_sort, pdq_sorter
#include <cpp-sort/sorters/split_sorter.h>      // for split_sort, split_sorter
#include <fmt/compile.h>
#include <parallel_hashmap/btree.h>  // for btree_iterator
#include <spdlog/spdlog.h>           // for info, warn

#include <algorithm>   // for max, fill, min, copy, clamp
#include <atomic>      // for atomic
#include <cassert>     // for assert
#include <chrono>      // for microseconds
#include <cmath>       // for log, round, exp, floor, sqrt
#include <cstdlib>     // for abs
#include <deque>       // for _Deque_iterator<>::_Self
#include <filesystem>  // for operator<<, path
#include <hictk/cooler.hpp>
#include <limits>       // for numeric_limits
#include <memory>       // for shared_ptr, unique_ptr, make...
#include <mutex>        // for mutex
#include <numeric>      // for iota
#include <stdexcept>    // for runtime_error
#include <string>       // for string
#include <string_view>  // for string_view
#include <utility>      // for make_pair, pair
#include <vector>       // for vector, vector<>::iterator

#include "modle/bigwig/bigwig.hpp"
#include "modle/common/common.hpp"             // for bp_t, contacts_t
#include "modle/common/dna.hpp"                // for dna::REV, dna::FWD
#include "modle/common/random.hpp"             // for bernoulli_trial, poisson_distribution
#include "modle/common/simulation_config.hpp"  // for Config
#include "modle/config/version.hpp"
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier, update_states
#include "modle/extrusion_factors.hpp"   // for Lef, ExtrusionUnit
#include "modle/genome.hpp"              // for Genome::iterator, GenomicInterval
#include "modle/io/contact_matrix_dense.hpp"
#include "modle/stats/descriptive.hpp"

namespace modle {

static void override_extrusion_barrier_occupancy(Genome& genome, const double barrier_occupied_stp,
                                                 const double barrier_not_occupied_stp) {
  for (auto& interval : genome) {
    auto& barriers = interval.barriers();
    std::transform(barriers.begin(), barriers.end(), barriers.begin(), [&](const auto& barrier) {
      return ExtrusionBarrier{barrier.pos, barrier_occupied_stp, barrier_not_occupied_stp,
                              barrier.blocking_direction.complement()};
    });
  }
}

static void issue_warnings_for_small_intervals(const Genome& genome, const bp_t diagonal_width) {
  std::vector<std::string> warnings;
  for (const auto& interval : genome) {
    if (interval.size() < diagonal_width) {
      warnings.emplace_back(fmt::format(FMT_STRING("{}: {};"), interval, interval.size()));
    }
  }
  if (!warnings.empty()) {
    spdlog::warn(FMT_STRING("simulated size for the following {} interval(s) is smaller than "
                            "the simulation diagonal width ({} bp). Is this intended?\n - {}"),
                 warnings.size(), diagonal_width, fmt::join(warnings, "\n - "));
  }
}

template <usize num_pixels_threshold = 250'000>
static void issue_warnings_for_small_matrices(const Genome& genome) {
  std::vector<std::string> warnings;
  for (const auto& interval : genome) {
    if (const auto npixels = interval.npixels(); npixels < num_pixels_threshold) {
      warnings.emplace_back(fmt::format(FMT_STRING("{}: {} pixels"), interval, npixels));
    }
  }
  if (!warnings.empty()) {
    spdlog::warn(FMT_STRING("contact matrix for the following {} interval(s) appears to be "
                            "really small (less than {} pixels). Is this intended?\n - {}"),
                 warnings.size(), num_pixels_threshold, fmt::join(warnings, "\n - "));
  }
}

Simulation::Simulation(Config config_, bool import_intervals)
    : _config(std::move(config_)),
      _genome(import_intervals
                  ? Genome(c().path_to_chrom_sizes, c().path_to_extr_barriers,
                           c().path_to_genomic_intervals, c().bin_size, c().diagonal_width,
                           c().barrier_occupied_stp, c().barrier_not_occupied_stp,
                           c().interpret_bed_name_field_as_barrier_not_occupied_stp)
                  : Genome{}),
      _ctx(Simulation::compute_num_worker_threads(c().nthreads, this->_genome.num_intervals(),
                                                  c().num_cells)) {
  if (c().override_extrusion_barrier_occupancy) {
    // Override barrier occupancies read from BED file
    override_extrusion_barrier_occupancy(this->_genome, c().barrier_occupied_stp,
                                         c().barrier_not_occupied_stp);
  }

  issue_warnings_for_small_intervals(this->_genome, c().diagonal_width);
  issue_warnings_for_small_matrices(this->_genome);
}

usize Simulation::size() const { return this->_genome.size(); }

usize Simulation::simulated_size() const { return this->_genome.simulated_size(); }

[[nodiscard]] static hictk::cooler::File init_cooler_file(
    const std::filesystem::path& path, bool force,
    const std::vector<std::shared_ptr<const Chromosome>>& chroms, bp_t bin_size,
    std::string_view assembly_name, std::string_view metadata) {
  std::vector<Chromosome> chroms_(chroms.size());
  std::transform(chroms.begin(), chroms.end(), chroms_.begin(),
                 [](const auto& chrom_ptr) { return *chrom_ptr; });

  return io::init_cooler_file<i32>(path, force, chroms_.begin(), chroms_.end(), bin_size,
                                   assembly_name, config::version::str_long(), metadata);
}

[[nodiscard]] static io::bigwig::Writer init_bigwig_writer(const std::filesystem::path& path,
                                                           const Genome& genome) {
  auto bwf = io::bigwig::Writer(path.string());
  std::vector<std::pair<std::string, u32>> chroms(genome.num_chromosomes());
  std::transform(genome.chromosomes().begin(), genome.chromosomes().end(), chroms.begin(),
                 [](const auto& chrom_ptr) {
                   return std::make_pair(chrom_ptr->name(), static_cast<u32>(chrom_ptr->size()));
                 });
  bwf.write_chromosomes(chroms);

  return bwf;
}

static void write_contact_matrix_to_cooler(hictk::cooler::File& cf,
                                           const GenomicInterval& interval) {
  if (!cf) {
    return;
  }

  const auto& matrix = interval.contacts();
  if (matrix.npixels() != 0) {
    const auto t0 = absl::Now();
    spdlog::info(FMT_STRING("[io] writing contacts for {} to file \"{}\"..."), interval, cf.uri());
    const auto missing_interactions = matrix.get_fraction_of_missed_updates();
    if (missing_interactions >= 0.01) {
      spdlog::warn(
          FMT_STRING(
              "[io] {:.2f}% missing interactions for {}! Please make sure this is intended."),
          missing_interactions * 100, interval);
    }

    io::append_contact_matrix_to_cooler(cf, interval.chrom().name(), matrix, interval.start());
    spdlog::info(
        FMT_STRING("[io]: written {} contacts for {} to file \"{}\" in {} ({:.2f}M nnz out of "
                   "{:.2f}M pixels)."),
        matrix.get_tot_contacts(), interval, cf.uri(), absl::FormatDuration(absl::Now() - t0),
        static_cast<double>(matrix.get_nnz()) / 1.0e6,
        static_cast<double>(matrix.npixels()) / 1.0e6);
  }
}

static void write_lef_occupancy_to_bwig(io::bigwig::Writer& bw, const GenomicInterval& interval,
                                        const bp_t resolution, std::vector<float>& buff) {
  if (!bw || interval.npixels() == 0) {
    return;
  }

  try {
    const auto t0 = absl::Now();
    spdlog::info(FMT_STRING("[io]: writing 1D LEF occupancy profile for {} to file {}..."),
                 interval, bw.path());
    buff.resize(interval.lef_1d_occupancy().size());
    const auto max_element =
        std::max_element(interval.lef_1d_occupancy().begin(), interval.lef_1d_occupancy().end())
            ->load();

    std::transform(interval.lef_1d_occupancy().begin(), interval.lef_1d_occupancy().end(),
                   buff.begin(), [&](const auto& n) {
                     return static_cast<float>(static_cast<double>(n.load()) /
                                               static_cast<double>(max_element));
                   });
    bw.write_range(interval.chrom().name(), absl::MakeSpan(buff), resolution, resolution,
                   interval.start());
    spdlog::info(FMT_STRING("[io]: writing 1D LEF occupancy profile for {} to file {} took {}"),
                 interval, bw.path(), absl::FormatDuration(absl::Now() - t0));
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("An error occurred while writing LEF 1D occupancy profile to file {}: {}"),
        bw.path(), e.what()));
  }
}

// GCC8 doesn't compile init_output_file_writers if we use std::pair<> as return type
struct CoolerBigwigPair {
  hictk::cooler::File clr{};
  io::bigwig::Writer bw{};
};

[[nodiscard]] static CoolerBigwigPair init_output_file_writers(const Genome& genome,
                                                               const Config& c) {
  if (c.skip_output) {
    return {};
  }

  return {init_cooler_file(c.path_to_output_file_cool, false, genome.chromosomes(), c.bin_size,
                           c.assembly_name, c.args_json),
          c.track_1d_lef_position ? init_bigwig_writer(c.path_to_lef_1d_occupancy_bw_file, genome)
                                  : io::bigwig::Writer{}};
}

void Simulation::simulate_io(std::chrono::milliseconds wait_time) {
  spdlog::info("spawning IO thread...");
  assert(!std::filesystem::exists(c().path_to_output_file_cool));
  if (c().track_1d_lef_position) {
    assert(!std::filesystem::exists(c().path_to_lef_1d_occupancy_bw_file));
  }

  Task task{};

  try {
    auto [cf, bw] = init_output_file_writers(this->_genome, c());
    auto current_interval = this->_genome.begin();
    auto last_interval = this->_genome.end();

    auto ctok = this->_ctx.register_consumer<Task::Status::COMPLETED>();
    phmap::btree_map<GenomicInterval*, usize> task_map{};
    std::vector<float> bw_buff{};
    while (!!this->_ctx && current_interval != last_interval) {
      if (auto it = task_map.find(&(*current_interval)); it != task_map.end() && it->second == 0) {
        // All tasks for the current interval successfully completed!
        write_contact_matrix_to_cooler(cf, *it->first);
        write_lef_occupancy_to_bwig(bw, *it->first, c().bin_size, bw_buff);
        it->first->deallocate();
        task_map.erase(it);
        ++current_interval;
        continue;
      }

      const auto task_opt = this->_ctx.wait_dequeue_task<Task::Status::COMPLETED>(ctok, wait_time);
      if (task_opt.has_value()) {
        task = *task_opt;
        // Add interval to task_map if not already present.
        // Then decrement the number of tasks still pending.
        auto [it, _] = task_map.try_emplace(task.interval, c().num_cells);
        --it->second;
      }
    }
    assert(this->_ctx.num_tasks_submitted() == this->_ctx.num_tasks_completed());
  } catch (const std::exception& e) {
    if (task.interval) {
      this->_ctx.throw_exception(std::runtime_error(fmt::format(
          FMT_STRING("exception encountered while writing interactions for {} to file {}: {}"),
          *task.interval, c().path_to_output_file_cool, e.what())));
    }
    this->_ctx.throw_exception(std::runtime_error(
        fmt::format(FMT_STRING("exception encountered while writing interactions to file {}: {}"),
                    c().path_to_output_file_cool, e.what())));
  } catch (...) {
    this->_ctx.throw_exception(
        std::runtime_error("unhandled exception caught in IO thread! This should never "
                           "happen! Please file an issue on GitHub."));
  }
}

// Generate moves for extrusion units moving in the same direction
template <class MoveGeneratorT>
static void generate_moves_helper(const absl::Span<const Lef> lefs, const absl::Span<bp_t> moves,
                                  const double avg_extr_speed, const double extr_speed_std,
                                  random::PRNG_t& rand_eng) {
  const bp_t move_int = static_cast<bp_t>(std::round(avg_extr_speed));
  auto generate_move = [&](const auto& lef) -> bp_t {
    // Don't bother generating move for inactive LEFs
    if (!lef.is_bound()) {
      return 0;
    }

    // When std == 0 always return the avg. extrusion speed
    if (extr_speed_std == 0.0) {
      return move_int;
    }
    // NOTE: on my laptop generating doubles from a normal distribution, rounding them, then
    // casting double to uint is a lot faster (~4x) than drawing uints directly from a Poisson
    // distr.
    return static_cast<bp_t>(
        std::round(std::max(0.0, MoveGeneratorT{avg_extr_speed, extr_speed_std}(rand_eng))));
  };

  for (usize i = 0; i < lefs.size(); ++i) {
    moves[i] = generate_move(lefs[i]);
  }
}

void Simulation::generate_moves(const GenomicInterval& interval, const absl::Span<const Lef> lefs,
                                const absl::Span<const usize> rev_lef_ranks,
                                const absl::Span<const usize> fwd_lef_ranks,
                                const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
                                const bool burnin_completed, random::PRNG_t& rand_eng,
                                bool adjust_moves_) const noexcept(utils::ndebug_defined()) {
  {
    assert(lefs.size() == fwd_lef_ranks.size());
    assert(lefs.size() == rev_lef_ranks.size());
    assert(lefs.size() == fwd_moves.size());
    assert(lefs.size() == rev_moves.size());
  }

  const auto rev_extr_speed = static_cast<double>(
      burnin_completed ? c().rev_extrusion_speed : c().rev_extrusion_speed_burnin);
  const auto fwd_extr_speed = static_cast<double>(
      burnin_completed ? c().fwd_extrusion_speed : c().fwd_extrusion_speed_burnin);

  generate_moves_helper<lef_move_generator_t>(lefs, rev_moves, rev_extr_speed,
                                              c().rev_extrusion_speed_std, rand_eng);
  generate_moves_helper<lef_move_generator_t>(lefs, fwd_moves, fwd_extr_speed,
                                              c().fwd_extrusion_speed_std, rand_eng);

  if (adjust_moves_) {  // Adjust moves of consecutive extr. units to make LEF behavior more
    // realistic See comments in adjust_moves_of_consecutive_extr_units for more
    // details on what this entails
    Simulation::adjust_moves_of_consecutive_extr_units(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                       rev_moves, fwd_moves);
  }

  Simulation::clamp_moves(interval, lefs, rev_moves, fwd_moves);
}

void Simulation::clamp_moves(const GenomicInterval& interval, const absl::Span<const Lef> lefs,
                             const absl::Span<bp_t> rev_moves,
                             const absl::Span<bp_t> fwd_moves) noexcept {
  for (usize i = 0; i < lefs.size(); ++i) {
    if (!lefs[i].is_bound()) {
      assert(rev_moves[i] == 0);
      assert(fwd_moves[i] == 0);
      continue;
    }

    assert(lefs[i].rev_unit.pos() >= interval.start());
    assert(lefs[i].fwd_unit.pos() < interval.end());
    rev_moves[i] = std::min(rev_moves[i], lefs[i].rev_unit.pos() - interval.start());
    fwd_moves[i] = std::min(fwd_moves[i], interval.end() - lefs[i].fwd_unit.pos() - 1);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Simulation::adjust_moves_of_consecutive_extr_units(
    [[maybe_unused]] const GenomicInterval& interval, absl::Span<const Lef> lefs,
    absl::Span<const usize> rev_lef_ranks, absl::Span<const usize> fwd_lef_ranks,
    absl::Span<bp_t> rev_moves, absl::Span<bp_t> fwd_moves) noexcept(utils::ndebug_defined()) {
  assert(!lefs.empty());

  // Loop over pairs of consecutive extr. units.
  // Extr. units moving in rev direction are processed in 3'-5' order, while units moving in fwd
  // direction are processed in 5'-3' direction
  for (usize i = lefs.size() - 1; i > 0; --i) {
    const auto i1 = rev_lef_ranks[i - 1];
    const auto i2 = rev_lef_ranks[i];

    if (lefs[i1].is_bound() && lefs[i2].is_bound()) {
      // Don't bother processing moves that will cause LEFs to go past chromosomal boundaries: we
      // have a post-processing pass at the end of generate_moves() to handle this case
      if (lefs[i1].rev_unit.pos() <= interval.start() + rev_moves[i1] ||
          lefs[i2].rev_unit.pos() <= interval.start() + rev_moves[i2]) {
        continue;
      }

      const auto pos1 = lefs[i1].rev_unit.pos() - rev_moves[i1];
      const auto pos2 = lefs[i2].rev_unit.pos() - rev_moves[i2];

      // If moving extr. units 1 and 2 by their respective moves would cause unit 1 (which comes
      // first in 3'-5' order) to be surpassed by unit 2 (which comes second in 3'-5' order),
      // increment the move for unit 1 so that after moving both units, unit 1 will be located one
      // bp upstream of unit 2 in 5'-3' direction.
      // This mimics what would probably happen in a real system, where extr. unit 2 would most
      // likely push extr. unit 1, temporarily increasing extr. speed of unit 1.
      if (pos2 <= pos1) {
        rev_moves[i1] += (pos1 - pos2) + 1;
      }
    }
  }

  // Process fwd units
  for (usize i = 1; i < lefs.size(); ++i) {
    const auto i1 = fwd_lef_ranks[i - 1];
    const auto i2 = fwd_lef_ranks[i];

    // See above for detailed comments. The code logic is the same used on rev units (but
    // mirrored!)
    if (lefs[i1].is_bound() && lefs[i2].is_bound()) {
      if (lefs[i1].fwd_unit.pos() + fwd_moves[i1] > interval.end() - 1 ||
          lefs[i2].fwd_unit.pos() + fwd_moves[i2] > interval.end() - 1) {
        continue;
      }

      const auto pos1 = lefs[i1].fwd_unit.pos() + fwd_moves[i1];
      const auto pos2 = lefs[i2].fwd_unit.pos() + fwd_moves[i2];

      if (pos1 >= pos2) {
        fwd_moves[i2] += (pos1 - pos2) + 1;
      }
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Simulation::rank_lefs(const absl::Span<const Lef> lefs,
                           const absl::Span<usize> rev_lef_rank_buff,
                           const absl::Span<usize> fwd_lef_rank_buff,
                           bool ranks_are_partially_sorted,
                           bool init_buffers) noexcept(utils::ndebug_defined()) {
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());

  auto rev_comparator = [&](const auto r1, const auto r2) constexpr noexcept {
    assert(r1 < lefs.size());
    assert(r2 < lefs.size());
    return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
  };

  auto fwd_comparator = [&](const auto r1, const auto r2) constexpr noexcept {
    assert(r1 < lefs.size());
    assert(r2 < lefs.size());
    return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
  };

  if (MODLE_UNLIKELY(init_buffers)) {  // Init rank buffers
    std::iota(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), 0);
    std::iota(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), 0);
  }

  if (MODLE_LIKELY(ranks_are_partially_sorted)) {
    cppsort::split_adapter<cppsort::pdq_sorter> sorter;
    sorter(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), rev_comparator);
    sorter(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), fwd_comparator);
  } else {
    // Fallback to pattern-defeating quicksort we have no information regarding the level of
    // pre-sortedness of LEFs
    cppsort::pdq_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), rev_comparator);
    cppsort::pdq_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), fwd_comparator);
  }

  // TODO Figure out a better way to deal with ties
  usize begin{};
  usize end{};
  for (usize i = 1; i < rev_lef_rank_buff.size(); ++i) {
    const auto& r1 = rev_lef_rank_buff[i - 1];
    const auto& r2 = rev_lef_rank_buff[i];
    if (MODLE_UNLIKELY(lefs[r1].rev_unit.pos() == lefs[r2].rev_unit.pos())) {
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
                                assert(r11 < lefs.size());
                                assert(r22 < lefs.size());
                                return lefs[r11].binding_epoch < lefs[r22].binding_epoch;
                              });
      // begin = end;
    }
  }

  for (usize i = 1; i < fwd_lef_rank_buff.size(); ++i) {
    const auto& r1 = fwd_lef_rank_buff[i - 1];
    const auto& r2 = fwd_lef_rank_buff[i];
    if (MODLE_UNLIKELY(lefs[r1].fwd_unit.pos() == lefs[r2].fwd_unit.pos())) {
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
                                assert(r11 < lefs.size());
                                assert(r22 < lefs.size());
                                return lefs[r22].binding_epoch < lefs[r11].binding_epoch;
                              });
      // begin = end;
    }
  }
}

void Simulation::extrude([[maybe_unused]] const GenomicInterval& interval,
                         const absl::Span<Lef> lefs, const absl::Span<const bp_t> rev_moves,
                         const absl::Span<const bp_t> fwd_moves) noexcept(utils::ndebug_defined()) {
  assert(lefs.size() == rev_moves.size());
  assert(lefs.size() == fwd_moves.size());

  for (usize i = 0; i < lefs.size(); ++i) {
    auto& lef = lefs[i];
    const auto& rev_move = rev_moves[i];
    const auto& fwd_move = fwd_moves[i];
    if (MODLE_UNLIKELY(!lef.is_bound())) {  // Do not process inactive LEFs
      continue;
    }
    assert(lef.rev_unit.pos() >= interval.start() + rev_move);
    assert(lef.fwd_unit.pos() + fwd_move < interval.end());

    // Extrude rev unit
    lef.rev_unit._pos -= rev_move;  // Advance extr. unit in 3'-5' direction

    // Extrude fwd unit
    lef.fwd_unit._pos += fwd_move;  // Advance extr. unit in 5'-3' direction
    assert(lef.rev_unit.pos() <= lef.fwd_unit.pos());
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
  assert(collision_pos <= rev_pos);
  if constexpr (utils::ndebug_not_defined()) {
    [[maybe_unused]] const auto collision_pos_ =
        static_cast<double>(rev_pos) - static_cast<double>(rev_speed) * time_to_collision;
    assert(abs(static_cast<double>(collision_pos) - collision_pos_) < 1.0);
  }
  if (MODLE_UNLIKELY(collision_pos == fwd_pos)) {
    assert(collision_pos >= fwd_pos);
    assert(collision_pos + 1 <= rev_pos);
    return std::make_pair(collision_pos + 1, collision_pos);
  }
  assert(collision_pos > 0);
  assert(collision_pos - 1 >= fwd_pos);
  return std::make_pair(collision_pos, collision_pos - 1);
}

usize Simulation::release_lefs(const absl::Span<Lef> lefs, const ExtrusionBarriers& barriers,
                               const absl::Span<const CollisionT> rev_collisions,
                               const absl::Span<const CollisionT> fwd_collisions,
                               random::PRNG_t& rand_eng,
                               const bool burnin_completed) const noexcept {
  auto compute_lef_unloader_affinity = [&](const auto j) {
    assert(lefs[j].is_bound());

    auto is_lef_bar_hard_collision = [&](const auto c, dna::Direction d) {
      assert(d == dna::REV || d == dna::FWD);
      if (c.collision_occurred(CollisionT::LEF_BAR)) {
        assert(c.decode_index() <= barriers.size());
        return barriers.direction(c.decode_index()) == d;
      }
      return false;
    };

    const auto num_hard_collisions =
        static_cast<u8>(is_lef_bar_hard_collision(rev_collisions[j], dna::REV) +
                        is_lef_bar_hard_collision(fwd_collisions[j], dna::FWD));

    switch (num_hard_collisions) {
      case 0:
        return 1.0;
      case 1:
        return 1.0 / c().soft_stall_lef_stability_multiplier;
      case 2:
        return 1.0 / c().hard_stall_lef_stability_multiplier;
      default:
        MODLE_UNREACHABLE_CODE;
    }
  };

  const auto& base_prob_lef_release =
      burnin_completed ? c().prob_of_lef_release : c().prob_of_lef_release_burnin;

  usize lefs_released = 0;
  for (usize i = 0; i < lefs.size(); ++i) {
    if (MODLE_LIKELY(lefs[i].is_bound())) {
      const auto p = compute_lef_unloader_affinity(i) * base_prob_lef_release;
      assert(p >= 0 && p <= 1);
      if (MODLE_UNLIKELY(random::bernoulli_trial{p}(rand_eng))) {
        ++lefs_released;
        lefs[i].release();
      }
    }
  }
  return lefs_released;
}

void Simulation::State::resize_buffers(usize new_size) {
  if (new_size == (std::numeric_limits<usize>::max)()) {
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
}

void Simulation::State::reset_buffers() {  // TODO figure out which resets are redundant
  std::for_each(lef_buff.begin(), lef_buff.end(), [](auto& lef) { lef.reset(); });
  std::iota(rank_buff1.begin(), rank_buff1.end(), 0);
  std::copy(rank_buff1.begin(), rank_buff1.end(), rank_buff2.begin());
  std::fill(moves_buff1.begin(), moves_buff1.end(), 0);
  std::fill(moves_buff2.begin(), moves_buff2.end(), 0);
  std::for_each(collision_buff1.begin(), collision_buff1.end(), [&](auto& c) { c.clear(); });
  std::for_each(collision_buff2.begin(), collision_buff2.end(), [&](auto& c) { c.clear(); });
  cfx_of_variation_buff.clear();
  avg_loop_size_buff.clear();
}

absl::Span<Lef> Simulation::State::get_lefs(usize size) noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  assert(size <= this->lef_buff.size());
  return absl::MakeSpan(this->lef_buff.data(), size);
}
absl::Span<usize> Simulation::State::get_rev_ranks(usize size) noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->rank_buff1.data(), size);
}
absl::Span<usize> Simulation::State::get_fwd_ranks(usize size) noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->rank_buff2.data(), size);
}
absl::Span<bp_t> Simulation::State::get_rev_moves(usize size) noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->moves_buff1.data(), size);
}
absl::Span<bp_t> Simulation::State::get_fwd_moves(usize size) noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->moves_buff2.data(), size);
}
absl::Span<usize> Simulation::State::get_idx_buff(usize size) noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->idx_buff.data(), size);
}
auto Simulation::State::get_rev_collisions(usize size) noexcept -> absl::Span<CollisionT> {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeSpan(this->collision_buff1.data(), size);
}
auto Simulation::State::get_fwd_collisions(usize size) noexcept -> absl::Span<CollisionT> {
  if (size == State::npos) {
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

absl::Span<const Lef> Simulation::State::get_lefs(usize size) const noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->lef_buff.data(), size);
}

absl::Span<const usize> Simulation::State::get_rev_ranks(usize size) const noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->rank_buff1.data(), size);
}
absl::Span<const usize> Simulation::State::get_fwd_ranks(usize size) const noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->rank_buff2.data(), size);
}
absl::Span<const bp_t> Simulation::State::get_rev_moves(usize size) const noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->moves_buff1.data(), size);
}
absl::Span<const bp_t> Simulation::State::get_fwd_moves(usize size) const noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->moves_buff2.data(), size);
}
absl::Span<const usize> Simulation::State::get_idx_buff(usize size) const noexcept {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->idx_buff.data(), size);
}
auto Simulation::State::get_rev_collisions(usize size) const noexcept
    -> absl::Span<const CollisionT> {
  if (size == State::npos) {
    size = this->num_active_lefs;
  }
  return absl::MakeConstSpan(this->collision_buff1.data(), size);
}
auto Simulation::State::get_fwd_collisions(usize size) const noexcept
    -> absl::Span<const CollisionT> {
  if (size == State::npos) {
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

Simulation::State& Simulation::State::operator=(const Task& task) {
  this->epoch = 0;
  this->num_burnin_epochs = 0;
  this->burnin_completed = false;
  this->num_active_lefs = 0;
  this->num_contacts = 0;

  this->id = task.id;
  this->interval = task.interval;
  this->cell_id = task.cell_id;
  this->num_target_epochs = task.num_target_epochs;
  this->num_target_contacts = task.num_target_contacts;
  this->num_lefs = task.num_lefs;
  if (task.interval) {
    this->barriers =
        ExtrusionBarriers{task.interval->barriers().begin(), task.interval->barriers().end()};
  }
  this->rand_eng = task.rand_eng;

  return *this;
}

std::pair<usize, usize> Simulation::process_collisions(
    const GenomicInterval& interval, const absl::Span<Lef> lefs, ExtrusionBarriers& barriers,
    const absl::Span<usize> rev_lef_ranks, const absl::Span<usize> fwd_lef_ranks,
    const absl::Span<bp_t> rev_moves, const absl::Span<bp_t> fwd_moves,
    const absl::Span<CollisionT> rev_collisions, const absl::Span<CollisionT> fwd_collisions,
    random::PRNG_t& rand_eng) const noexcept(utils::ndebug_defined()) {
  const auto& [num_rev_units_at_5prime, num_fwd_units_at_3prime] =
      Simulation::detect_units_at_interval_boundaries(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                                      rev_moves, fwd_moves, rev_collisions,
                                                      fwd_collisions);

  this->detect_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves,
                                  barriers, rev_collisions, fwd_collisions, rand_eng,
                                  num_rev_units_at_5prime, num_fwd_units_at_3prime);

  this->detect_primary_lef_lef_collisions(lefs, barriers, rev_lef_ranks, fwd_lef_ranks, rev_moves,
                                          fwd_moves, rev_collisions, fwd_collisions, rand_eng,
                                          num_rev_units_at_5prime, num_fwd_units_at_3prime);
  Simulation::correct_moves_for_lef_bar_collisions(lefs, barriers, rev_moves, fwd_moves,
                                                   rev_collisions, fwd_collisions);

  Simulation::correct_moves_for_primary_lef_lef_collisions(
      lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, rev_collisions, fwd_collisions);
  this->process_secondary_lef_lef_collisions(
      interval, lefs, rev_lef_ranks, fwd_lef_ranks, rev_moves, fwd_moves, rev_collisions,
      fwd_collisions, rand_eng, num_rev_units_at_5prime, num_fwd_units_at_3prime);
  Simulation::fix_secondary_lef_lef_collisions(interval, lefs, rev_lef_ranks, fwd_lef_ranks,
                                               rev_moves, fwd_moves, rev_collisions, fwd_collisions,
                                               num_rev_units_at_5prime, num_fwd_units_at_3prime);
  return std::make_pair(num_rev_units_at_5prime, num_fwd_units_at_3prime);
}

void Simulation::compute_loop_size_stats(const absl::Span<const Lef> lefs,
                                         std::deque<double>& cfx_of_variation_buff,
                                         std::deque<double>& avg_loop_size_buff,
                                         const usize buff_capacity) noexcept {
  assert(cfx_of_variation_buff.size() == avg_loop_size_buff.size());
  if (lefs.empty()) {
    cfx_of_variation_buff.clear();
    avg_loop_size_buff.clear();
    return;
  }

  const auto avg_loop_size =
      stats::mean(lefs.begin(), lefs.end(), [](const auto& lef) { return lef.loop_size(); });

  const auto avg_loop_size_std = stats::standard_dev(
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
                                 const usize buff_capacity, const usize window_size) noexcept {
  assert(cfx_of_variation_buff.size() == avg_loop_size_buff.size());
  assert(cfx_of_variation_buff.size() <= buff_capacity);

  if (cfx_of_variation_buff.size() != buff_capacity) {
    return false;
  }

  // Count the number of adjacent pairs of values where the first value is larger than the second.
  // This visually corresponds to a local dip in the avg loop size plot.
  // When we are in a stable state the above line should be relatively flat.
  // For this reason, we expect n to be roughly equal to 1/2 of all the measurements
  usize n = 0;
  assert(window_size < buff_capacity);
  for (auto it1 = cfx_of_variation_buff.begin() + 1,
            it2 = cfx_of_variation_buff.begin() + static_cast<isize>(window_size + 1);
       it2 != cfx_of_variation_buff.end(); ++it1, ++it2) {
    const auto n1 = stats::mean(it1 - 1, it2 - 1);
    const auto n2 = stats::mean(it1, it2);
    n += static_cast<usize>(n1 > n2);
  }
  const auto r1 = static_cast<double>(n) / static_cast<double>(buff_capacity - window_size - n);

  const auto cfx_of_variation_is_stable = r1 >= 0.95 && r1 <= 1.05;
  if (!cfx_of_variation_is_stable) {
    return false;
  }

  n = 0;
  for (auto it1 = avg_loop_size_buff.begin() + 1,
            it2 = avg_loop_size_buff.begin() + static_cast<isize>(window_size + 1);
       it2 != avg_loop_size_buff.end(); ++it1, ++it2) {
    const auto n1 = stats::mean(it1 - 1, it2 - 1);
    const auto n2 = stats::mean(it1, it2);
    n += static_cast<usize>(n1 > n2);
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
          modle::random::poisson_distribution<usize>{lef_binding_rate_burnin}(s.rand_eng);
      s.num_active_lefs = std::min(s.num_active_lefs + num_lefs_to_bind, s.num_lefs);
    } else {
      Simulation::compute_loop_size_stats(s.get_lefs(), s.get_cfx_of_variation(),
                                          s.get_avg_loop_sizes(), c().burnin_history_length);

      s.burnin_completed =
          Simulation::evaluate_burnin(s.get_cfx_of_variation(), s.get_avg_loop_sizes(),
                                      c().burnin_history_length, c().burnin_smoothing_window_size);
      s.burnin_completed &= s.epoch > c().min_burnin_epochs;

      if (!s.burnin_completed && s.epoch >= c().max_burnin_epochs) {
        s.burnin_completed = true;
        s.num_active_lefs = s.num_lefs;
        spdlog::warn(
            FMT_STRING("The burn-in phase for {} from cell #{} was terminated earlier than "
                       "it should've as its simulation instance failed to stabilize "
                       "within --max-burnin-epochs={} epochs."),
            *s.interval, s.cell_id, c().max_burnin_epochs);
      }
    }
    // Guard against the rare occasion where the poisson prng samples 0 in the first epoch
  } while (s.num_active_lefs == 0);
}

void Simulation::simulate_one_cell([[maybe_unused]] u64 tid, State& s) const {
  assert(s.epoch == 0);
  assert(s.num_burnin_epochs == 0);
  assert(!s.burnin_completed);
  assert(s.num_active_lefs == 0);
  assert(s.num_target_epochs != (std::numeric_limits<usize>::max)());

  try {
    // Seed is computed based on chrom. name, interval start/end and cellid

    const auto lef_binding_rate_burnin =
        static_cast<double>(s.num_lefs) /
        static_cast<double>(c().burnin_target_epochs_for_lef_activation);

    const auto sampling_events_per_epoch = this->compute_contacts_per_epoch(s.num_lefs);

    // Init the local counter for the number of contacts generated by the current simulation
    // instance
    s.num_contacts = 0;

    // Generate initial extr. barrier states, so that they are already at or close to
    // equilibrium
    s.barriers.init_states(s.rand_eng);

    if (c().skip_burnin) {
      s.num_active_lefs = s.num_lefs;
      s.burnin_completed = true;
    }

    auto stop_condition = [this, &s]() {
      if (c().target_contact_density >= 0) {
        assert(s.num_target_contacts != 0);
        return s.num_contacts >= s.num_target_contacts;
      }
      return s.epoch - s.num_burnin_epochs >= s.num_target_epochs;
    };

    for (; this->_ctx && !stop_condition(); ++s.epoch) {
      if (!s.burnin_completed) {
        this->Simulation::run_burnin(s, lef_binding_rate_burnin);
      }

      ////////////////////////
      // Simulate one epoch //
      ////////////////////////

      // Select inactive LEFs and bind them
      Simulation::select_and_bind_lefs(s);
      if (s.burnin_completed) {  // Register contacts
        this->sample_and_register_contacts(s, sampling_events_per_epoch);
        if (s.num_target_contacts != 0 && s.num_contacts >= s.num_target_contacts) {
          return;  // Enough contacts have been generated. Yay!
        }
      }

      this->generate_moves(*s.interval, s.get_lefs(), s.get_rev_ranks(), s.get_fwd_ranks(),
                           s.get_rev_moves(), s.get_fwd_moves(), s.burnin_completed, s.rand_eng);

      s.barriers.next_state(s.rand_eng);

      // Reset collision masks
      std::for_each(s.get_rev_collisions().begin(), s.get_rev_collisions().end(),
                    [&](auto& c) { c.clear(); });
      std::for_each(s.get_fwd_collisions().begin(), s.get_fwd_collisions().end(),
                    [&](auto& c) { c.clear(); });

      // Detect collision and correct moves
      Simulation::process_collisions(*s.interval, s.get_lefs(), s.barriers, s.get_rev_ranks(),
                                     s.get_fwd_ranks(), s.get_rev_moves(), s.get_fwd_moves(),
                                     s.get_rev_collisions(), s.get_fwd_collisions(), s.rand_eng);
      // Advance LEFs
      Simulation::extrude(*s.interval, s.get_lefs(), s.get_rev_moves(), s.get_fwd_moves());

      // Log model internal state
      if (MODLE_UNLIKELY(c().log_model_internal_state)) {
        assert(s.model_state_logger);
        Simulation::dump_stats(s.id, s.epoch, s.cell_id, !s.burnin_completed, *s.interval,
                               s.get_lefs(), s.barriers, s.get_rev_collisions(),
                               s.get_fwd_collisions(), *s.model_state_logger);
      }

      // Select LEFs to be released in the current epoch and release them
      this->release_lefs(s.get_lefs(), s.barriers, s.get_rev_collisions(), s.get_fwd_collisions(),
                         s.rand_eng, s.burnin_completed);
    }
  } catch (const std::exception& err) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Caught an exception while simulating epoch #{} in cell #{} of {}:\n{}"),
        s.epoch, s.cell_id, *s.interval, err.what()));
  }
}

void Simulation::select_and_bind_lefs(State& s) noexcept(utils::ndebug_defined()) {
  const auto lef_mask = s.get_idx_buff();
  Simulation::select_lefs_to_bind(s.get_lefs(), lef_mask);
  Simulation::bind_lefs(*s.interval, s.get_lefs(), s.get_rev_ranks(), s.get_fwd_ranks(), lef_mask,
                        s.rand_eng, s.epoch);
}

void Simulation::dump_stats(const usize task_id, const usize epoch, const usize cell_id,
                            const bool burnin, const GenomicInterval& interval,
                            const absl::Span<const Lef> lefs, const ExtrusionBarriers& barriers,
                            const absl::Span<const CollisionT> rev_collisions,
                            const absl::Span<const CollisionT> fwd_collisions,
                            compressed_io::Writer& log_writer) const
    noexcept(utils::ndebug_defined()) {
  assert(c().log_model_internal_state);
  assert(log_writer);

  const auto barriers_occupied = barriers.count_active();
  const auto effective_barrier_occupancy =
      static_cast<double>(barriers_occupied) / static_cast<double>(barriers.size());

  const auto lefs_stalled_rev =
      std::count_if(rev_collisions.begin(), rev_collisions.end(),
                    [](const auto collision) { return collision.collision_occurred(); });
  const auto lefs_stalled_fwd =
      std::count_if(fwd_collisions.begin(), fwd_collisions.end(),
                    [](const auto collision) { return collision.collision_occurred(); });
  const auto lefs_stalled_both = [&]() {
    usize stalls = 0;
    for (usize i = 0; i < lefs.size(); ++i) {
      // NOLINTNEXTLINE(readability-implicit-bool-conversion)
      stalls += rev_collisions[i].collision_occurred() && fwd_collisions[i].collision_occurred();
    }
    return stalls;
  }();

  const auto count_collisions = [&](const auto collision_type) {
    auto filter_fx = [&](const auto c) { return c.collision_occurred(collision_type); };
    const auto n1 = std::count_if(rev_collisions.begin(), rev_collisions.end(), filter_fx);
    const auto n2 = std::count_if(fwd_collisions.begin(), fwd_collisions.end(), filter_fx);
    return static_cast<usize>(n1 + n2);
  };

  const auto lef_bar_collisions = count_collisions(CollisionT::LEF_BAR);
  const auto lef_lef_primary_collisions = count_collisions(CollisionT::LEF_LEF_PRIMARY);
  const auto lef_lef_secondary_collisions = count_collisions(CollisionT::LEF_LEF_SECONDARY);

  // TODO log number of units at interval boundaries

  const auto avg_loop_size =
      stats::mean(lefs.begin(), lefs.end(), [](const auto& lef) { return lef.loop_size(); });

  // clang-format off
  log_writer.write(fmt::format(
      FMT_COMPILE("{}\t{}\t{}\t"
                  "{}\t{}\t{}\t"
                  "{}\t{}\t"
                  "{}\t{}\t{}\t"
                  "{}\t{}\t{}\t"
                  "{}\t{}\n"),
      task_id, epoch, cell_id,
      interval.chrom().name(), interval.start(), interval.end(),
      burnin ? "True" : "False", effective_barrier_occupancy,
      lefs.size(), lefs_stalled_rev, lefs_stalled_fwd,
      lefs_stalled_both, lef_bar_collisions, lef_lef_primary_collisions,
      lef_lef_secondary_collisions, avg_loop_size));
  // clang-format on
}

usize Simulation::compute_tot_target_epochs(usize nlefs, usize npixels) const noexcept {
  using SC = Config::StoppingCriterion;
  if (c().stopping_criterion == SC::simulation_epochs) {
    assert(c().target_simulation_epochs != 0 &&
           c().target_simulation_epochs != std::numeric_limits<usize>::max());
    return c().num_cells * c().target_simulation_epochs;
  }

  assert(c().stopping_criterion == SC::contact_density);
  const auto tot_target_contacts =
      std::max(1.0, std::round(c().target_contact_density * static_cast<double>(npixels)));

  const auto new_contacts_per_epoch = this->compute_contacts_per_epoch(nlefs);
  return static_cast<usize>(
      std::round(tot_target_contacts / static_cast<double>(new_contacts_per_epoch)));
}

usize Simulation::compute_contacts_per_epoch(usize nlefs) const noexcept {
  const auto extrusion_speed =
      static_cast<double>(c().rev_extrusion_speed + c().fwd_extrusion_speed);
  const auto prob_of_contact_sampling =
      extrusion_speed / static_cast<double>(c().contact_sampling_interval);

  return static_cast<usize>(
      std::max(1.0, std::round(static_cast<double>(nlefs) * prob_of_contact_sampling)));
}

usize Simulation::compute_num_lefs(const usize size_bp) const noexcept {
  const auto size_mbp = static_cast<double>(size_bp) / Mbp;
  return std::max(usize(1), static_cast<usize>(std::round(c().number_of_lefs_per_mbp * size_mbp)));
}

void Simulation::print_status_update(const Task& t) const noexcept {
  assert(t.interval);
  auto tot_target_epochs = this->compute_tot_target_epochs(t.num_lefs, t.interval->npixels());
  spdlog::info(FMT_STRING("begin processing {}: simulating ~{} epochs across {} cells using {} "
                          "LEFs and {} barriers (~{} epochs per cell)..."),
               *t.interval, tot_target_epochs, c().num_cells, t.num_lefs,
               t.interval->barriers().size(), tot_target_epochs / c().num_cells);
}

usize Simulation::compute_num_worker_threads(usize num_threads, usize num_intervals,
                                             usize num_cells) noexcept {
  return std::min(num_threads, num_intervals * num_cells);
}
}  // namespace modle
