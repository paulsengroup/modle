#pragma once

#include <absl/container/flat_hash_map.h>
#include <absl/strings/str_cat.h>
#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration, Time
#include <absl/types/span.h>
#include <fmt/format.h>  // for FMT_STRING
#include <moodycamel/blockingconcurrentqueue.h>
#ifdef USE_XOSHIRO
#include <XoshiroCpp.hpp>
#endif

#include <cpp-sort/sorters/drop_merge_sorter.h>
#include <cpp-sort/sorters/insertion_sorter.h>

#include <algorithm>  // for max, min, sort, transform, all_of, for_each, remove_if
#include <atomic>     // for atomic, memory_order_relaxed
#include <boost/asio/thread_pool.hpp>
#include <boost/dynamic_bitset.hpp>
#include <cassert>
#include <chrono>  // for seconds
#include <cmath>   // for llround, floor, lround, sqrt
#include <condition_variable>
#include <cstdint>  // for uint*_t, UINT*_MAX
#include <cstdio>   // for stderr
#include <deque>
#include <filesystem>  // for create_directories
#include <functional>  // for greater, hash
#include <iosfwd>      // for size_t
#include <mutex>       // for mutex, unique_lock, scoped_lock
#include <numeric>     // for accumulate, partial_sum
#include <random>  // for mt19937,  bernoulli_distribution, seed_seq, discrete_distribution, uniform_int_distribution
#include <stdexcept>  // for runtime_error
#include <thread>
#include <type_traits>  // for declval

#include "modle/bed.hpp"        // for BED, Parser, BED::BED6, BED::Standard
#include "modle/chr_sizes.hpp"  // for ChrSize, Parser
#include "modle/common.hpp"
#include "modle/config.hpp"    // for config
#include "modle/contacts.hpp"  // for ContactMatrix
#include "modle/cooler.hpp"
#include "modle/extrusion_barriers.hpp"
#include "modle/extrusion_factors.hpp"
#include "modle/suppress_compiler_warnings.hpp"

namespace modle {

#ifdef USE_XOSHIRO
using PRNG = XoshiroCpp::Xoshiro256PlusPlus;
using seeder = XoshiroCpp::SplitMix64;
#else
using PRNG = std::mt19937_64;
using seeder = std::seed_seq;
#endif

Genome::Genome(const config& c, bool import_chroms)
    : _path_to_chrom_sizes(c.path_to_chr_sizes),
      _path_to_chrom_subranges(c.path_to_chr_subranges),
      _path_to_extr_barriers(c.path_to_extr_barriers_bed),
      _bin_size(c.bin_size),
      _diagonal_width(c.diagonal_width),
      _avg_lef_lifetime(c.average_lef_lifetime),
      _nlefs_per_mbp(c.number_of_lefs_per_mbp),
      _probability_of_barrier_block(c.probability_of_extrusion_barrier_block),
      _probability_of_lef_rebind(c.probability_of_lef_rebind),
      _probability_of_extr_unit_bypass(c.probability_of_extrusion_unit_bypass),
      _soft_stall_multiplier(c.soft_stall_multiplier),
      _hard_stall_multiplier(c.hard_stall_multiplier),
      _allow_lef_lifetime_extension(c.allow_lef_lifetime_extension),
      _sampling_interval(c.contact_sampling_interval),
      _randomize_contact_sampling(c.randomize_contact_sampling_interval),
      _nthreads(std::min(std::thread::hardware_concurrency(), c.nthreads)),
      _seed(c.seed),
      _chromosomes(import_chroms ? import_chromosomes(_path_to_chrom_sizes, _path_to_extr_barriers,
                                                      _path_to_chrom_subranges)
                                 : Chromosomes{}) {}

std::size_t Genome::size() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
      [](auto accumulator, const auto& chrom) { return accumulator + chrom.size(); });
}

std::size_t Genome::simulated_size() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
                         [](auto accumulator, const auto& chrom) {
                           return accumulator + (chrom.end_pos() - chrom.start_pos());
                         });
}

// TODO Add flag to import skip chrom without barriers
Genome::Chromosomes Genome::import_chromosomes(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_extr_barriers,
    const std::filesystem::path& path_to_chrom_subranges) {
  assert(!path_to_chrom_sizes.empty());
  assert(!path_to_extr_barriers.empty());
  Chromosomes chroms;

  {
    // Parse chrom subranges from BED. We parse everything at once to deal with duplicate entries
    absl::flat_hash_map<std::string, std::pair<uint64_t, uint64_t>> chrom_ranges;
    if (!path_to_chrom_subranges.empty()) {
      for (auto&& record : modle::bed::Parser(path_to_chrom_subranges).parse_all()) {
        chrom_ranges.emplace(std::move(record.chrom),
                             std::make_pair(record.chrom_start, record.chrom_end));
      }
    }

    // Parse chrom. sizes and build the set of chromosome to be simulated.
    // When the BED file with the chrom. subranges is not provided, all the chromosome in the
    // chrom.sizes file will be selected and returned. When a BED file with the chrom. subranges is
    // available, then only chromosomes that are present in both files will be selected. Furthermore
    // we are also checking that the subrange lies within the coordinates specified in the chrom.
    // sizes file
    for (auto&& chrom : modle::chr_sizes::Parser(path_to_chrom_sizes).parse_all()) {
      if (auto match = chrom_ranges.find(chrom.name); match != chrom_ranges.end()) {
        const auto& range_start = match->second.first;
        const auto& range_end = match->second.second;
        if (range_start < chrom.start || range_end > chrom.end) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("According to the chrom.sizes file {}, chromosome '{}' should have a size "
                         "of '{}', but the subrange specified in BED file {} extends past this "
                         "region: range {}:{}-{} does not fit in range {}:{}-{}"),
              path_to_chrom_sizes, chrom.name, chrom.end, path_to_chrom_subranges, chrom.name,
              range_start, range_end, chrom.name, chrom.start, chrom.end));
        }
        chrom.start = range_start;
        chrom.end = range_end;
        chroms.emplace(std::move(chrom));
      } else if (chrom_ranges.empty()) {
        chroms.emplace(std::move(chrom));
      }
    }
  }

  // Parse all the records from the BED file. parse_all() will throw in case of duplicates.
  // This for loop selects extrusion barriers that fall within the chromosomes to be simulated
  for (auto&& record : modle::bed::Parser(path_to_extr_barriers).parse_all()) {
    if (record.score < 0 || record.score > 1) {
      throw std::runtime_error(
          fmt::format("Invalid score field detected for record {}[{}-{}]: expected a score "
                      "between 0 and 1, got {:.4g}.",
                      record.chrom, record.chrom_start, record.chrom_end, record.score));
    }
    if (auto match = chroms.find(record.chrom); match != chroms.end()) {
      match->add_extrusion_barrier(record);
    }
  }
  return chroms;
}

std::vector<ExtrusionBarrier> Genome::allocate_barriers(const Chromosome* const chrom,
                                                        double default_prob_of_block) {
  std::vector<ExtrusionBarrier> barriers;
  std::size_t barriers_skipped = 0;
  for (const auto& b : chrom->get_barriers()) {
    if (b.strand == '+' || b.strand == '-') {
      const auto pos =
          static_cast<Bp>(std::round(static_cast<double>(b.chrom_start + b.chrom_end) / 2.0));
      const auto pblock = b.score != 0 ? b.score : default_prob_of_block;
      barriers.emplace_back(pos, pblock, b.strand);
    } else {
      ++barriers_skipped;
    }
  }

  fmt::print(stderr,
             FMT_STRING("Instantiated {} extr. barriers for '{}' ({} barriers were skipped)\n"),
             barriers.size(), chrom->name(), barriers_skipped);

  return barriers;
}

void Genome::simulate_extrusion(const std::filesystem::path& output_path, uint32_t ncells,
                                uint32_t simulation_rounds, double target_contact_density) {
  auto tpool = Genome::instantiate_thread_pool(this->_nthreads + 1, false);
  std::atomic<bool> end_of_simulation = false;

  std::mutex progress_mutex;
  std::mutex barrier_mutex;
  std::deque<std::pair<Chromosome*, std::size_t>> progress_queue;
  std::deque<std::shared_ptr<std::vector<ExtrusionBarrier>>> barriers;

  moodycamel::BlockingConcurrentQueue<Genome::Task> task_queue(ncells * 2, 1, 0);
  moodycamel::ProducerToken ptok(task_queue);

  boost::asio::post(tpool, [&]() {
    Chromosome* chrom_to_be_written = nullptr;
    auto c = output_path.empty() ? nullptr
                                 : std::make_unique<cooler::Cooler>(
                                       output_path, cooler::Cooler::WRITE_ONLY, this->_bin_size);
    while (true) {
      std::this_thread::sleep_for(std::chrono::milliseconds(25));
      {
        std::scoped_lock l(progress_mutex);
        if (auto& [chrom, count] = progress_queue.front(); chrom == nullptr) {
          end_of_simulation = true;
          return;
        } else if (count == ncells) {
          chrom_to_be_written = chrom;
          progress_queue.pop_front();
        } else {
          assert(count < ncells);
          continue;
        }
      }
      try {
        if (c) {
          fmt::print(stderr, "Writing contacts for '{}' to file {}...\n",
                     chrom_to_be_written->name(), c->get_path());
          c->write_or_append_cmatrix_to_file(
              chrom_to_be_written->contacts(), chrom_to_be_written->name(),
              chrom_to_be_written->start_pos(), chrom_to_be_written->end_pos(),
              chrom_to_be_written->size(), true);
          fmt::print(stderr, "Written {} contacts for '{}' in {:.2f}M pixels to file {}.\n",
                     chrom_to_be_written->contacts().get_tot_contacts(),
                     chrom_to_be_written->name(),
                     static_cast<double>(chrom_to_be_written->contacts().npixels()) / 1.0e6,
                     c->get_path());
        }
      } catch (const std::runtime_error& err) {
        throw std::runtime_error(fmt::format(
            FMT_STRING(
                "The following error occurred while writing contacts for '{}' to file {}: {}"),
            chrom_to_be_written->name(), c->get_path(), err.what()));
      }
      chrom_to_be_written->deallocate_contacts();
      std::scoped_lock l(barrier_mutex);
      barriers.pop_front();
    }
  });

  for (auto i = 0UL; i < this->_nthreads; ++i) {
    moodycamel::ConsumerToken ctok(task_queue);
    boost::asio::post(tpool, [&, ctok = std::move(ctok),
                              progress = progress_queue.begin()]() mutable {
      fmt::print(stderr, "Spawning simulation thread {}...\n", std::this_thread::get_id());
      try {
        std::array<Task, 256> task_buff;
        std::vector<Lef> lefs;
        // TODO use overload with timeout. Check conditional var after timeout expires
        while (true) {
          const auto avail_tasks = task_queue.wait_dequeue_bulk_timed(
              ctok, task_buff.begin(), task_buff.size(), std::chrono::milliseconds(10));
          if (avail_tasks == 0) {
            if (end_of_simulation) {
              return;
            }
            continue;
          }
          auto tasks = absl::MakeSpan(task_buff.data(), avail_tasks);
          for (auto& task : tasks) {
            if (task.id == std::numeric_limits<decltype(task.id)>::max()) {
              return;
            }
            lefs.resize(task.nlefs);
            std::for_each(lefs.begin(), lefs.end(), [](auto& lef) { lef.reset(); });

            if (target_contact_density != 0) {
              task.nrounds = static_cast<std::size_t>(std::round(
                  (target_contact_density * static_cast<double>(this->_sampling_interval *
                                                                task.chrom->contacts().npixels())) /
                  static_cast<double>(ncells * task.nlefs)));
            }

            if (task.cell_id == 0) {
              fmt::print(stderr, FMT_STRING("Simulating {} epochs for '{}' across {} cells...\n"),
                         task.nrounds, task.chrom->name(), ncells);
            }

            Genome::simulate_extrusion_kernel(task.chrom, task.cell_id, task.nrounds, lefs,
                                              task.barriers);
            std::scoped_lock l(progress_mutex);
            if (progress == progress_queue.end() || progress->first != task.chrom) {
              progress = std::find_if(progress_queue.begin(), progress_queue.end(),
                                      [&](const auto& p) { return task.chrom == p.first; });
              assert(progress != progress_queue.end());
            }
            if (++progress->second == ncells - 1) {
              fmt::print(stderr, "Simulation for '{}' successfully completed.\n",
                         task.chrom->name());
            }
          }
        }
      } catch (const std::exception& err) {
        fmt::print(stderr, FMT_STRING("Detected an error in thread {}:\n{}\n"),
                   std::this_thread::get_id(), err.what());
#ifndef BOOST_STACKTRACE_USE_NOOP
        const auto* st = boost::get_error_info<modle::utils::traced>(err);
        if (st) {
          std::cerr << *st << '\n';
        } else {
          fmt::print(stderr, "Stack trace not available!\n");
        }
#endif
        throw;
      }
    });
  }

  std::vector<Chromosome*> chroms(static_cast<std::size_t>(this->_chromosomes.size()));
  std::transform(this->_chromosomes.begin(), this->_chromosomes.end(), chroms.begin(),
                 [](auto& chrom) { return &chrom; });

  std::sort(chroms.begin(), chroms.end(),
            [](const auto* const c1, const auto* const c2) { return c1->size() > c2->size(); });

  std::shared_ptr<std::vector<ExtrusionBarrier>> extr_barriers_buff{nullptr};
  std::vector<Task> tasks(ncells);
  std::size_t taskid = 0;
  for (auto* chrom : chroms) {
    if (!chrom->ok()) {
      continue;
    }
    chrom->allocate_contacts(this->_bin_size, this->_diagonal_width);
    {
      std::scoped_lock l(barrier_mutex, progress_mutex);
      progress_queue.emplace_back(chrom, 0UL);
      barriers.emplace_back(std::make_shared<std::vector<ExtrusionBarrier>>(
          Genome::allocate_barriers(chrom, this->_probability_of_barrier_block)));
      extr_barriers_buff = barriers.back();
    }

    tasks.resize(ncells);
    const auto nlefs = static_cast<std::size_t>(
        std::round(this->_nlefs_per_mbp * (static_cast<double>(chrom->simulated_size()) / 1.0e6)));
    std::generate(tasks.begin(), tasks.end(), [&, i = 0UL]() mutable {
      return Task{taskid++, chrom, i++, simulation_rounds, nlefs, extr_barriers_buff};
    });

    while (
        !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), tasks.size())) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  }
  progress_queue.emplace_back(nullptr, 0UL);  // Signal end of simulation
  tpool.join();
}

void Genome::simulate_extrusion(const std::filesystem::path& output_path, uint32_t ncells,
                                uint32_t simulation_rounds) {
  this->simulate_extrusion(output_path, ncells, simulation_rounds, 0.0);
}

void Genome::simulate_extrusion(const std::filesystem::path& output_path, uint32_t ncells,
                                double target_contact_density) {
  this->simulate_extrusion(output_path, ncells, 0, target_contact_density);
}

void Genome::simulate_extrusion_kernel(
    Chromosome* chrom, std::size_t cell_id, std::size_t simulation_rounds,
    std::vector<Lef> lef_buff,
    std::shared_ptr<const std::vector<ExtrusionBarrier>> extr_barrier_buff) {
  // Bind all lefs
  const auto seed = this->_seed + std::hash<std::string_view>{}(chrom->name()) +
                    std::hash<std::size_t>{}(chrom->size()) + std::hash<std::size_t>{}(cell_id);
#ifdef USE_XOSHIRO
  modle::seeder seeder_(seed);
  modle::PRNG rand_eng(seeder_.generateSeedSequence<4>());
#else
  modle::seeder seeder_{seed};
  modle::PRNG rand_eng(seeder_);
#endif

  std::vector<std::size_t> lef_initial_loading_round(lef_buff.size());
  {
    absl::btree_multiset<std::size_t> sorted_rounds;
    std::uniform_int_distribution<std::size_t> round_gen{
        0, (4 * this->_avg_lef_lifetime) / this->_bin_size};

    std::generate_n(std::inserter(sorted_rounds, sorted_rounds.begin()), lef_buff.size(),
                    [&]() { return round_gen(rand_eng); });
    std::transform(sorted_rounds.begin(), sorted_rounds.end(), lef_initial_loading_round.rbegin(),
                   [&](const auto n) { return n - *sorted_rounds.begin(); });
  }
  const auto burnin_rounds = lef_initial_loading_round.front();
  simulation_rounds += burnin_rounds;

  std::vector<std::size_t> rev_lef_rank_buff(lef_buff.size());
  std::iota(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), 0);
  std::vector<std::size_t> fwd_lef_rank_buff{rev_lef_rank_buff};
  boost::dynamic_bitset<> mask{lef_buff.size()};

  // this->bind_all_lefs(chrom, lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rand_eng);

  std::vector<collision_t> rev_lef_collision_buff(lef_buff.size());
  std::vector<collision_t> fwd_lef_collision_buff(lef_buff.size());
  std::bernoulli_distribution sample_contacts{1.0 / this->_sampling_interval};

  auto lefs = absl::MakeSpan(lef_buff);
  const auto barriers = absl::MakeConstSpan(*extr_barrier_buff);
  auto rev_lef_ranks = absl::MakeSpan(rev_lef_rank_buff);
  auto fwd_lef_ranks = absl::MakeSpan(fwd_lef_rank_buff);
  auto rev_lef_collisions = absl::MakeSpan(rev_lef_collision_buff);
  auto fwd_lef_collisions = absl::MakeSpan(fwd_lef_collision_buff);

  for (auto round = 0UL; round < simulation_rounds; ++round) {
    const auto register_contacts =
        round > burnin_rounds &&
        (this->_randomize_contact_sampling ? sample_contacts(rand_eng)
                                           : round % this->_sampling_interval == 0);

    if (!lef_initial_loading_round.empty()) {
      while (!lef_initial_loading_round.empty() && lef_initial_loading_round.back() == round) {
        lef_initial_loading_round.pop_back();
      }
      const auto nlefs = lef_buff.size() - lef_initial_loading_round.size();

      lefs = absl::MakeSpan(lef_buff.data(), nlefs);
      rev_lef_ranks = absl::MakeSpan(rev_lef_rank_buff.data(), nlefs);
      fwd_lef_ranks = absl::MakeSpan(fwd_lef_rank_buff.data(), nlefs);
      rev_lef_collisions = absl::MakeSpan(rev_lef_collision_buff.data(), nlefs);
      fwd_lef_collisions = absl::MakeSpan(fwd_lef_collision_buff.data(), nlefs);
    }

    this->select_lefs_to_bind(lefs, mask);
    this->bind_lefs(chrom, lefs, rev_lef_ranks, fwd_lef_ranks, rand_eng, mask);

    this->check_lef_lef_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, rev_lef_collisions,
                                   fwd_lef_collisions);
    this->apply_lef_lef_stalls(lefs, rev_lef_collisions, fwd_lef_collisions, rev_lef_ranks,
                               fwd_lef_ranks, rand_eng, this->_probability_of_extr_unit_bypass);

    this->check_lef_bar_collisions(lefs, rev_lef_ranks, fwd_lef_ranks, barriers, rev_lef_collisions,
                                   fwd_lef_collisions);
    this->apply_lef_bar_stalls(lefs, rev_lef_collisions, fwd_lef_collisions, barriers, rand_eng,
                               this->_soft_stall_multiplier, this->_hard_stall_multiplier);

    if (register_contacts) {
      this->register_contacts(chrom, lefs);
    }

    this->extrude(chrom, lefs);
  }
}

template <typename MaskT>
void Genome::bind_lefs(const Chromosome* const chrom, absl::Span<Lef> lefs,
                       absl::Span<std::size_t> rev_lef_rank_buff,
                       absl::Span<std::size_t> fwd_lef_rank_buff, modle::PRNG& rand_eng,
                       MaskT& mask) {
  static_assert(std::is_same_v<boost::dynamic_bitset<>, MaskT> ||
                    std::is_integral_v<std::decay_t<decltype(std::declval<MaskT&>().operator[](
                        std::declval<std::size_t>()))>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");

  assert(lefs.size() <= mask.size() || mask.empty());
  assert(chrom);
  chrom_pos_generator_t pos_generator{chrom->start_pos(), chrom->end_pos()};
  lef_lifetime_generator_t lifetime_generator{1.0 / static_cast<double>(this->_avg_lef_lifetime)};
  for (auto i = 0UL; i < lefs.size(); ++i) {
    // Here we only consider even indexes, as binding a single extr. unit is not
    // allowed: either we bind all extr. units belonging to a LEF, or we bind none
    if (mask.empty() || mask[i]) {
      lefs[i].bind_at_pos(pos_generator(rand_eng), lifetime_generator(rand_eng));
    }
  }
  this->rank_lefs(lefs, rev_lef_rank_buff, fwd_lef_rank_buff);
}

void Genome::bind_all_lefs(const Chromosome* chrom, absl::Span<Lef> lefs,
                           absl::Span<std::size_t> rev_lef_rank_buff,
                           absl::Span<std::size_t> fwd_lef_rank_buff, PRNG& rand_eng) {
  chrom_pos_generator_t chrom_pos_generator{chrom->start_pos(), chrom->end_pos()};
  lef_lifetime_generator_t lef_lifetime_generator{1.0 /
                                                  static_cast<double>(this->_avg_lef_lifetime)};
  absl::btree_multiset<Bp> positions;
  std::generate_n(std::inserter(positions, positions.begin()), lefs.size(),
                  [&]() { return chrom_pos_generator(rand_eng); });
  std::size_t i = 0;
  for (auto& pos : positions) {
    lefs[i].bind_at_pos(pos, lef_lifetime_generator(rand_eng));
    rev_lef_rank_buff[i] = i;
    fwd_lef_rank_buff[i] = i;
    ++i;
  }
}

void Genome::rank_lefs(absl::Span<Lef> lefs, absl::Span<std::size_t> rev_lef_rank_buff,
                       absl::Span<std::size_t> fwd_lef_rank_buff, bool init_buffers) {
  if (init_buffers) {
    std::iota(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), 0);
    std::copy(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), rev_lef_rank_buff.begin());
  }
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());

#if 0  // 55.77s
  std::sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
  });

  std::sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(), [&](const auto r1, const auto r2) {
    return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
  });
#else  // 35.50s
  if (lefs.size() < 20) {  // TODO Come up with a reasonable threshold
    cppsort::insertion_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                            [&](const auto r1, const auto r2) {
                              return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                            });

    cppsort::insertion_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                            [&](const auto r1, const auto r2) {
                              return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                            });
  } else {
    cppsort::drop_merge_sort(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                             [&](const auto r1, const auto r2) {
                               return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                             });

    cppsort::drop_merge_sort(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                             [&](const auto r1, const auto r2) {
                               return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                             });
  }
#endif
}

void Genome::extrude(const Chromosome* chrom, absl::Span<Lef> lefs) {
  for (auto& lef : lefs) {
    if (lef.lifetime == 0) {
      continue;
    }
    if (lef.rev_unit.stalled()) {
      lef.rev_unit.decrement_stalls(this->_bin_size);
    } else if (lef.rev_unit._pos < chrom->start_pos() + this->_bin_size) {
      lef.rev_unit._pos = chrom->start_pos();
      lef.rev_unit._nstalls_lef_lef = std::numeric_limits<Bp>::max();
    } else {
      lef.rev_unit._pos -= this->_bin_size;
    }

    if (lef.fwd_unit.stalled()) {
      lef.fwd_unit.decrement_stalls(this->_bin_size);
    } else if (lef.fwd_unit._pos + this->_bin_size > chrom->end_pos() - 1) {
      lef.fwd_unit._pos = chrom->end_pos() - 1;
      lef.fwd_unit._nstalls_lef_lef = std::numeric_limits<Bp>::max();
    } else {
      lef.fwd_unit._pos += this->_bin_size;
    }

    if ((lef.lifetime -= std::min(lef.lifetime, 2 * this->_bin_size)) == 0) {
      lef.release();
    }
  }
}

void Genome::check_lef_lef_collisions(absl::Span<const Lef> lefs,
                                      absl::Span<const std::size_t> rev_lef_rank_buff,
                                      absl::Span<const std::size_t> fwd_lef_rank_buff,
                                      absl::Span<collision_t> rev_collision_buff,
                                      absl::Span<collision_t> fwd_collision_buff,
                                      Bp dist_threshold) {
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());
  assert(lefs.size() == fwd_collision_buff.size());
  assert(lefs.size() == rev_collision_buff.size());
  assert(std::is_sorted(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                        }));
  assert(std::is_sorted(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                        }));

  // Clear buffers
  std::fill(fwd_collision_buff.begin(), fwd_collision_buff.end(), 0);
  std::fill(rev_collision_buff.begin(), rev_collision_buff.end(), 0);

  /* Loop over lefs, using a procedure similar to merge in mergesort
   * The idea here is that if we have a way to visit extr. units in 5'-3' order, then detecting
   * LEF-LEF collisions boils down to:
   *  - Starting from the first fwd extr. unit in 5'-3' order:
   *    - Look for rev extr. units that are located downstream of the fwd unit that is being
   *      processed
   *    - Continue looking for the next rev unit, until the distance between the fwd and rev units
   *      being considered is less than a certain threshold (e.g. the bin size)
   *    - While doing so, increase the number of LEF-LEF collisions for the fwd/rev unit that are
   *      being processed
   */
  std::size_t i = 0, j = 0;
  while (i < lefs.size() && j < lefs.size()) {
    const auto& fwd_idx = fwd_lef_rank_buff[i];      // index of the ith fwd unit in 5'-3' order
    const auto& fwd = lefs[fwd_idx].fwd_unit.pos();  // pos of the ith unit

    for (auto k = j; k < lefs.size(); ++k) {
      const auto& rev_idx = rev_lef_rank_buff[k];      // index of the kth rev unit in 5'-3' order
      const auto& rev = lefs[rev_idx].rev_unit.pos();  // pos of the kth unit
      if (rev < fwd) {  // Look for the first rev unit that comes after the current fwd unit
        ++j;  // NOTE k is increased by the for loop. Increasing j does not affect the current for
              // loop
        continue;
      }

      // Note: this delta between uint is safe to do, because at this point rev_pos >= fwd_pos
      const auto delta = rev - fwd;

      // Increment the # of collisions if delta is less than the threshold, and if the two ext.
      // units do not belong to the same LEF
      fwd_collision_buff[fwd_idx] += delta <= dist_threshold && rev_idx != fwd_idx;
      rev_collision_buff[rev_idx] += delta <= dist_threshold && rev_idx != fwd_idx;

      /*
      fmt::print(stderr,
                 "i={}; j={}; fwd_idx={}; rev_idx={}; rev={}; fwd={}; delta={}; fwd_mask={}; "
                 "rev_mask={};\n",
                 i, j, fwd_idx, rev_idx, rev, fwd, delta, fwd_collision_buff[fwd_idx],
                 rev_collision_buff[rev_idx]);
      */
      // Break out of the loop if the units being processed are too far from each other, or if k
      // points to the last rev unit
      if (delta > dist_threshold || k == lefs.size() - 1) {
        // We always move to the next fwd unit, because when this branch is executed, we have
        // already detected all possible collisions for the current fwd unit
        ++i;

        // If we already know that rev == fwd, we can spare one loop iteration by increasing j
        // directly here
        j += rev == fwd;
        break;
      }
    }
  }
}

void Genome::check_lef_lef_collisions(absl::Span<const Lef> lefs,
                                      absl::Span<const std::size_t> rev_lef_rank_buff,
                                      absl::Span<const std::size_t> fwd_lef_rank_buff,
                                      absl::Span<collision_t> rev_collision_buff,
                                      absl::Span<collision_t> fwd_collision_buff) {
  this->check_lef_lef_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, rev_collision_buff,
                                 fwd_collision_buff, this->_bin_size);
}

void Genome::apply_lef_lef_stalls(absl::Span<Lef> lefs,
                                  absl::Span<const collision_t> rev_collision_buff,
                                  absl::Span<const collision_t> fwd_collision_buff,
                                  absl::Span<const std::size_t> rev_lef_rank_buff,
                                  absl::Span<const std::size_t> fwd_lef_rank_buff, PRNG& rand_eng,
                                  double prob_of_bypass) {
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());
  assert(lefs.size() == fwd_collision_buff.size());
  assert(lefs.size() == rev_collision_buff.size());
  assert(std::is_sorted(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                        }));
  assert(std::is_sorted(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                        }));
  assert(std::accumulate(fwd_collision_buff.begin(), fwd_collision_buff.end(), 0UL) ==
         std::accumulate(rev_collision_buff.begin(), rev_collision_buff.end(), 0UL));

  lef_lef_stall_generator_t rng(prob_of_bypass);
  std::size_t i = 0, j = 0;
  auto fwd_idx = fwd_lef_rank_buff[i];
  auto rev_idx = rev_lef_rank_buff[j];

  std::size_t collisions_fwd = fwd_collision_buff[fwd_idx];
  std::size_t collisions_rev = rev_collision_buff[rev_idx];

  /* Loop over the two buffers storing the number of collision for each LEF. We use the fwd/rev
   * ranks to iterate through collision events in 5'-3' direction. Draw and apply the appropriate
   * number of stalls
   */
  while (true) {
    while (collisions_fwd == 0) {  // Find the first collision event for fwd units
      if (++i >= lefs.size()) {
        assert(collisions_rev == 0);
        return;  // All collisions have been processed
      }
      fwd_idx = fwd_lef_rank_buff[i];
      collisions_fwd = fwd_collision_buff[fwd_idx];
      // Reset the number of stalls if the number of collisions for the current extr. unit is 0
      lefs[fwd_idx].fwd_unit._nstalls_lef_lef *= collisions_fwd != 0;
    }

    while (collisions_rev == 0) {  // Find the first collision event for rev units
      if (++j >= lefs.size()) {
        assert(collisions_fwd == 0);
        return;  // All collisions have been processed
      }
      rev_idx = rev_lef_rank_buff[j];
      collisions_rev = rev_collision_buff[rev_idx];
      // Reset the number of stalls if the number of collisions for the current extr. unit is 0
      lefs[rev_idx].rev_unit._nstalls_lef_lef *= collisions_rev != 0;
    }

    // Each iteration consumes a pair of collision events. When fwd or rev collision events (or
    // both) reach 0, this while loop is terminated, and we go back to the outer while loop (where
    // we look for the next non-zero collision event)
    while (collisions_fwd > 0 && collisions_rev > 0) {
      const auto nstalls = rng(rand_eng);
      lefs[fwd_idx].fwd_unit._nstalls_lef_lef += nstalls;
      lefs[rev_idx].rev_unit._nstalls_lef_lef += nstalls;
      --collisions_fwd;
      --collisions_rev;
    }
  }
}

void Genome::check_lef_bar_collisions(absl::Span<const Lef> lefs,
                                      absl::Span<const std::size_t> rev_lef_rank_buff,
                                      absl::Span<const std::size_t> fwd_lef_rank_buff,
                                      absl::Span<const ExtrusionBarrier> extr_barriers,
                                      absl::Span<collision_t> rev_collision_buff,
                                      absl::Span<collision_t> fwd_collision_buff,
                                      Bp dist_threshold) {
  assert(lefs.size() == fwd_lef_rank_buff.size());
  assert(lefs.size() == rev_lef_rank_buff.size());
  assert(lefs.size() == fwd_collision_buff.size());
  assert(lefs.size() == rev_collision_buff.size());
  assert(std::is_sorted(fwd_lef_rank_buff.begin(), fwd_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].fwd_unit.pos() < lefs[r2].fwd_unit.pos();
                        }));
  assert(std::is_sorted(rev_lef_rank_buff.begin(), rev_lef_rank_buff.end(),
                        [&](const auto r1, const auto r2) {
                          return lefs[r1].rev_unit.pos() < lefs[r2].rev_unit.pos();
                        }));

  // Clear buffers (-1 indicates no collisions)
  std::fill(fwd_collision_buff.begin(), fwd_collision_buff.end(),
            std::numeric_limits<std::size_t>::max());
  std::fill(rev_collision_buff.begin(), rev_collision_buff.end(),
            std::numeric_limits<std::size_t>::max());

  /* Loop over lefs, using a procedure similar to merge in mergesort
   * The idea here is that if we have a way to visit extr. units in 5'-3' order, then detecting
   * LEF-LEF collisions boils down to:
   *  - Starting from the first fwd extr. unit in 5'-3' order:
   *    - Look for rev extr. units that are located downstream of the fwd unit that is being
   *      processed
   *    - Continue looking for the next rev unit, until the distance between the fwd and rev units
   *      being considered is less than a certain threshold (e.g. the bin size)
   *    - While doing so, increase the number of LEF-LEF collisions for the fwd/rev unit that are
   *      being processed
   * The fwd and rev_collision_buff will contain the index corresponding to the ext. barrier that
   * caused the collision, or -1 in case of no collisions
   */
  std::size_t j = 0;
  for (auto i = 0UL; i < lefs.size(); ++i) {
    const auto& fwd_idx = fwd_lef_rank_buff[i];  // index of the ith fwd unit in 5'-3' order
    if (lefs[fwd_idx].fwd_unit.stalled()) {
      continue;
    }
    const auto& fwd_pos = lefs[fwd_idx].fwd_unit.pos();  // pos of the ith unit
    auto extr_barr_pos = extr_barriers[j].pos();

    while (extr_barr_pos < fwd_pos) {
      if (++j == extr_barriers.size()) {
        goto endloop1;
      }
      extr_barr_pos = extr_barriers[j].pos();
    }

    const auto delta = extr_barr_pos - fwd_pos;
    fwd_collision_buff[fwd_idx] = delta <= dist_threshold ? j : -1UL;
  }
endloop1:

  j = extr_barriers.size() - 1;
  for (auto i = lefs.size() - 1; i > 0; --i) {
    const auto& rev_idx = rev_lef_rank_buff[i];  // index of the ith fwd unit in 5'-3' order
    if (lefs[rev_idx].rev_unit.stalled()) {
      continue;
    }
    const auto& rev_pos = lefs[rev_idx].rev_unit.pos();  // pos of the ith unit
    auto extr_barr_pos = extr_barriers[j].pos();

    while (extr_barr_pos > rev_pos) {
      if (j-- == 0) {
        return;
      }
      extr_barr_pos = extr_barriers[j].pos();
    }

    const auto delta = rev_pos - extr_barr_pos;
    rev_collision_buff[rev_idx] = delta <= dist_threshold ? j : -1UL;
  }
}

void Genome::check_lef_bar_collisions(absl::Span<const Lef> lefs,
                                      absl::Span<const std::size_t> rev_lef_rank_buff,
                                      absl::Span<const std::size_t> fwd_lef_rank_buff,
                                      absl::Span<const ExtrusionBarrier> extr_barriers,
                                      absl::Span<collision_t> rev_collision_buff,
                                      absl::Span<collision_t> fwd_collision_buff) {
  this->check_lef_bar_collisions(lefs, rev_lef_rank_buff, fwd_lef_rank_buff, extr_barriers,
                                 rev_collision_buff, fwd_collision_buff, this->_bin_size);
}

void Genome::apply_lef_bar_stalls(absl::Span<Lef> lefs,
                                  absl::Span<const collision_t> rev_collision_buff,
                                  absl::Span<const collision_t> fwd_collision_buff,
                                  absl::Span<const ExtrusionBarrier> extr_barriers, PRNG& rand_eng,
                                  double soft_stall_multiplier, double hard_stall_multiplier) {
  assert(lefs.size() == fwd_collision_buff.size());
  assert(lefs.size() == rev_collision_buff.size());
  constexpr auto NO_COLLISION =
      std::numeric_limits<std::remove_pointer_t<decltype(rev_collision_buff.data())>>::max();

  for (auto i = 0UL; i < lefs.size(); ++i) {
    if ((rev_collision_buff[i] == NO_COLLISION && fwd_collision_buff[i] == NO_COLLISION) ||
        (lefs[i].rev_unit.stalled() && lefs[i].fwd_unit.stalled())) {
      continue;
    }
    auto& lef = lefs[i];

    if (const auto barr_idx = rev_collision_buff[i]; barr_idx != NO_COLLISION) {
      assert(!lef.rev_unit.stalled());
      const auto& barrier = extr_barriers[barr_idx];
      auto nstalls = lef_bar_stall_generator_t{1.0 - barrier.pblock()}(rand_eng);
      if (barrier.blocking_direction() != dna::rev) {
        assert(barrier.blocking_direction() == dna::fwd);
        nstalls = static_cast<Bp>(std::round(soft_stall_multiplier * static_cast<double>(nstalls)));
      }
      lef.rev_unit._nstalls_lef_bar = nstalls;
    }

    if (const auto barr_idx = fwd_collision_buff[i]; barr_idx != NO_COLLISION) {
      assert(!lef.fwd_unit.stalled());
      const auto& barrier = extr_barriers[barr_idx];
      auto nstalls = lef_bar_stall_generator_t{1.0 - barrier.pblock()}(rand_eng);
      if (barrier.blocking_direction() != dna::fwd) {
        assert(barrier.blocking_direction() == dna::rev);
        nstalls = static_cast<Bp>(std::round(soft_stall_multiplier * static_cast<double>(nstalls)));
      }
      lef.fwd_unit._nstalls_lef_bar = nstalls;
    }

    const auto rev_barr_idx = rev_collision_buff[i];
    const auto fwd_barr_idx = fwd_collision_buff[i];

    if (lef.rev_unit.stalled() && lef.fwd_unit.stalled() && rev_barr_idx != NO_COLLISION &&
        fwd_barr_idx != NO_COLLISION &&
        extr_barriers[rev_barr_idx].blocking_direction() == dna::rev &&
        extr_barriers[fwd_barr_idx].blocking_direction() == dna::fwd) {
      const auto nstalls = static_cast<Bp>(std::round(
          hard_stall_multiplier * static_cast<double>(std::min(lef.rev_unit._nstalls_lef_bar,
                                                               lef.fwd_unit._nstalls_lef_bar))));
      lef.rev_unit._nstalls_lef_bar += nstalls;
      lef.fwd_unit._nstalls_lef_bar += nstalls;

      if (this->_allow_lef_lifetime_extension) {
        assert(std::numeric_limits<Bp>::max() - nstalls > lef.lifetime);  // NOLINT
        lef.lifetime += nstalls;
      }
    }
  }
}

void Genome::register_contacts(Chromosome* chrom, absl::Span<const Lef> lefs) {
  for (const auto& lef : lefs) {
    if (lef.is_bound() && lef.rev_unit.pos() != chrom->start_pos() &&
        lef.rev_unit.pos() != chrom->end_pos() - 1 && lef.fwd_unit.pos() != chrom->start_pos() &&
        lef.fwd_unit.pos() != chrom->end_pos() - 1) {
      chrom->increment_contacts(lef.rev_unit.pos(), lef.fwd_unit.pos(), this->_bin_size);
    }
  }
}

template <typename MaskT>
void Genome::select_lefs_to_bind(absl::Span<const Lef> lefs, MaskT& mask) {
  static_assert(std::is_same_v<boost::dynamic_bitset<>, MaskT> ||
                    std::is_integral_v<std::decay_t<decltype(std::declval<MaskT&>().operator[](
                        std::declval<std::size_t>()))>>,
                "mask should be a vector of integral numbers or a boost::dynamic_bitset.");
  assert(lefs.size() <= mask.size());
  for (auto i = 0UL; i < lefs.size(); ++i) {
    mask[i] = lefs[i].lifetime == 0;
  }
}

boost::asio::thread_pool Genome::instantiate_thread_pool() const {
  return boost::asio::thread_pool(this->_nthreads);
}

template <typename I>
boost::asio::thread_pool Genome::instantiate_thread_pool(I nthreads, bool clamp_nthreads) {
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
