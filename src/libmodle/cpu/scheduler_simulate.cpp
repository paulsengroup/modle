// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "modle/simulation.hpp"
// clang-format on
#include <fmt/format.h>     // for make_format_args, vformat_to, FMT_STRING
#include <spdlog/spdlog.h>  // for info
#include <xxhash.h>

#include <algorithm>   // for max, copy, min, find_if, generate
#include <cassert>     // for assert
#include <chrono>      // for microseconds, milliseconds
#include <cmath>       // for round
#include <deque>       // for deque, operator-, operator!=, _Deque_ite...
#include <exception>   // for exception_ptr, exception, current_exception
#include <filesystem>  // for path
#include <limits>      // for numeric_limits
#include <stdexcept>   // for runtime_error
#include <thread>      // IWYU pragma: keep for sleep_for
#include <utility>     // for pair
#include <vector>      // for vector

#include "modle/common/common.hpp"  // for u64
#include "modle/genome.hpp"         // for GenomicInterval, Genome

namespace modle {

void Simulation::spawn_worker_threads(const usize num_workers, const usize batch_size) {
  assert(num_workers != 0);
  assert(batch_size != 0);
  for (usize i = 0; i < num_workers; ++i) {
    this->_ctx.spawn_worker_thread(
        [this, batch_size, tid = i]() { this->simulate_worker(tid, batch_size); });
  }
}

void Simulation::spawn_io_threads() {
  this->_ctx.spawn_io_thread([this]() { this->simulate_io(); });
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Simulation::run_simulate() {
  if (!c().skip_output) {
    assert(std::filesystem::exists(c().path_to_output_prefix.parent_path()));
    if (c().force) {
      std::filesystem::remove(c().path_to_output_file_cool);
      if (c().track_1d_lef_position) {
        std::filesystem::remove(c().path_to_lef_1d_occupancy_bw_file);
      }
    }
  }

  // These are the threads spawned by run_simulate:
  // - 1 thread to write contacts to disk. This thread pops GenomicInterval* from a std::deque once
  // the
  //   simulation on the GenomicInterval* has been ultimated. This thread is also responsible of
  //   freeing memory as soon as it is not needed by any other thread.
  // - The main thread (i.e. this thread), which loops over chromosomes, instantiates data
  //   structures that are shared across simulation threads (e.g. the vector of extr. barriers),
  //   then it submits simulation tasks to a concurrent queue
  // - N simulation threads. These threads pop tasks from the concurrent queue fed by the main
  //   thread and carry out the actual simulation
  //
  // Switching to this thread-layout allowed us to significantly reduce memory allocations, as by
  // running few simulation threads for the whole simulation allows us to recycle buffers more
  // efficiently. Most of the time, the only data moved are tasks, which are implemented as
  // light-weight structs

  const auto task_batch_size = [this]() -> usize {
    if (c().num_cells <= c().nthreads) {
      return 1;
    }
    return std::min(usize(16), c().num_cells / c().nthreads);
  }();

  // Queue used to submit simulation tasks to the thread pool
  auto ptok_pending = this->_ctx.register_producer<Task::Status::PENDING>();
  auto ptok_finished = this->_ctx.register_producer<Task::Status::COMPLETED>();

  // std::mutex model_state_logger_mtx;  // Protect rw access to the log file located at
  // this->path_to_model_state_log_file
  if (c().log_model_internal_state && !c().skip_output) {
    compressed_io::Writer(c().path_to_model_state_log_file).write(model_internal_state_log_header);
  }

  try {
    this->spawn_io_threads();
    this->spawn_worker_threads(
        std::min(c().nthreads, this->_genome.num_intervals() * c().num_cells), task_batch_size);

    // The remaining code submits simulation tasks to the queue. Then it waits until all the tasks
    // have been completed and contacts have been written to disk

    usize taskid = 0;
    const auto max_sleep_time =
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::milliseconds(250));

    std::unique_ptr<XXH3_state_t, utils::XXH3_Deleter> xxh_state{XXH3_createState()};

    // Loop over chromosomes
    for (auto& interval : this->_genome) {
      this->_ctx.check_exceptions();
      auto sleep_time = std::chrono::microseconds(50);

      auto rand_eng = random::PRNG(interval.hash(*xxh_state, c().seed));

      // Don't bother simulating intervals without barriers
      if (!c().simulate_chromosomes_wo_barriers && interval.num_barriers() == 0) {
        spdlog::info(FMT_STRING("{} has 0 barriers... SKIPPING!"), interval);
        Task t{};
        t.interval = &interval;
        for (usize cellid = 0; cellid < c().num_cells; ++cellid) {
          while (!this->_ctx.try_enqueue<Task::Status::COMPLETED>(Task{t}, ptok_finished)) {
            this->_ctx.check_exceptions();
            sleep_time = std::min(max_sleep_time, sleep_time * 2);
            std::this_thread::sleep_for(sleep_time);
          }
          rand_eng.jump();
        }
        continue;
      }

      // Compute # of LEFs to be simulated based on interval.sizes
      const auto nlefs = this->compute_num_lefs(interval.size());

      const auto npixels = interval.npixels();
      const auto tot_target_contacts =
          static_cast<usize>(std::round(static_cast<double>(npixels) * c().target_contact_density));
      const auto target_epochs = this->compute_tot_target_epochs(nlefs, npixels);

      const auto target_contacts_per_cell =
          (tot_target_contacts + c().num_cells - 1) / c().num_cells;

      usize tot_target_contacts_rolling_count = 0;
      for (usize cellid = 0; cellid < c().num_cells; ++cellid) {
        // This is needed to not overshoot the target contact density
        const auto effective_target_contacts = std::min(
            target_contacts_per_cell, tot_target_contacts - tot_target_contacts_rolling_count);
        tot_target_contacts_rolling_count += effective_target_contacts;

        Task t{taskid++,
               &interval,
               cellid,
               target_epochs,
               effective_target_contacts,
               nlefs,
               rand_eng,
               Task::Status::PENDING};
        spdlog::debug(FMT_STRING("[main]: submitting task #{} ({} cell #{})..."), t.id, *t.interval,
                      t.cell_id);
        while (!this->_ctx.try_enqueue<Task::Status::PENDING>(std::move(t), ptok_pending)) {
          this->_ctx.check_exceptions();
          sleep_time = std::min(max_sleep_time, sleep_time * 2);
          std::this_thread::sleep_for(sleep_time);
        }
        rand_eng.jump();
      }
    }

    this->_ctx.shutdown();
    assert(this->_ctx.num_submitted() == this->_ctx.num_completed());
    assert(!this->_ctx.exception_thrown());
  } catch (...) {
    this->_ctx.set_exception_main(std::current_exception());
    this->_ctx.check_exceptions();
    throw;  // should not be necessary, but better safe than sorry
  }
}

static std::string format_rand_eng(const random::PRNG_t& rand_eng) {
  const auto state = rand_eng.serialize();
  static_assert(state.size() == 4);
  static_assert(sizeof(state[0]) == sizeof(u64));
  return fmt::format(FMT_STRING("0x{:016x}{:016x}{:016x}{:016x}"), state[0], state[1], state[2],
                     state[3]);
}

void Simulation::simulate_worker(const u64 tid, const usize task_batch_size) {
  spdlog::info(FMT_STRING("Spawning simulation thread {}..."), tid);

  auto ctok = this->_ctx.register_consumer<Task::Status::PENDING>();
  auto ptok = this->_ctx.register_producer<Task::Status::COMPLETED>();
  std::vector<Task> task_buff(task_batch_size);  // Tasks are dequeued in batch.
  // This is to reduce contention when accessing the queue

  // This state object owns all the buffers and PRNG + seed required in order to simulate loop
  // extrusion for a single cell. Buffers are allocated once and resized, cleared and reused
  // throughout the simulation
  Simulation::State local_state;

  const auto tmp_model_internal_state_log_path = [&]() {
    auto p = c().path_to_model_state_log_file;
    if (c().log_model_internal_state && !c().skip_output) {
      p.replace_extension(fmt::format(FMT_STRING("{}{}"), tid, p.extension().string()));
      local_state.model_state_logger = std::make_unique<compressed_io::Writer>(p);
    }
    return p;
  }();

  try {
    while (!!this->_ctx) {
      this->_ctx.wait_dequeue<Task::Status::PENDING>(ctok, task_buff);
      if (task_buff.empty() && !this->_ctx.shutdown_signal_sent()) {
        spdlog::debug(FMT_STRING("[W{}]: no tasks available: sleeping for a bit..."), tid);
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        continue;
      }
      if (task_buff.empty() && this->_ctx.shutdown_signal_sent()) {
        // Reached end of simulation (i.e. all tasks have been processed)
        // if (this->log_model_internal_state && !this->skip_output) {
        //   assert(!tmp_model_internal_state_log_path.empty());
        //   // Ensure all log records have been written to disk
        //   local_state.model_state_logger = nullptr;
        //   std::scoped_lock<std::mutex> lck(model_state_logger_mtx);
        //   utils::concatenate_files<true>(this->path_to_model_state_log_file,
        //                                  tmp_model_internal_state_log_path);
        // }
        spdlog::debug(FMT_STRING("[W{}]: all tasks have been processed: returning!"), tid);
        return;
      }

      // Loop over new tasks
      for (const auto& task : task_buff) {
        if (!this->_ctx) {
          return;
        }
        local_state = task;            // Set simulation state based on task data
        local_state.resize_buffers();  // This resizes buffers based on the nlefs to be simulated
        local_state.reset_buffers();   // Clear all buffers

        if (task.cell_id == 0) {
          this->print_status_update(task);
        }

        // When simulating using the target contact density as stopping criterion across a large
        // number of cells, if the total number of contacts to be simulated is low (e.g. because a
        // chromosome is small, the resolution is low, or the target contact density is very small),
        // some of the tasks may have num_target_contacts = 0.
        // These tasks can be safely skipped as they won't have any effect on the output contact
        // matrix.
        if (MODLE_LIKELY(task.num_target_epochs != (std::numeric_limits<usize>::max)() ||
                         task.num_target_contacts != 0)) {
          local_state.status = State::Status::RUNNING;
          spdlog::debug(FMT_STRING("[W{}]: begin processing task {} ({} cell #{}, {})..."), tid,
                        task.id, *task.interval, task.cell_id, format_rand_eng(task.rand_eng));
          Simulation::simulate_one_cell(tid, local_state);
          spdlog::debug(
              FMT_STRING("[W{}]: finished processing task {} ({} cell #{}, {}): collected {} "
                         "interactions throughout {} epochs ({} burnin epochs)"),
              tid, local_state.id, *local_state.interval, local_state.cell_id,
              format_rand_eng(task.rand_eng), local_state.num_contacts, local_state.epoch,
              local_state.num_burnin_epochs);
        }

        local_state.status = State::Status::COMPLETED;
        // Update progress for the current chrom
        while (!this->_ctx.try_enqueue<Task::Status::COMPLETED>(static_cast<Task>(local_state),
                                                                ptok)) {
          if (!this->_ctx) {
            return;
          }
        }
      }
    }
  } catch (const std::exception& e) {
    const auto excp = std::runtime_error(
        fmt::format(FMT_STRING("Exception raised in worker thread {}:\n   {}\n   {}"), tid,
                    local_state, e.what()));
    this->_ctx.set_exception_worker(tid, std::make_exception_ptr(excp));
  } catch (...) {
    this->_ctx.set_exception_worker(tid, std::current_exception());
  }
}

}  // namespace modle
