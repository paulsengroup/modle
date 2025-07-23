// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "modle/simulation.hpp"
// clang-format on
#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include <xxhash.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <deque>
#include <exception>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <thread>  // IWYU pragma: keep for sleep_for
#include <utility>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/genome.hpp"

namespace modle {

void Simulation::spawn_worker_threads(const usize batch_size) {
  assert(_ctx.num_worker_threads() != 0);
  assert(batch_size != 0);
  for (usize tid = 0; tid < _ctx.num_worker_threads(); ++tid) {
    _ctx.spawn_worker_thread([this, batch_size, tid]() { simulate_worker(tid, batch_size); });
  }
}

void Simulation::spawn_io_threads() {
  _ctx.spawn_io_thread([this]() { simulate_io(); });
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
      if (c().log_model_internal_state) {
        std::filesystem::remove(c().path_to_model_state_log_file);
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
  auto ptok_pending = _ctx.register_producer<Task::Status::PENDING>();
  auto ptok_finished = _ctx.register_producer<Task::Status::COMPLETED>();

  // std::mutex model_state_logger_mtx;  // Protect rw access to the log file located at
  // path_to_model_state_log_file
  if (c().log_model_internal_state && !c().skip_output) {
    _ctx.init_model_state_logger(c().path_to_model_state_log_file, model_internal_state_log_header);
  }

  try {
    spawn_io_threads();
    spawn_worker_threads(task_batch_size);

    // The remaining code submits simulation tasks to the queue. Then it waits until all the tasks
    // have been completed and contacts have been written to disk

    usize taskid = 0;
    const auto max_sleep_time =
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::milliseconds(250));

    std::unique_ptr<XXH3_state_t, utils::XXH3_Deleter> xxh_state{XXH3_createState()};

    // Loop over chromosomes
    for (auto& interval : _genome) {
      _ctx.check_exceptions();
      auto sleep_time = std::chrono::microseconds(50);

      auto rand_eng = random::PRNG(interval.hash(*xxh_state, c().seed));

      // Don't bother simulating intervals without barriers
      if (!c().simulate_chromosomes_wo_barriers && interval.num_barriers() == 0) {
        SPDLOG_INFO("{} has 0 barriers... SKIPPING!", interval);
        Task t{};
        t.interval = &interval;
        for (usize cellid = 0; cellid < c().num_cells; ++cellid) {
          while (!_ctx.try_enqueue_task<Task::Status::COMPLETED>(t, ptok_finished)) {
            _ctx.check_exceptions();
            sleep_time = std::min(max_sleep_time, sleep_time * 2);
            std::this_thread::sleep_for(sleep_time);
          }
          rand_eng.jump();
        }
        continue;
      }

      // Compute # of LEFs to be simulated based on interval.sizes
      const auto nlefs = compute_num_lefs(interval.size());

      const auto npixels = interval.npixels();
      const auto tot_target_contacts =
          static_cast<usize>(std::round(static_cast<double>(npixels) * c().target_contact_density));

      const auto target_contacts_per_cell =
          (tot_target_contacts + c().num_cells - 1) / c().num_cells;

      usize tot_target_contacts_rolling_count = 0;
      for (usize cellid = 0; cellid < c().num_cells; ++cellid) {
        // This is needed to not overshoot the target contact density
        const auto num_target_contacts = std::min(
            target_contacts_per_cell, tot_target_contacts - tot_target_contacts_rolling_count);
        tot_target_contacts_rolling_count += num_target_contacts;

        Task t{taskid++,
               &interval,
               cellid,
               c().target_simulation_epochs,
               num_target_contacts,
               nlefs,
               rand_eng,
               Task::Status::PENDING};
        SPDLOG_DEBUG("[main]: submitting task #{} ({} cell #{})...", t.id, *t.interval, t.cell_id);

        while (!_ctx.try_enqueue_task<Task::Status::PENDING>(t, ptok_pending)) {
          _ctx.check_exceptions();
          sleep_time = std::min(max_sleep_time, sleep_time * 2);
          std::this_thread::sleep_for(sleep_time);
        }
        rand_eng.jump();
      }
    }

    _ctx.shutdown();
    assert(_ctx.num_tasks_submitted() == _ctx.num_tasks_completed());
    assert(!!_ctx);
  } catch (...) {
    _ctx.set_exception_main(std::current_exception());
    _ctx.check_exceptions();
    throw;
  }
}

[[maybe_unused]] static std::string format_rand_eng(const random::PRNG_t& rand_eng) {
  const auto state = rand_eng.serialize();
  static_assert(state.size() == 4);
  static_assert(sizeof(state[0]) == sizeof(u64));
  return fmt::format("0x{:016x}{:016x}{:016x}{:016x}", state[0], state[1], state[2], state[3]);
}

[[nodiscard]] static std::unique_ptr<compressed_io::Writer> init_local_model_state_logger(
    const Config& c, const u64 task_id) {
  if (!c.log_model_internal_state) {
    return nullptr;
  }

  auto path = c.path_to_model_state_log_file;
  path.replace_extension(fmt::format("{}{}", task_id, path.extension().string()));
  return std::make_unique<compressed_io::Writer>(path);
}

void Simulation::simulate_worker(const u64 tid, const usize task_batch_size) {
  SPDLOG_INFO("spawning worker thread W{}...", tid);
  // This state object owns all the buffers and PRNG + seed required in order to simulate loop
  // extrusion for a single cell. Buffers are allocated once and resized, cleared and reused
  // throughout the simulation
  Simulation::State local_state{};

  try {
    auto ctok = _ctx.register_consumer<Task::Status::PENDING>();
    auto ptok = _ctx.register_producer<Task::Status::COMPLETED>();
    std::vector<Task> task_buff(task_batch_size);  // Tasks are dequeued in batch.
    // This is to reduce contention when accessing the queue

    while (!!_ctx) {
      _ctx.wait_dequeue_tasks<Task::Status::PENDING>(ctok, task_buff);
      if (task_buff.empty() && !_ctx.shutdown_signal_sent()) {
        SPDLOG_DEBUG("[W{}]: no tasks available: sleeping for a bit...", tid);
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        continue;
      }
      if (task_buff.empty() && _ctx.shutdown_signal_sent()) {
        SPDLOG_DEBUG("[W{}]: all tasks have been processed: returning!", tid);
        return;
      }

      // Loop over new tasks
      for (const auto& task : task_buff) {
        if (!_ctx) {
          return;
        }
        local_state = task;            // Set simulation state based on task data
        local_state.resize_buffers();  // This resizes buffers based on the nlefs to be simulated
        local_state.reset_buffers();   // Clear all buffers

        if (task.cell_id == 0) {
          print_status_update(task);
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
          SPDLOG_DEBUG("[W{}]: begin processing task {} ({} cell #{}, {})...", tid, task.id,
                       *task.interval, task.cell_id, format_rand_eng(task.rand_eng));
          local_state.model_state_logger = init_local_model_state_logger(c(), local_state.id);
          simulate_one_cell(tid, local_state);
          if (local_state.model_state_logger) {
            const auto path = local_state.model_state_logger->path();
            local_state.model_state_logger = nullptr;  // Flush data to disk
            _ctx.append_to_model_state_log(path, true);
          }
          SPDLOG_DEBUG(
              "[W{}]: finished processing task {} ({} cell #{}, {}): collected {} "
              "interactions throughout {} epochs ({} burnin epochs)",
              tid, local_state.id, *local_state.interval, local_state.cell_id,
              format_rand_eng(task.rand_eng), local_state.num_contacts, local_state.epoch,
              local_state.num_burnin_epochs);
        }

        local_state.status = Task::Status::COMPLETED;
        // Update progress for the current chrom
        while (!_ctx.try_enqueue_task<Task::Status::COMPLETED>(
            static_cast<const Task&>(local_state), ptok)) {
          if (!_ctx) {
            return;
          }
        }
      }
    }
  } catch (const std::exception& e) {
    _ctx.throw_exception(std::runtime_error(fmt::format(
        "Exception raised in worker thread {}:\n   {}\n   {}", tid, local_state, e.what())));
  } catch (...) {
    _ctx.throw_exception(
        std::runtime_error(fmt::format("Unhandled exception raised in worker thread {}!", tid)));
  }
}

}  // namespace modle
