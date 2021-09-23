// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "modle/simulation.hpp"
// clang-format on

#include <absl/container/fixed_array.h>          // for FixedArray
#include <absl/types/span.h>                     // for Span, MakeConstSpan
#include <fmt/format.h>                          // for format, make_format_args, vformat_to
#include <moodycamel/blockingconcurrentqueue.h>  // for BlockingConcurrentQueue
#include <moodycamel/concurrentqueue.h>          // for ConsumerToken, ProducerToken
#include <spdlog/spdlog.h>                       // for info

#include <algorithm>                        // for max, min
#include <array>                            // for array, array<>::value_type
#include <atomic>                           // for atomic
#include <boost/filesystem/operations.hpp>  // for exists
#include <cassert>                          // for assert
#include <chrono>                           // for microseconds, milliseconds
#include <cstddef>                          // for size_t
#include <cstdint>                          // for uint64_t
#include <exception>                        // for exception_ptr, exception, current_exception
#include <iterator>                         // for move_iterator, make_move_iterator
#include <mutex>                            // for mutex, scoped_lock
#include <stdexcept>                        // for runtime_error
#include <string>                           // for string
#include <thread>                           // for sleep_for
#include <thread_pool/thread_pool.hpp>      // for thread_pool
#include <vector>                           // for vector

#include "modle/common/common.hpp"  // for contacts_t, bp_t
#include "modle/compressed_io.hpp"  // for Reader, Writer
#include "modle/contacts.hpp"       // for ContactMatrix

namespace modle {

void Simulation::run_replay() {
  assert(boost::filesystem::exists(this->path_to_task_file));  // NOLINT
  compressed_io::Reader task_reader(this->path_to_task_file);

  const size_t task_batch_size_enq = 32;  // NOLINTNEXTLINE
  moodycamel::BlockingConcurrentQueue<TaskPW> task_queue(this->nthreads * 2, 1, 0);
  moodycamel::ProducerToken ptok(task_queue);
  std::array<TaskPW, task_batch_size_enq> tasks;

  std::mutex cooler_mutex;

  try {
    this->_tpool.reset(this->nthreads);
    for (uint64_t tid = 0; tid < this->nthreads; ++tid) {  // Start simulation threads
      this->_tpool.push_task([&, tid]() { this->replay_worker(tid, task_queue, cooler_mutex); });
    }

    const auto task_filter = import_task_filter(this->path_to_task_filter_file);

    std::string task_definition;
    size_t num_tasks = 0;
    while (task_reader.getline(task_definition)) {
      if (!this->ok()) {
        this->handle_exceptions();
      }
      if (num_tasks == tasks.size()) {  // Enqueue a batch of tasks
        auto sleep_us = 100;            // NOLINT(readability-magic-numbers)
        while (
            !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), num_tasks)) {
          if (!this->ok()) {
            this->handle_exceptions();
          }  // NOLINTNEXTLINE(readability-magic-numbers)
          sleep_us = std::min(100000, sleep_us * 2);
          std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
        }
        num_tasks = 0;
      }  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
      const auto& t = (tasks[num_tasks] = TaskPW::from_string(task_definition, this->_genome));
      if (task_filter.empty() || task_filter.contains(t.id)) {
        ++num_tasks;
      }
    }

    if (num_tasks != 0) {
      while (
          !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), num_tasks)) {
        if (!this->ok()) {
          this->handle_exceptions();
        }  // NOLINTNEXTLINE(readability-magic-numbers)
        std::this_thread::sleep_for(std::chrono::microseconds(100));
      }
    }

    this->_end_of_simulation = true;
    this->_tpool.wait_for_tasks();     // Wait on simulate_worker threads
    assert(!this->_exception_thrown);  // NOLINT
  } catch (...) {
    this->_exception_thrown = true;
    this->_tpool.paused = true;
    this->_tpool.wait_for_tasks();
    throw;
  }
}

void Simulation::replay_worker(const uint64_t tid,
                               moodycamel::BlockingConcurrentQueue<Simulation::TaskPW>& task_queue,
                               std::mutex& cooler_mutex, const size_t task_batch_size) {
  spdlog::info(FMT_STRING("Spawning simulation thread {}..."), tid);
  moodycamel::ConsumerToken ctok(task_queue);

  absl::FixedArray<TaskPW> task_buff(task_batch_size);  // Tasks are dequeue in batch.

  Simulation::State local_state{};
  local_state.contacts = std::make_shared<ContactMatrix<contacts_t>>();
  compressed_io::Writer null_stream{};

  try {
    while (this->ok()) {  // Try to dequeue a batch of tasks
      const auto avail_tasks = task_queue.wait_dequeue_bulk_timed(
          ctok, task_buff.begin(), task_buff.size(), std::chrono::milliseconds(10));
      // Check whether dequeue operation timed-out before any task became available
      if (avail_tasks == 0) {
        if (this->_end_of_simulation) {
          // Reached end of simulation (i.e. all tasks have been processed)
          return;
        }
        // Keep waiting until one or more tasks become available
        continue;
      }

      // Loop over new tasks
      for (const auto& task : absl::MakeConstSpan(task_buff.data(), avail_tasks)) {
        if (!this->ok()) {
          return;
        }
        assert(!task.barriers.empty());  // NOLINT
        assert(!task.feats1.empty());    // NOLINT
        assert(!task.feats2.empty());    // NOLINT
        assert(local_state.contacts);    // NOLINT

        local_state = task;  // Set simulation local_state based on task data
        local_state.contacts->unsafe_resize(local_state.window_end - local_state.window_start,
                                            this->diagonal_width, this->bin_size);

        Simulation::simulate_window(local_state, null_stream, cooler_mutex, true);
      }
    }
  } catch (const std::exception& e) {
    std::scoped_lock<std::mutex> l(this->_exceptions_mutex);
    this->_exceptions.emplace_back(std::make_exception_ptr(
        std::runtime_error(fmt::format(FMT_STRING("Detected an error in worker thread {}:\n{}\n{}"),
                                       tid, local_state.to_string(), e.what()))));
    this->_exception_thrown = true;

  } catch (...) {
    std::scoped_lock<std::mutex> l(this->_exceptions_mutex);
    this->_exceptions.emplace_back(std::current_exception());
    this->_exception_thrown = true;
  }
}
}  // namespace modle
