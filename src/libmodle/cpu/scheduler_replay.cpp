// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "modle/simulation.hpp"
// clang-format on

#include <absl/container/fixed_array.h>          // for FixedArray
#include <absl/container/flat_hash_set.h>        // for flat_hash_set, BitMask
#include <absl/types/span.h>                     // for Span, MakeConstSpan
#include <fmt/format.h>                          // for format, make_format_args, vformat_to
#include <moodycamel/blockingconcurrentqueue.h>  // for BlockingConcurrentQueue
#include <moodycamel/concurrentqueue.h>          // for ConsumerToken, ProducerToken
#include <spdlog/spdlog.h>                       // for info

#include <BS_thread_pool.hpp>  // for BS::thread_pool
#include <algorithm>           // for max, min
#include <array>               // for array, array<>::value_type
#include <atomic>              // for atomic
#include <cassert>             // for assert
#include <chrono>              // for microseconds, milliseconds
#include <exception>           // for exception_ptr, exception, current_exception
#include <filesystem>          // for exists
#include <iterator>            // for move_iterator, make_move_iterator
#include <memory>              // for shared_ptr, make_shared, __shared...
#include <mutex>               // for mutex, scoped_lock
#include <stdexcept>           // for runtime_error
#include <string>              // for string
#include <thread>              // IWYU pragma: keep for sleep_for
#include <vector>              // for vector

#include "modle/common/common.hpp"                // for bp_t, contacts_t, u64
#include "modle/compressed_io/compressed_io.hpp"  // for Reader, Writer
#include "modle/contact_matrix_dense.hpp"         // for ContactMatrixDense

namespace modle {

void Simulation::run_replay() {
  assert(std::filesystem::exists(this->path_to_task_file));
  compressed_io::Reader task_reader(this->path_to_task_file);

  const usize task_batch_size_enq = 32;
  moodycamel::BlockingConcurrentQueue<TaskPW> task_queue(this->nthreads * 2, 1, 0);
  moodycamel::ProducerToken ptok(task_queue);
  std::array<TaskPW, task_batch_size_enq> tasks;

  std::mutex cooler_mutex;

  try {
    this->_tpool.reset(utils::conditional_static_cast<BS::concurrency_t>(this->nthreads));
    for (u64 tid = 0; tid < this->nthreads; ++tid) {  // Start simulation threads
      this->_tpool.push_task([&, tid]() { this->replay_worker(tid, task_queue, cooler_mutex); });
    }

    const auto task_filter = import_task_filter(this->path_to_task_filter_file);

    std::string task_definition;
    usize num_tasks = 0;
    while (task_reader.getline(task_definition)) {
      if (!this->ok()) {
        this->handle_exceptions();
      }
      if (num_tasks == tasks.size()) {  // Enqueue a batch of tasks
        auto sleep_us = 100;
        while (
            !task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), num_tasks)) {
          if (!this->ok()) {
            this->handle_exceptions();
          }
          sleep_us = std::min(100000, sleep_us * 2);
          std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
        }
        num_tasks = 0;
      }
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
        }
        std::this_thread::sleep_for(std::chrono::microseconds(100));
      }
    }

    this->_end_of_simulation = true;
    this->_tpool.wait_for_tasks();  // Wait on simulate_worker threads
    assert(!this->_exception_thrown);
  } catch (...) {
    this->_exception_thrown = true;
    this->_tpool.pause();
    this->_tpool.wait_for_tasks();
    throw;
  }
}

void Simulation::replay_worker(const u64 tid,
                               moodycamel::BlockingConcurrentQueue<Simulation::TaskPW>& task_queue,
                               std::mutex& cooler_mtx, const usize task_batch_size) {
  spdlog::info(FMT_STRING("Spawning simulation thread {}..."), tid);
  moodycamel::ConsumerToken ctok(task_queue);

  absl::FixedArray<TaskPW> task_buff(task_batch_size);  // Tasks are dequeue in batch.

  Simulation::State local_state{};
  local_state.contacts = std::make_shared<ContactMatrixDense<contacts_t>>();
  compressed_io::Writer null_stream{};

  try {
    while (this->ok()) {  // Try to dequeue a batch of tasks
      const auto avail_tasks = this->consume_tasks_blocking(task_queue, ctok, task_buff);
      if (avail_tasks == 0) {
        assert(this->_end_of_simulation);
        // Reached end of simulation (i.e. all tasks have been processed)
        return;
      }

      // Loop over new tasks
      for (const auto& task : absl::MakeConstSpan(task_buff.data(), avail_tasks)) {
        if (!this->ok()) {
          return;
        }
        assert(!task.barriers.empty());
        assert(!task.feats1.empty());
        assert(!task.feats2.empty());
        assert(local_state.contacts);

        local_state = task;  // Set simulation local_state based on task data
        local_state.contacts->unsafe_resize(local_state.window_end - local_state.window_start,
                                            this->diagonal_width, this->bin_size);

        Simulation::simulate_window(local_state, null_stream, cooler_mtx, true);
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
