// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <moodycamel/blockingconcurrentqueue.h>

#include <BS_thread_pool.hpp>
#include <atomic>
#include <chrono>
#include <exception>
#include <optional>
#include <vector>

#include "modle/common/common.hpp"

namespace modle {

template <typename Task>
class ContextManager {
  using Status = typename Task::Status;
  using QueueT = moodycamel::BlockingConcurrentQueue<Task>;
  QueueT _pending{};
  QueueT _finished{};
  BS::thread_pool _worker_tpool;
  BS::thread_pool _io_tpool;
  std::atomic<usize> _num_in{};
  std::atomic<usize> _num_out{};

  std::vector<std::exception_ptr> _exception_ptrs_io{};
  std::vector<std::exception_ptr> _exception_ptrs_workers{};
  std::atomic<bool> _exception_thrown{false};
  std::atomic<bool> _queue_closed{false};
  std::atomic<bool> _shutdown_requested{false};

 public:
  ContextManager() = delete;
  explicit ContextManager(usize num_worker_threads, usize num_io_threads = 1);

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] bool shutdown_signal_sent() const noexcept;

  template <Status s>
  [[nodiscard]] bool try_enqueue(Task&& t, moodycamel::ProducerToken& ptok);
  template <Status s>
  [[nodiscard]] std::optional<Task> try_dequeue(moodycamel::ConsumerToken& ctok);
  template <Status s, typename TimeT = std::chrono::milliseconds>
  [[nodiscard]] std::optional<Task> wait_dequeue(moodycamel::ConsumerToken& ctok,
                                                 TimeT timeout = std::chrono::milliseconds(100));
  template <Status s>
  void try_dequeue(moodycamel::ConsumerToken& ctok, std::vector<Task>& buff);
  template <Status s, typename TimeT = std::chrono::milliseconds>
  void wait_dequeue(moodycamel::ConsumerToken& ctok, std::vector<Task>& buff,
                    TimeT timeout = std::chrono::milliseconds(100));

  void set_exception_worker(usize idx, std::exception_ptr e);
  void set_exception_io(usize idx, std::exception_ptr e);
  void set_exception_main(std::exception_ptr e);
  [[nodiscard]] bool exception_thrown() const noexcept;

  [[nodiscard]] usize num_submitted() const noexcept;
  [[nodiscard]] usize num_completed() const noexcept;

  template <Status s>
  [[nodiscard]] moodycamel::ProducerToken register_producer();
  template <Status s>
  [[nodiscard]] moodycamel::ConsumerToken register_consumer();
  void check_exceptions();
  [[nodiscard]] usize num_threads() const noexcept;
  [[nodiscard]] usize num_worker_threads() const noexcept;
  [[nodiscard]] usize num_io_threads() const noexcept;

  template <typename TaskLambda>
  void spawn_worker_thread(TaskLambda lambda);
  template <typename TaskLambda>
  void spawn_io_thread(TaskLambda lambda);
  void shutdown();

 private:
  void close_queue() noexcept;
  [[nodiscard]] bool queue_is_closed() const noexcept;
  [[noreturn]] void rethrow_exceptions() const;
};

}  // namespace modle

#include "../../context_manager_impl.hpp"