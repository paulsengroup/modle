// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <moodycamel/blockingconcurrentqueue.h>

#include <BS_thread_pool.hpp>
#include <atomic>
#include <chrono>
#include <exception>
#include <filesystem>
#include <future>
#include <optional>
#include <string_view>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/state_logger.hpp"

namespace modle {

template <typename Task>
class ContextManager {
  using Status = typename Task::Status;
  using QueueT = moodycamel::BlockingConcurrentQueue<Task>;
  std::unique_ptr<StateLoggerAggregator> _state_logger_ptr{};
  QueueT _pending{};
  QueueT _finished{};
  BS::light_thread_pool _worker_tpool;
  BS::light_thread_pool _io_tpool;
  std::atomic<std::size_t> _num_in{};
  std::atomic<std::size_t> _num_out{};

  mutable std::vector<std::future<void>> _futures{};
  mutable std::atomic<bool> _exception_thrown{false};
  std::atomic<bool> _queue_closed{false};
  std::atomic<bool> _shutdown_requested{false};

 public:
  ContextManager() = delete;
  explicit ContextManager(std::size_t num_worker_threads, std::size_t num_io_threads = 1);

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] bool shutdown_signal_sent() const noexcept;

  template <Status s>
  [[nodiscard]] bool try_enqueue_task(Task t, moodycamel::ProducerToken& ptok);
  template <Status s>
  [[nodiscard]] std::optional<Task> try_dequeue_task(moodycamel::ConsumerToken& ctok);
  template <Status s, typename TimeT = std::chrono::milliseconds>
  [[nodiscard]] std::optional<Task> wait_dequeue_task(
      moodycamel::ConsumerToken& ctok, TimeT timeout = std::chrono::milliseconds(100));
  template <Status s>
  void try_dequeue_tasks(moodycamel::ConsumerToken& ctok, std::vector<Task>& buff);
  template <Status s, typename TimeT = std::chrono::milliseconds>
  void wait_dequeue_tasks(moodycamel::ConsumerToken& ctok, std::vector<Task>& buff,
                          TimeT timeout = std::chrono::milliseconds(100));

  template <typename Exception>
  [[noreturn]] void throw_exception(Exception exception) const;
  void set_exception_main(std::exception_ptr e);

  [[nodiscard]] std::size_t num_tasks_submitted() const noexcept;
  [[nodiscard]] std::size_t num_tasks_completed() const noexcept;

  template <Status s>
  [[nodiscard]] moodycamel::ProducerToken register_producer();
  template <Status s>
  [[nodiscard]] moodycamel::ConsumerToken register_consumer();
  void check_exceptions();
  [[nodiscard]] std::size_t num_threads() const noexcept;
  [[nodiscard]] std::size_t num_worker_threads() const noexcept;
  [[nodiscard]] std::size_t num_io_threads() const noexcept;

  template <typename TaskLambda>
  void spawn_worker_thread(TaskLambda lambda);
  template <typename TaskLambda>
  void spawn_io_thread(TaskLambda lambda);
  void shutdown();

  void init_model_state_logger(std::filesystem::path path, std::string_view header);
  void append_to_model_state_log(const std::filesystem::path& path_,
                                 bool remove_file_after_append = false);

 private:
  template <Status s>
  auto get_queue() noexcept -> QueueT&;
  void close_queue() noexcept;
  [[nodiscard]] bool queue_is_closed() const noexcept;
  [[noreturn]] void rethrow_exceptions() const;
};

}  // namespace modle

#include "../../context_manager_impl.hpp"
