// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
#include <cassert>
#include <exception>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/compressed_io/compressed_io.hpp"

namespace modle {

template <typename Task>
inline ContextManager<Task>::ContextManager(usize num_worker_threads, usize num_io_threads)
    : _pending(2 * num_worker_threads, 1, 0),
      _finished(2 * num_worker_threads, num_worker_threads + 1, 0),
      _worker_tpool(num_worker_threads),
      _io_tpool(num_io_threads) {}

template <typename Task>
inline ContextManager<Task>::operator bool() const noexcept {
  return !this->_exception_thrown.load();
}

template <typename Task>
template <typename ContextManager<Task>::Status s>
inline auto ContextManager<Task>::get_queue() noexcept -> QueueT& {
  static_assert(s == Status::PENDING || s == Status::COMPLETED);
  if constexpr (s == Status::PENDING) {
    return this->_pending;
  } else {
    return this->_finished;
  }
}

template <typename Task>
inline void ContextManager<Task>::close_queue() noexcept {
  this->_queue_closed = true;
}

template <typename Task>
inline bool ContextManager<Task>::queue_is_closed() const noexcept {
  return this->_queue_closed.load();
}

template <typename Task>
inline bool ContextManager<Task>::shutdown_signal_sent() const noexcept {
  return this->_shutdown_requested.load();
}

template <typename Task>
template <typename ContextManager<Task>::Status s>
inline bool ContextManager<Task>::try_enqueue_task(Task t, moodycamel::ProducerToken& ptok) {
  const auto successful = this->get_queue<s>().try_enqueue(ptok, std::move(t));

  this->_num_in += successful;
  return successful;
}

template <typename Task>
template <typename ContextManager<Task>::Status s>
inline std::optional<Task> ContextManager<Task>::try_dequeue_task(moodycamel::ConsumerToken& ctok) {
  Task t;  // NOLINT
  if (this->get_queue<s>().try_dequeue(ctok, t)) {
    ++this->_num_out;
    return {t};
  }
  return std::optional<Task>{};
}

template <typename Task>
template <typename ContextManager<Task>::Status s, typename TimeT>
inline std::optional<Task> ContextManager<Task>::wait_dequeue_task(moodycamel::ConsumerToken& ctok,
                                                                   TimeT timeout) {
  Task t;  // NOLINT
  if (this->get_queue<s>().wait_dequeue_timed(ctok, t, timeout)) {
    ++this->_num_out;
    return {t};
  }
  return std::optional<Task>{};
}

template <typename Task>
template <typename ContextManager<Task>::Status s>
inline void ContextManager<Task>::try_dequeue_tasks(moodycamel::ConsumerToken& ctok,
                                                    std::vector<Task>& buff) {
  assert(buff.capacity() != 0);
  buff.clear();
  const auto num_tasks =
      this->get_queue<s>().try_dequeue_bulk(ctok, std::back_inserter(buff), buff.capacity());
  assert(num_tasks == buff.size());
  this->_num_out += num_tasks;
}

template <typename Task>
template <typename ContextManager<Task>::Status s, typename TimeT>
inline void ContextManager<Task>::wait_dequeue_tasks(moodycamel::ConsumerToken& ctok,
                                                     std::vector<Task>& buff, TimeT timeout) {
  static_assert(s == Status::PENDING || s == Status::COMPLETED);
  assert(buff.capacity() != 0);
  buff.clear();
  const auto num_tasks = this->get_queue<s>().wait_dequeue_bulk_timed(
      ctok, std::back_inserter(buff), buff.capacity(), timeout);
  assert(num_tasks == buff.size());
  this->_num_out += num_tasks;
}

template <typename Task>
template <typename Exception>
inline void ContextManager<Task>::throw_exception(Exception except) const {
  this->_exception_thrown = true;
  throw std::runtime_error(except);
}

template <typename Task>
inline void ContextManager<Task>::set_exception_main(std::exception_ptr e) {
  if (e != nullptr) {
    this->_exception_thrown = true;
  }
}

template <typename Task>
inline usize ContextManager<Task>::num_tasks_submitted() const noexcept {
  return this->_num_in.load();
}
template <typename Task>
inline usize ContextManager<Task>::num_tasks_completed() const noexcept {
  return this->_num_out.load();
}

template <typename Task>
template <typename ContextManager<Task>::Status s>
inline moodycamel::ProducerToken ContextManager<Task>::register_producer() {
  return moodycamel::ProducerToken(this->get_queue<s>());
}

template <typename Task>
template <typename ContextManager<Task>::Status s>
inline moodycamel::ConsumerToken ContextManager<Task>::register_consumer() {
  return moodycamel::ConsumerToken(this->get_queue<s>());
}

template <typename Task>
inline void ContextManager<Task>::rethrow_exceptions() const {
  std::vector<std::string> messages;
  for (auto& fut : this->_futures) {
    try {
      fut.get();
    } catch (const std::exception& e) {
      messages.emplace_back(e.what());
    } catch (...) {
      messages.emplace_back(
          "An unhandled exception was caught! This should never happen! If you see "
          "this message, please file an issue on GitHub.");
    }
  }

  throw std::runtime_error(fmt::format(
      FMT_STRING("the following error(s) occurred while simulating loop extrusion:\n - {}"),
      fmt::join(messages, "\n - ")));
}

template <typename Task>
inline void ContextManager<Task>::check_exceptions() {
  if (MODLE_UNLIKELY(this->_exception_thrown.load())) {
    SPDLOG_ERROR("MoDLE encountered an exception. Shutting down worker threads...");
    this->shutdown();
  }
}

template <typename Task>
inline usize ContextManager<Task>::num_threads() const noexcept {
  return this->num_worker_threads() + this->num_io_threads();
}
template <typename Task>
inline usize ContextManager<Task>::num_worker_threads() const noexcept {
  return utils::conditional_static_cast<usize>(this->_worker_tpool.get_thread_count());
}
template <typename Task>
inline usize ContextManager<Task>::num_io_threads() const noexcept {
  return utils::conditional_static_cast<usize>(this->_io_tpool.get_thread_count());
}

template <typename Task>
template <typename TaskLambda>
inline void ContextManager<Task>::spawn_worker_thread(TaskLambda lambda) {
  this->_futures.emplace_back(this->_worker_tpool.submit_task(lambda));
}

template <typename Task>
template <typename TaskLambda>
inline void ContextManager<Task>::spawn_io_thread(TaskLambda lambda) {
  this->_futures.emplace_back(this->_io_tpool.submit_task(lambda));
}

template <typename Task>
inline void ContextManager<Task>::shutdown() {
  this->_shutdown_requested = true;
  this->close_queue();

  SPDLOG_DEBUG("waiting for worker threads to return...");
  this->_worker_tpool.wait();
  SPDLOG_DEBUG("waiting for io threads to return...");
  this->_io_tpool.wait();

  SPDLOG_DEBUG("all background threads returned! Checking if any exception have been raised...");
  if (this->_exception_thrown) {
    this->rethrow_exceptions();
  } else {
    for (auto& fut : this->_futures) {
      fut.get();
    }
  }

  const auto num_completed = this->num_tasks_completed();
  const auto num_submitted = this->num_tasks_submitted();
  if (num_completed != num_submitted) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("ContextManager: not all tasks have been processed! {} tasks have "
                               "been submitted but only {} successfully completed.\n"
                               "This should never happen!\n"
                               "Please file an issue on GitHub."),
                    num_submitted, num_completed));
  }
  SPDLOG_DEBUG("context manager shutdown was successful.");
}

template <typename Task>
inline void ContextManager<Task>::init_model_state_logger(std::filesystem::path path,
                                                          std::string_view header) {
  assert(!this->_state_logger_ptr);
  assert(!std::filesystem::exists(path));
  SPDLOG_DEBUG(FMT_STRING("[io] initializing model state logger at {}..."), path);
  compressed_io::Writer(path).write(header);
  this->_state_logger_ptr = std::make_unique<StateLoggerAggregator>(std::move(path));
}

template <typename Task>
inline void ContextManager<Task>::append_to_model_state_log(const std::filesystem::path& path,
                                                            bool remove_file_after_append) {
  assert(!!this->_state_logger_ptr);
  SPDLOG_DEBUG(FMT_STRING("[io] appending {} to log file at {}..."), path,
               this->_state_logger_ptr->path());
  this->_state_logger_ptr->append(path, remove_file_after_append);
}

}  // namespace modle
