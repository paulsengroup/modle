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
      _worker_tpool(utils::conditional_static_cast<BS::concurrency_t>(num_worker_threads)),
      _io_tpool(utils::conditional_static_cast<BS::concurrency_t>(num_io_threads)),
      _exception_ptrs_io(num_io_threads),
      _exception_ptrs_workers(num_worker_threads) {}

template <typename Task>
inline ContextManager<Task>::operator bool() const noexcept {
  return !this->exception_thrown();
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
inline bool ContextManager<Task>::try_enqueue_task(Task&& t, moodycamel::ProducerToken& ptok) {
  static_assert(s == Status::PENDING || s == Status::COMPLETED);
  bool successful{};
  if constexpr (s == Status::PENDING) {
    successful = this->_pending.try_enqueue(ptok, std::move(t));

  } else {
    successful = this->_finished.try_enqueue(ptok, std::move(t));
  }

  this->_num_in += successful;
  return successful;
}

template <typename Task>
template <typename ContextManager<Task>::Status s>
inline std::optional<Task> ContextManager<Task>::try_dequeue_task(moodycamel::ConsumerToken& ctok) {
  static_assert(s == Status::PENDING || s == Status::COMPLETED);
  Task t;  // NOLINT
  const auto successful = [&]() {
    if constexpr (s == Status::PENDING) {
      return this->_pending.try_dequeue(ctok, t);
    } else {
      return this->_finished.try_dequeue(ctok, t);
    }
  }();

  if (successful) {
    ++this->_num_out;
    return {t};
  }
  return std::optional<Task>{};
}

template <typename Task>
template <typename ContextManager<Task>::Status s, typename TimeT>
inline std::optional<Task> ContextManager<Task>::wait_dequeue_task(moodycamel::ConsumerToken& ctok,
                                                                   TimeT timeout) {
  static_assert(s == Status::PENDING || s == Status::COMPLETED);
  Task t;  // NOLINT
  const auto successful = [&]() {
    if constexpr (s == Status::PENDING) {
      return this->_pending.wait_dequeue_timed(ctok, t, timeout);
    } else {
      return this->_finished.wait_dequeue_timed(ctok, t, timeout);
    }
  }();

  if (successful) {
    ++this->_num_out;
    return {t};
  }
  return std::optional<Task>{};
}

template <typename Task>
template <typename ContextManager<Task>::Status s>
inline void ContextManager<Task>::try_dequeue_tasks(moodycamel::ConsumerToken& ctok,
                                                    std::vector<Task>& buff) {
  static_assert(s == Status::PENDING || s == Status::COMPLETED);
  assert(buff.capacity() != 0);
  buff.clear();
  const auto num_tasks = [&]() {
    if constexpr (s == Status::PENDING) {
      return this->_pending.try_dequeue_bulk(ctok, std::back_inserter(buff.begin()),
                                             buff.capacity());
    } else {
      return this->_finished.try_dequeue_bulk(ctok, std::back_inserter(buff.begin()),
                                              buff.capacity());
    }
  }();

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
  [[maybe_unused]] const auto num_tasks = [&]() {
    if constexpr (s == Status::PENDING) {
      return this->_pending.wait_dequeue_bulk_timed(ctok, std::back_inserter(buff), buff.capacity(),
                                                    timeout);
    } else {
      return this->_finished.wait_dequeue_bulk_timed(ctok, std::back_inserter(buff),
                                                     buff.capacity(), timeout);
    }
  }();

  assert(num_tasks == buff.size());
  this->_num_out += num_tasks;
}
template <typename Task>
inline void ContextManager<Task>::set_exception_main(std::exception_ptr e) {
  this->_exception_thrown = e != nullptr;
}
template <typename Task>
inline void ContextManager<Task>::set_exception_io(usize idx, std::exception_ptr e) {
  assert(idx < this->_exception_ptrs_io.size());
  assert(this->_exception_ptrs_io[idx] == nullptr);
  this->_exception_thrown = e != nullptr;
  this->_exception_ptrs_io[idx] = e;
}
template <typename Task>
inline void ContextManager<Task>::set_exception_worker(usize idx, std::exception_ptr e) {
  assert(idx < this->_exception_ptrs_workers.size());
  assert(this->_exception_ptrs_workers[idx] == nullptr);
  this->_exception_thrown = e != nullptr;
  this->_exception_ptrs_workers[idx] = e;
}

template <typename Task>
inline bool ContextManager<Task>::exception_thrown() const noexcept {
  return this->_exception_thrown.load();
}

template <typename Task>
inline usize ContextManager<Task>::num_submitted() const noexcept {
  return this->_num_in.load();
}
template <typename Task>
inline usize ContextManager<Task>::num_completed() const noexcept {
  return this->_num_out.load();
}

template <typename Task>
template <typename ContextManager<Task>::Status s>
inline moodycamel::ProducerToken ContextManager<Task>::register_producer() {
  static_assert(s == Status::PENDING || s == Status::COMPLETED);
  if constexpr (s == Status::PENDING) {
    return moodycamel::ProducerToken(this->_pending);
  } else {
    return moodycamel::ProducerToken(this->_finished);
  }
}

template <typename Task>
template <typename ContextManager<Task>::Status s>
inline moodycamel::ConsumerToken ContextManager<Task>::register_consumer() {
  static_assert(s == Status::PENDING || s == Status::COMPLETED);
  if constexpr (s == Status::PENDING) {
    return moodycamel::ConsumerToken(this->_pending);
  } else {
    return moodycamel::ConsumerToken(this->_finished);
  }
}

template <typename Task>
inline void ContextManager<Task>::rethrow_exceptions() const {
  std::vector<std::string> messages;
  auto process_exceptions = [&](const auto& exceptions) {
    bool unhandled_except_thrown = false;
    for (const auto& except : exceptions) {
      if (except == nullptr) {
        continue;
      }
      try {
        std::rethrow_exception(except);
      } catch (const std::exception& e) {
        messages.emplace_back(e.what());
      } catch (...) {
        unhandled_except_thrown = true;
      }
    }
    if (unhandled_except_thrown) {
      messages.emplace_back(
          "An unhandled exception was caught! This should never happen! If you see "
          "this message, please file an issue on GitHub.");
    }
  };

  process_exceptions(this->_exception_ptrs_workers);
  process_exceptions(this->_exception_ptrs_io);

  throw std::runtime_error(fmt::format(
      FMT_STRING("The following error(s) occurred while simulating loop extrusion:\n - {}"),
      fmt::join(messages, "\n - ")));
}

template <typename Task>
inline void ContextManager<Task>::check_exceptions() {
  if (MODLE_UNLIKELY(this->exception_thrown())) {
    spdlog::error(FMT_STRING("MoDLE encountered an exception. Shutting down worker threads..."));
    this->shutdown();
    this->rethrow_exceptions();
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
  this->_worker_tpool.push_task(lambda);
}

template <typename Task>
template <typename TaskLambda>
inline void ContextManager<Task>::spawn_io_thread(TaskLambda lambda) {
  this->_io_tpool.push_task(lambda);
}

template <typename Task>
inline void ContextManager<Task>::shutdown() {
  this->_shutdown_requested = true;
  this->close_queue();
  this->_worker_tpool.pause();
  this->_io_tpool.pause();
  this->_worker_tpool.wait_for_tasks();
  this->_io_tpool.wait_for_tasks();
}

template <typename Task>
inline void ContextManager<Task>::init_model_state_logger(std::filesystem::path path,
                                                          std::string_view header) {
  assert(!this->_state_logger_ptr);
  assert(!std::filesystem::exists(path));
  spdlog::debug(FMT_STRING("[io] initializing model state logger at {}..."), path);
  compressed_io::Writer(path).write(header);
  this->_state_logger_ptr = std::make_unique<StateLoggerAggregator>(std::move(path));
}

template <typename Task>
inline void ContextManager<Task>::append_to_model_state_log(const std::filesystem::path& path,
                                                            bool remove_file_after_append) {
  assert(!!this->_state_logger_ptr);
  spdlog::debug(FMT_STRING("[io] appending {} to log file at {}..."), path,
                this->_state_logger_ptr->path());
  this->_state_logger_ptr->append(path, remove_file_after_append);
}

}  // namespace modle
