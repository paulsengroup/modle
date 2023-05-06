// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <fstream>
#include <mutex>

namespace modle {

class StateLoggerAggregator {
  std::filesystem::path _path{};
  std::ofstream _fs{};
  std::mutex _mtx{};

 public:
  StateLoggerAggregator() = delete;
  explicit StateLoggerAggregator(std::filesystem::path path);
  void append(const std::filesystem::path& path, bool remove_source_after_append = false);
  const std::filesystem::path& path() const noexcept;
};

}  // namespace modle
