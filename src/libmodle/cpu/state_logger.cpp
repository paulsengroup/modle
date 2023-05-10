// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/state_logger.hpp"

#include <fmt/format.h>
#include <fmt/std.h>

#include <cassert>
#include <exception>
#include <fstream>
#include <ios>
#include <stdexcept>
#include <string>

namespace modle {

StateLoggerAggregator::StateLoggerAggregator(std::filesystem::path path)
    : _path(std::move(path)),
      _fs(_path,
          std::ios_base::binary | std::ios_base::in | std::ios_base::out | std::ios_base::ate) {
  _fs.exceptions(_fs.exceptions() | std::ios_base::failbit | std::ios_base::badbit);
  if (!_path.empty() && !_fs) {
    throw fmt::system_error(errno, FMT_STRING("failed to open file {} for writing"), _path);
  }
}

void StateLoggerAggregator::append(const std::filesystem::path &path,
                                   bool remove_source_after_append) {
  try {
    assert(std::filesystem::exists(path));
    {
      auto input = std::ifstream(path, std::ios_base::binary);
      const auto lck = std::scoped_lock(this->_mtx);
      this->_fs << input.rdbuf();
    }
    if (remove_source_after_append) {
      std::filesystem::remove(path);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("the following error occurred while appending file {} to file {}: {}"), path,
        this->_path, std::string{e.what()}));
  }
}

const std::filesystem::path &StateLoggerAggregator::path() const noexcept { return this->_path; }

}  // namespace modle
