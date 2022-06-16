// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/strings/match.h>  // for EndsWith, StartsWith

#include <atomic>                           // for atomic
#include <boost/filesystem/operations.hpp>  // for unique_path
#include <filesystem>
#include <utility>  // for move

namespace modle::test {
// The point of this class is to provide a reliable way to create a directory that automatically
// deletes istelf and its content by default.
// This can be prevented by setting the internal flag to true.
// The default constructor will create a unique, randomly named directory under the system tmpdir
class SelfDeletingFolder {
  std::filesystem::path _path;
  std::atomic<bool> _delete_on_destruction{true};

 public:
  [[maybe_unused]] SelfDeletingFolder() {
    std::filesystem::path tmpdir{};
    try {
      tmpdir = std::filesystem::temp_directory_path();
    } catch (const std::filesystem::filesystem_error& e) {
      // Workaround spurious CI failures due to missing /tmp folder exception
      assert(absl::StartsWith(e.what(), "std::filesystem::temp_directory_path: Not a directory"));
      assert(absl::EndsWith(e.what(), "\"/tmp\""));
      tmpdir = "test/data/unit_tests/scratch";
    }

    _path = tmpdir / boost::filesystem::unique_path().string();
    std::filesystem::create_directories(_path);
  }

  [[maybe_unused]] explicit SelfDeletingFolder(std::filesystem::path path,
                                               bool delete_on_destruction = true)
      : _path(std::move(path)), _delete_on_destruction(delete_on_destruction) {
    std::filesystem::create_directories(_path);
  }

  [[maybe_unused]] explicit SelfDeletingFolder(bool delete_on_destruction) : SelfDeletingFolder() {
    this->set_delete_on_destruction(delete_on_destruction);
  }

  SelfDeletingFolder(const SelfDeletingFolder& other) = delete;
  SelfDeletingFolder(SelfDeletingFolder&& other) = delete;

  ~SelfDeletingFolder() {
    if (this->get_delete_on_destruction()) {
      std::filesystem::remove_all(this->_path);
    }
  }

  [[nodiscard]] const std::filesystem::path& operator()() const noexcept { return this->_path; }
  [[maybe_unused]] [[nodiscard]] bool get_delete_on_destruction() const noexcept {
    return this->_delete_on_destruction;
  }

  [[maybe_unused]] void set_delete_on_destruction(const bool flag) noexcept {
    this->_delete_on_destruction = flag;
  }

  SelfDeletingFolder& operator=(const SelfDeletingFolder& other) = delete;
  SelfDeletingFolder& operator=(SelfDeletingFolder&& other) = delete;
};
}  // namespace modle::test
