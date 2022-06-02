// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <atomic>                           // for atomic
#include <boost/filesystem/operations.hpp>  // for create_directories, remove_all, temp_director...
#include <boost/filesystem/path.hpp>        // for path, operator/
#include <utility>                          // for move

namespace modle::test {
// The point of this class is to provide a reliable way to create a directory that automatically
// deletes istelf as well as its content by default.
// This can be prevented by setting the internal flag to true.
// The default constructor will create a unique, randomly named directory under the system tmpdir
class SelfDeletingFolder {
  boost::filesystem::path _path;
  std::atomic<bool> _delete_on_destruction{true};

 public:
  [[maybe_unused]] SelfDeletingFolder()
      : _path(boost::filesystem::temp_directory_path() / boost::filesystem::unique_path()) {
    boost::filesystem::create_directories(_path);
  }
  [[maybe_unused]] explicit SelfDeletingFolder(boost::filesystem::path path,
                                               bool delete_on_destruction = true)
      : _path(std::move(path)), _delete_on_destruction(delete_on_destruction) {
    boost::filesystem::create_directories(_path);
  }

  [[maybe_unused]] explicit SelfDeletingFolder(bool delete_on_destruction) : SelfDeletingFolder() {
    this->set_delete_on_destruction(delete_on_destruction);
  }

  SelfDeletingFolder(const SelfDeletingFolder& other) = delete;
  SelfDeletingFolder(SelfDeletingFolder&& other) = delete;

  ~SelfDeletingFolder() {
    if (this->get_delete_on_destruction()) {
      boost::filesystem::remove_all(this->_path);
    }
  }

  [[nodiscard]] const boost::filesystem::path& operator()() const noexcept { return this->_path; }
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
