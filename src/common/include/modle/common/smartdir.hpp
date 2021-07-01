#pragma once

#include <atomic>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>

namespace modle {
// The point of this class is to provide a reliable way to create a directory that automatically
// deletes istelf as well as its content by default.
// This can be prevented by setting the internal flag to true.
// The default constructor will create a unique, randomly named directory under the system tmpdir
class SmartDir {  // NOLINT
 public:
  [[maybe_unused]] inline SmartDir() {
    boost::filesystem::path path{
        (boost::filesystem::temp_directory_path() / boost::filesystem::unique_path()).string()};
    boost::filesystem::create_directories(path);
    _path = path;
  }
  [[maybe_unused]] explicit inline SmartDir(boost::filesystem::path path,
                                            bool delete_on_destruction = true)
      : _path(std::move(path)), _delete_on_destruction(delete_on_destruction) {
    boost::filesystem::create_directories(_path);
  }

  [[maybe_unused]] explicit inline SmartDir(bool delete_on_destruction) : SmartDir() {
    this->set_delete_on_destruction(delete_on_destruction);
  }

  inline ~SmartDir() {
    if (this->get_delete_on_destruction()) {
      boost::filesystem::remove_all(this->_path);
    }
  }

  [[nodiscard]] inline const boost::filesystem::path& operator()() const noexcept {
    return this->_path;
  }
  [[maybe_unused]] [[nodiscard]] inline bool get_delete_on_destruction() const noexcept {
    return this->_delete_on_destruction;
  }

  [[maybe_unused]] inline void set_delete_on_destruction(const bool flag) noexcept {
    this->_delete_on_destruction = flag;
  }

 private:
  boost::filesystem::path _path;
  std::atomic<bool> _delete_on_destruction{true};
};
}  // namespace modle
