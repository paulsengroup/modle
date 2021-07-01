#pragma once

#ifndef MODLE_WITH_BOOST_FILESYSTEM
#include <boost/filesystem.hpp>
#include <boost/filesystem/file_status.hpp>  // for file_type (is_fifo replacement)
#include <boost/filesystem/operations.hpp>   // for file_size
#include <string>                            // for string
#include <string_view>                       // for string_view
#include <system_error>                      // for system_error
#else
#include <filesystem>
#endif

namespace modle {
#ifndef MODLE_WITH_BOOST_FILESYSTEM
namespace filesystem = boost::filesystem;

[[nodiscard]] inline bool path_is_fifo(const boost::filesystem::path& p) {
  return boost::filesystem::status(p).type() == boost::filesystem::fifo_file;
}

[[nodiscard]] inline size_t file_size(const boost::filesystem::path& p) { return p.size(); }
[[nodiscard]] inline size_t file_size(std::string_view s) {
  return boost::filesystem::file_size(std::string{s});
}

inline void remove_path_recursive(const boost::filesystem::path& p) noexcept {
  std::system_error c;
}

#else
namespace filesystem = std::filesystem;
[[nodiscard]] inline bool path_is_fifo(const boost::filesystem::path& p) {
  return boost::filesystem::is_fifo(p);
}

[[nodiscard]] inline size_t file_size(const boost::filesystem::path& p) {
  return boost::filesystem::file_size(p);
}
#endif
}  // namespace modle
