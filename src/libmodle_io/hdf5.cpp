#include "modle/hdf5.hpp"

// IWYU pragma: no_include <H5DaccProp.h>
// IWYU pragma: no_include <H5DataSet.h>
// IWYU pragma: no_include <H5DataSpace.h>
// IWYU pragma: no_include <H5DcreatProp.h>
// IWYU pragma: no_include <H5Exception.h>
// IWYU pragma: no_include <H5File.h>
// IWYU pragma: no_include <H5Group.h>
// IWYU pragma: no_include <H5PredType.h>
// IWYU pragma: no_include <H5Public.h>
// IWYU pragma: no_include "H5Public.h"
// IWYU pragma: no_include <H5Ppublic.h>
// IWYU pragma: no_include "H5Ppublic.h"
// IWYU pragma: no_include <H5SPublic.h>
// IWYU pragma: no_include <H5StrType.h>

#include <H5Cpp.h>                 // IWYU pragma: keep
#include <absl/strings/str_cat.h>  // for StrCat
#include <absl/strings/strip.h>    // for ConsumePrefix, StripPrefix, StripSuffix
#include <fcntl.h>                 // for SEEK_END, SEEK_SET
#include <fmt/format.h>            // for format, FMT_STRING

#include <algorithm>    // for max
#include <cassert>      // for assert
#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for int64_t
#include <cstdio>       // for fclose, fseek, tmpfile, ferror, fread, ftell, FILE
#include <memory>       // for unique_ptr
#include <stdexcept>    // for runtime_error
#include <string>       // for string, basic_string
#include <string_view>  // for string_view
#include <vector>       // for vector

namespace modle::hdf5 {

std::string construct_error_stack() {
  std::string buff;
  auto fp = std::unique_ptr<FILE, decltype(&fclose)>(std::tmpfile(), &fclose);
  if (fp) {  // TODO: Make this portable
    H5::Exception::printErrorStack(fp.get());
    fseek(fp.get(), 0L, SEEK_END);
    const auto buff_capacity = ftell(fp.get());
    if (buff_capacity < 0) {
      throw std::runtime_error(
          "h5pp::construct_error_stack(): unable to determine buffer size required to store an "
          "error message");
    }
    auto buff_size = static_cast<size_t>(buff_capacity);
    buff.resize(buff_size);
    if (fseek(fp.get(), 0L, SEEK_SET) != 0) {
      throw std::runtime_error(
          "h5pp::construct_error_stack(): failed to seek to the beginning to a temporary file");
    }
    /* Read the entire file into memory. */
    buff_size = fread(buff.data(), sizeof(char), buff_size, fp.get());
    buff.resize(buff_size);
    if (ferror(fp.get()) != 0) {
      throw std::runtime_error(
          "h5pp::construct_error_stack(): failed to read error message from temporary file");
    }
  } else {
    throw std::runtime_error("h5pp::construct_error_stack(): unable to create temporary file");
  }
  H5::Exception::clearErrorStack();
  return buff;
}

hsize_t read_str(const H5::DataSet &dataset, std::string &buff, hsize_t file_offset) {
  try {
    H5::Exception::dontPrint();

    constexpr hsize_t DIMS{1};
    constexpr hsize_t RANK{1};

    auto type = dataset.getDataType();
    auto file_space = dataset.getSpace();

    auto mem_space{H5::DataSpace(RANK, &DIMS)};

    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &RANK, &file_offset);
    dataset.read(buff, type, mem_space, file_space);

    return file_offset + 1;

  } catch ([[maybe_unused]] const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to read a string from dataset '{}' at offset {} using a "
                               "buffer size of 1:\n{}"),
                    dataset.getObjName(), file_offset, construct_error_stack()));
  }
}

hsize_t read_strings(const H5::DataSet &dataset, std::vector<std::string> &buff,
                     hsize_t file_offset) {
  try {
    H5::Exception::dontPrint();
    constexpr hsize_t DIMS = 1;
    constexpr hsize_t RANK = 1;
    auto type = dataset.getDataType();
    auto file_space = dataset.getSpace();

    const auto n_strings = static_cast<hsize_t>(file_space.getSimpleExtentNpoints());
    buff.resize(n_strings);

    auto mem_space{H5::DataSpace(RANK, &DIMS)};

    for (hsize_t i = 0U; i < n_strings; ++i) {
      file_space = dataset.getSpace();
      file_space.selectHyperslab(H5S_SELECT_SET, &RANK, &i);
      dataset.read(buff[i], type, mem_space, file_space);
    }

    return n_strings;

  } catch ([[maybe_unused]] const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to read strings from dataset '{}' at offset {} using a "
                               "buffer size of {}:\n{}"),
                    dataset.getObjName(), file_offset, buff.size(), construct_error_stack()));
  }
}

std::string read_str(const H5::DataSet &dataset, hsize_t file_offset) {
  std::string buff;
  (void)read_str(dataset, buff, file_offset);
  return buff;
}
std::vector<std::string> read_strings(const H5::DataSet &dataset, hsize_t file_offset) {
  std::vector<std::string> buff;
  (void)read_strings(dataset, buff, file_offset);
  return buff;
}

std::string read_attribute_str(H5::H5File &f, std::string_view attr_name, std::string_view path) {
  std::string buff;
  read_attribute(f, attr_name, buff, path);
  return buff;
}

std::string read_attribute_str(std::string_view path_to_file, std::string_view attr_name,
                               std::string_view path) {
  std::string buff;
  H5::H5File f({path_to_file.data(), path_to_file.size()}, H5F_ACC_RDONLY);
  read_attribute(f, attr_name, buff, path);
  return buff;
}

int64_t read_attribute_int(H5::H5File &f, std::string_view attr_name, std::string_view path) {
  int64_t buff;  // NOLINT
  read_attribute(f, attr_name, buff, path);
  return buff;
}

int64_t read_attribute_int(std::string_view path_to_file, std::string_view attr_name,
                           std::string_view path) {
  int64_t buff;  // NOLINT
  H5::H5File f({path_to_file.data(), path_to_file.size()}, H5F_ACC_RDONLY);
  read_attribute(f, attr_name, buff, path);
  return buff;
}

bool has_attribute(const H5::Group &g, std::string_view attr_name) {
  absl::ConsumePrefix(&attr_name, "/");
  return g.attrExists(std::string{attr_name});
}

bool has_attribute(const H5::DataSet &d, std::string_view attr_name) {
  absl::ConsumePrefix(&attr_name, "/");
  return d.attrExists(std::string{attr_name});
}

bool has_attribute(H5::H5File &f, std::string_view attr_name, std::string_view path) {
  auto g = f.openGroup(std::string{path});

  return has_attribute(g, attr_name);
}

bool has_group(H5::H5File &f, std::string_view name, std::string_view root_path) {
  H5O_info2_t info;
  const auto path =
      absl::StrCat(absl::StripSuffix(root_path, "/"), "/", absl::StripPrefix(name, "/"));
  assert(!path.empty());
  size_t pos = 0;
  do {
    pos = path.find_first_of('/', pos + 1);
    if (!f.nameExists(std::string{path.substr(0, pos)})) {
      return false;
    }
  } while (pos != std::string::npos);
  f.getObjinfo(path, info);
  if (info.type == H5O_TYPE_GROUP) {
    return true;
  }
  throw std::runtime_error(
      fmt::format(FMT_STRING("'{}' exists but is not a group"), absl::StrCat(root_path, name)));
}

bool has_dataset(H5::H5File &f, std::string_view name, std::string_view root_path) {
  H5O_info2_t info;
  const auto path =
      absl::StrCat(absl::StripSuffix(root_path, "/"), "/", absl::StripPrefix(name, "/"));
  assert(!path.empty());
  size_t pos = 0;
  do {
    pos = path.find_first_of('/', pos + 1);
    if (!f.nameExists(std::string{path.substr(0, pos)})) {
      return false;
    }
  } while (pos != std::string::npos);
  f.getObjinfo(path, info);

  return info.type == H5O_TYPE_DATASET;
}

H5::H5File open_file_for_reading(std::string_view path_to_file) {
  try {
    H5::Exception::dontPrint();
    H5::H5File f(std::string{path_to_file.data(), path_to_file.size()}, H5F_ACC_RDONLY);
    return f;
  } catch ([[maybe_unused]] const H5::Exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Failed to open file {} for reading:\n{}"),
                                         path_to_file, construct_error_stack()));
  }
}
}  // namespace modle::hdf5
