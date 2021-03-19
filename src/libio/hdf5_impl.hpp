#pragma once

// IWYU pragma: private, include "modle/hdf5.hpp"

#include <H5Cpp.h>                 // IWYU pragma: keep
#include <absl/strings/str_cat.h>  // for StrCat
#include <absl/strings/strip.h>    // for ConsumePrefix, StripPrefix
#include <fcntl.h>                 // for SEEK_END, SEEK_SET
#include <fmt/format.h>            // for format, FMT_STRING, to_string

#include <algorithm>    // for max
#include <cassert>      // for assert
#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for int64_t, int32_t, int16_t
#include <cstdio>       // for fclose, fseek, tmpfile, ferror
#include <limits>       // for numeric_limits
#include <memory>       // for unique_ptr
#include <stdexcept>    // for runtime_error, logic_error
#include <string>       // for string, basic_string
#include <string_view>  // for string_view
#include <type_traits>  // for decay_t, declval, remove_poi...
#include <variant>      // for visit, variant
#include <vector>       // for vector

#include "modle/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE...
#include "modle/utils.hpp"                       // for throw_with_trace, get_printable_ty...

namespace modle::hdf5 {

// Predeclare internal functions
using attr_types = std::variant<uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t,
                                int64_t, float, double, long double>;

attr_types getCpp_type(const H5::IntType &h5_type);
attr_types getCpp_type(const H5::FloatType &h5_type);

std::string construct_error_stack() {
  std::string buff;
  auto fp = std::unique_ptr<FILE, decltype(&fclose)>(std::tmpfile(), &fclose);
  if (fp) {  // TODO: Make this portable
    H5::Exception::printErrorStack(fp.get());
    fseek(fp.get(), 0L, SEEK_END);
    const auto bufsize = ftell(fp.get());
    if (bufsize < 0) {
      throw std::runtime_error(
          "h5pp::construct_error_stack(): unable to determine buffer size required to store an "
          "error message");
    }
    buff.resize(static_cast<std::size_t>(bufsize));
    if (fseek(fp.get(), 0L, SEEK_SET) != 0) {
      throw std::runtime_error(
          "h5pp::construct_error_stack(): failed to seek to the beginning to a temporary file");
    }
    /* Read the entire file into memory. */
    fread(buff.data(), sizeof(char), static_cast<std::size_t>(bufsize), fp.get());
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

template <typename S>
hsize_t write_str(const S &str, const H5::DataSet &dataset, const H5::StrType &str_type,
                  hsize_t file_offset) {
  static_assert(std::is_convertible_v<S, H5std_string>,
                "S should be a type convertible to std::string.");
  H5::Exception::dontPrint();
  constexpr hsize_t RANK{1};  // i.e. number of dimensions
  constexpr hsize_t BUFF_SIZE{1};
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset

#ifndef NDEBUG
  if (str.size() > str_type.getSize()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while writing strings '{}' to dataset "
                               "'{}' at offset {}:\n string does not fit in the receiving dataset: "
                               "string length: {}; Max. string length: {}"),
                    str, dataset.getObjName(), file_offset, str.size(), str_type.getSize()));
  }
#endif

  try {
    auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};
    // This is probably not very efficient, but it should be fine, given that we don't expect to
    // write that many strings
    hsize_t file_size{file_offset + 1};
    dataset.extend(&file_size);
    auto file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    std::string buff{str};
    buff.resize(str_type.getSize(), '\0');  // Pad strings with '\0'
    dataset.write(buff, str_type, mem_space, file_space);

    return file_size;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while writing strings '{}' to dataset "
                               "'{}' at offset {}:\n{}"),
                    str, dataset.getObjName(), file_offset, construct_error_stack()));
  }
}

template <typename CS>
hsize_t write_strings(const CS &strings, const H5::DataSet &dataset, const H5::StrType &str_type,
                      hsize_t file_offset, bool write_empty_strings) {
  static_assert(std::is_convertible_v<decltype(*strings.begin()), H5std_string> &&
                    std::is_convertible_v<decltype(*strings.end()), H5std_string>,
                "CS should be collection whose elements have a type convertible to std::string.");
  if (!write_empty_strings && strings.empty()) {
    return file_offset;
  }
  for (const auto &s : strings) {
    (void)write_str(s, dataset, str_type, file_offset++);
  }
  return file_offset;
}

template <typename N>
hsize_t write_number(N &num, const H5::DataSet &dataset, hsize_t file_offset) {
  static_assert(std::is_arithmetic_v<N>, "num should be a numeric type.");
  H5::Exception::dontPrint();

  constexpr hsize_t RANK{1};
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};
  constexpr hsize_t BUFF_SIZE{1};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    const auto *cbuff = reinterpret_cast<N(*)[BUFF_SIZE]>(&num);  // NOLINT
    DISABLE_WARNING_POP

    auto file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(cbuff, getH5_type<N>(), mem_space, file_space);

    return file_size;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while writing number {} to dataset "
                               "'{}' at offset {}:\n{}"),
                    num, dataset.getObjName(), file_offset, construct_error_stack()));
  }
}

template <typename CN>
hsize_t write_numbers(CN &numbers, const H5::DataSet &dataset, hsize_t file_offset) {
  static_assert(std::is_arithmetic_v<std::remove_pointer_t<decltype(std::declval<CN &>().data())>>,
                "numbers does not have a suitable ::data() member function.");
  static_assert(std::is_convertible_v<decltype(std::declval<CN &>().size()), std::size_t>,
                "numbers does not have a suitable ::size() member function.");
  if (numbers.empty()) {
    return file_offset;
  }
  H5::Exception::dontPrint();

  using N = std::remove_pointer_t<decltype(numbers.data())>;
  constexpr hsize_t RANK{1};
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};
  const auto num_type = getH5_type<N>();
  const hsize_t BUFF_SIZE{numbers.size()};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    const auto *cbuff = reinterpret_cast<N(*)[BUFF_SIZE]>(numbers.data());  // NOLINT
    DISABLE_WARNING_POP

    auto file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(cbuff, getH5_type<N>(), mem_space, file_space);

    return file_size;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error occurred while writing a collection of {} numbers to dataset "
            "'{}' at offset {}:\n{}"),
        numbers.size(), file_offset, dataset.getObjName(), construct_error_stack()));
  }
}

template <typename N>
hsize_t read_number(const H5::DataSet &dataset, N &buff, hsize_t file_offset) {
  static_assert(std::is_integral<N>::value, "I should be a numeric type.");
  H5::Exception::dontPrint();
  constexpr hsize_t DIMS = 1;
  constexpr hsize_t RANK = 1;
  try {
    auto mem_space{H5::DataSpace(RANK, &DIMS)};

    const auto type = dataset.getDataType();
    auto file_space = dataset.getSpace();
    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &RANK, &file_offset);
    dataset.read(&buff, type, mem_space, file_space);

    return file_offset + 1;

  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Failed to read a number from dataset '{}' "
                                                    "at offset {} using a buffer size of 1:\n{}"),
                                         dataset.getObjName(), file_offset,
                                         construct_error_stack()));
  }
}

template <typename CN>
hsize_t read_numbers(const H5::DataSet &dataset, CN &buff, hsize_t file_offset) {
  static_assert(std::is_arithmetic_v<std::decay_t<decltype(*std::declval<CN &>().begin())>>,
                "numbers does not have a suitable ::begin() member function.");
  static_assert(std::is_arithmetic_v<std::decay_t<decltype(*std::declval<CN &>().end())>>,
                "numbers does not have a suitable ::end() member function.");

  using N = std::remove_pointer_t<decltype(buff.data())>;
  H5::Exception::dontPrint();

  constexpr hsize_t RANK{1};  // i.e. number of dimensions
  hsize_t BUFF_SIZE{buff.size()};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE)};
    auto file_space = dataset.getSpace();
    if (const auto n_items = static_cast<hsize_t>(file_space.getSimpleExtentNpoints());
        n_items < BUFF_SIZE || BUFF_SIZE == 0) {
      BUFF_SIZE = n_items;
      buff.resize(BUFF_SIZE);
    }

    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    auto *cbuff = reinterpret_cast<N(*)[BUFF_SIZE]>(buff.data());  // NOLINT
    DISABLE_WARNING_POP

    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.read(cbuff, getH5_type<N>(), mem_space, file_space);

    return file_offset + BUFF_SIZE;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to read numbers from dataset '{}' at offset {} using a "
                               "buffer size of {} ({} bytes):\n{}"),
                    dataset.getObjName(), file_offset, BUFF_SIZE, BUFF_SIZE * sizeof(N),
                    construct_error_stack()));
  }
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

  } catch (const H5::Exception &e) {
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

  } catch (const H5::Exception &e) {
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

template <typename T>
void read_attribute(std::string_view path_to_file, std::string_view attr_name, T &buff,
                    std::string_view path) {
  H5::H5File f(std::string{path_to_file.data(), path_to_file.size()}, H5F_ACC_RDONLY);
  read_attribute(f, attr_name, buff, path);
}

template <typename T>
void read_attribute(H5::H5File &f, std::string_view attr_name, T &buff, std::string_view path) {
  auto g = f.openGroup(std::string{path});
  if (!has_attribute(g, attr_name)) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Unable to find an attribute named '{}' in group '/{}'"), attr_name, path));
  }
  auto attr = g.openAttribute(std::string{attr_name});

  // Handle simple cases, where HDF5 type and T are consistent
  if constexpr (std::is_constructible_v<H5std_string, T>) {
    if (const auto type = attr.getDataType(); type.getClass() == H5T_STRING) {
      buff.clear();
      attr.read(attr.getStrType(), buff);
      return;
    }
  }

  if constexpr (std::is_unsigned_v<T>) {
    if (const auto type = attr.getDataType(); type.getClass() == H5T_INTEGER &&
                                              attr.getIntType().getSign() == H5T_SGN_NONE &&
                                              attr.getIntType().getSize() == sizeof(T)) {
      attr.read(attr.getIntType(), &buff);
      return;
    }
  }

  if constexpr (std::is_integral_v<T>) {
    if (const auto type = attr.getDataType(); type.getClass() == H5T_INTEGER &&
                                              attr.getIntType().getSign() == H5T_SGN_2 &&
                                              attr.getIntType().getSize() == sizeof(T)) {
      attr.read(attr.getIntType(), &buff);
      return;
    }
  }

  if constexpr (std::is_floating_point_v<T>) {
    if (const auto type = attr.getDataType();
        type.getClass() == H5T_FLOAT && attr.getFloatType().getSize() == sizeof(T)) {
      attr.read(attr.getFloatType(), &buff);
      return;
    }
  }

  // Get HDF5 type info
  const auto type = attr.getDataType();
  const auto h5_class = type.getClass();

  if constexpr (std::is_arithmetic_v<T>) {
    if (type.getClass() == H5T_STRING) {
      std::string buff_;
      attr.read(attr.getStrType(), buff_);
      utils::parse_numeric_or_throw(buff_, buff);
      return;
    }
  }

  // Map HDF5 type to C++ type using std::variant
  assert(h5_class == H5T_INTEGER || h5_class == H5T_FLOAT);  // NOLINT
  auto vars =
      h5_class == H5T_INTEGER ? getCpp_type(attr.getIntType()) : getCpp_type(attr.getFloatType());
  // Visit the appropriate variant and read the attribute
  std::visit([&](auto &&var) { attr.read(attr.getDataType(), &var); }, vars);

  // Try to convert the variant to type T whenever it make sense
  if constexpr (std::is_constructible_v<H5std_string, T>) {
    std::visit([&](auto &&v) { buff = fmt::to_string(v); }, vars);
    return;
  }

  using VT = std::decay_t<decltype(vars)>;
  if constexpr (std::is_integral_v<T> != std::is_integral_v<VT> &&
                std::is_constructible_v<H5std_string, VT>) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Unable to store attribute '{}' of type {} into a buffer of type {}"), attr_name,
        utils::get_printable_type_name<T>(), utils::get_printable_type_name<VT>()));
  }
  std::visit(
      [&](auto &&v) {  // NOLINT unused-parameter
        if constexpr (std::is_arithmetic_v<T>) {
          DISABLE_WARNING_PUSH
          DISABLE_WARNING_IMPL_INT_TO_FLOAT
          DISABLE_WARNING_SIGN_COMPARE
          if (v >= std::numeric_limits<T>::min() && v <= std::numeric_limits<T>::max()) {  // NOLINT
            DISABLE_WARNING_POP
            buff = static_cast<T>(v);
          } else {
            throw std::runtime_error(fmt::format(
                FMT_STRING("Unable to store attribute '{}' with value {} in a buffer with "
                           "type {}. Reason: value does not fit in the receiver"),
                attr_name, v, utils::get_printable_type_name<T>()));  // NOLINT
          }
        }
      },
      vars);
}

template <typename T>
void read_attribute(const H5::DataSet &d, std::string_view attr_name, T &buff) {
  if (!has_attribute(d, attr_name)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find an attribute named '{}' in dataset '{}'"), attr_name,
                    d.getObjName()));
  }

  auto attr = d.openAttribute(std::string{attr_name});
  if constexpr (std::is_constructible_v<H5std_string, T>) {  // string-like (scalar)
    buff.clear();
    attr.read(attr.getStrType(), buff);
  } else if constexpr (std::is_arithmetic_v<T>) {  // numeric (scalar)
    attr.read(attr.getDataType(), &buff);
  } else {  // array of numbers (non-scalar)
    using N = std::remove_pointer_t<decltype(buff.data())>;
    if constexpr (!std::is_arithmetic_v<N>) {
      throw std::logic_error(
          "buff should be a container that is contiguous in memory, and whose elements have a "
          "numeric type.");
    }
    hsize_t buff_size{0};  // Figure out the appropriate buffer size
    attr.getSpace().getSimpleExtentDims(&buff_size);
    buff.resize(buff_size);
    assert(buff_size > 1);
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    auto *cbuff = reinterpret_cast<N(*)[buff_size]>(buff.data());  // NOLINT
    DISABLE_WARNING_POP
    const auto dtype = attr.getArrayType();  // It is important to get the array type, and not
                                             // just the data type here
    attr.read(attr.getArrayType(), cbuff);
  }
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

template <typename T>
inline void write_or_create_attribute(H5::H5File &f, std::string_view attr_name, T &buff,
                                      std::string_view path) {
  // TODO handle array attributes
  H5::DataSpace attr_space(H5S_SCALAR);
  if (!hdf5::has_group(f, path)) {
    throw std::runtime_error(fmt::format(FMT_STRING("Unable to find group '{}'"), path));
  }

  auto g = f.openGroup(std::string{path});
  if constexpr (std::is_convertible_v<T, H5std_string>) {
    if (hdf5::has_attribute(f, std::string{attr_name}, std::string{path})) {
      auto attr = g.openAttribute(std::string{attr_name});
      attr.write(attr.getDataType(), buff);
    } else {
      auto attr = g.createAttribute(std::string{attr_name}, METADATA_STR_TYPE(), attr_space);
      attr.write(attr.getDataType(), buff);
    }
  } else {
    if (hdf5::has_attribute(f, std::string{attr_name}, std::string{path})) {
      auto attr = g.openAttribute(std::string{attr_name});
      attr.write(attr.getDataType(), &buff);
    } else {
      auto attr = g.createAttribute(std::string{attr_name}, getH5_type<T>(), attr_space);
      attr.write(attr.getDataType(), &buff);
    }
  }
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
  std::size_t pos = 0;
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
  std::size_t pos = 0;
  do {
    pos = path.find_first_of('/', pos + 1);
    if (!f.nameExists(std::string{path.substr(0, pos)})) {
      return false;
    }
  } while (pos != std::string::npos);
  f.getObjinfo(path, info);

  return info.type == H5O_TYPE_DATASET;
}

template <typename T>
bool check_dataset_type(const H5::DataSet &dataset, T type, bool throw_on_failure) {
  static_assert(std::is_base_of_v<H5::DataType, T>,
                "type should have a type that is derived from H5::DataType (such as PredType, "
                "IntType or StrType).");
  const auto actual_type = dataset.getTypeClass();
  const auto expected_type = type.getClass();
  if (actual_type != expected_type) {
    if (throw_on_failure) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("'{}' exists but has incorrect datatype"), dataset.getObjName()));
    }
    return false;
  }
  const auto actual_size = dataset.getDataType().getSize();
  const auto expected_size = type.getSize();
  if (actual_size != expected_size) {
    if (throw_on_failure) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("'{}' exists but has incorrect datasize (expected {} bytes, got {})"),
          dataset.getObjName(), expected_size, actual_size));
    }
    return false;
  }

  if constexpr (std::is_same_v<T, H5::StrType>) {
    const auto expected_charset = type.getCset();
    const auto actual_charset = dataset.getStrType().getCset();
    if (actual_charset != expected_charset) {
      if (throw_on_failure) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("'{}' exists but has incorrect CharSet"), dataset.getObjName()));
      }
      return false;
    }
    const auto expected_padding = type.getStrpad();
    const auto actual_padding = dataset.getStrType().getStrpad();
    if (actual_padding != expected_padding) {
      if (throw_on_failure) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("'{}' exists but has incorrect String padding"), dataset.getObjName()));
      }
      return false;
    }
  } else {
    if (actual_type == H5T_INTEGER) {
      const auto actual_sign = dataset.getIntType().getSign();
      if (auto t = H5::IntType(type); t.getSign() != actual_sign) {
        if (throw_on_failure) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("'{}' exists but has incorrect signedness"), dataset.getObjName()));
        }
        return false;
      }
    }

    if (actual_type == H5T_FLOAT) {
      const auto actual_precision = dataset.getFloatType().getPrecision();
      if (auto t = H5::FloatType(type); t.getPrecision() != actual_precision) {
        if (throw_on_failure) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("'{}' exists but has incorrect precision"), dataset.getObjName()));
        }
        return false;
      }
    }
  }
  return true;
}

H5::H5File open_file_for_reading(std::string_view path_to_file) {
  try {
    H5::Exception::dontPrint();
    H5::H5File f(std::string{path_to_file.data(), path_to_file.size()}, H5F_ACC_RDONLY);
    return f;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Failed to open file {} for reading:\n{}"),
                                         path_to_file, construct_error_stack()));
  }
}

// Internal functions

template <typename DataType>
H5::PredType getH5_type() {
  // Taken from github.com/DavidAce/h5pp:
  // https://github.com/DavidAce/h5pp/blob/master/h5pp/include/h5pp/details/h5ppUtils.h
  // clang-format off
  using DecayType = typename std::decay_t<DataType>;
  if constexpr (std::is_same_v<DecayType, short>) return H5::PredType::NATIVE_SHORT;               // NOLINT
  if constexpr (std::is_same_v<DecayType, int>) return H5::PredType::NATIVE_INT;                   // NOLINT
  if constexpr (std::is_same_v<DecayType, long>) return H5::PredType::NATIVE_LONG;                 // NOLINT
  if constexpr (std::is_same_v<DecayType, long long>) return H5::PredType::NATIVE_LLONG;           // NOLINT
  if constexpr (std::is_same_v<DecayType, unsigned short>) return H5::PredType::NATIVE_USHORT;     // NOLINT
  if constexpr (std::is_same_v<DecayType, unsigned int>) return H5::PredType::NATIVE_UINT;         // NOLINT
  if constexpr (std::is_same_v<DecayType, unsigned long>) return H5::PredType::NATIVE_ULONG;       // NOLINT
  if constexpr (std::is_same_v<DecayType, unsigned long long>) return H5::PredType::NATIVE_ULLONG; // NOLINT
  if constexpr (std::is_same_v<DecayType, double>) return H5::PredType::NATIVE_DOUBLE;             // NOLINT
  if constexpr (std::is_same_v<DecayType, long double>) return H5::PredType::NATIVE_LDOUBLE;       // NOLINT
  if constexpr (std::is_same_v<DecayType, float>) return H5::PredType::NATIVE_FLOAT;               // NOLINT
  if constexpr (std::is_same_v<DecayType, bool>) return H5::PredType::NATIVE_UINT8;                // NOLINT
  if constexpr (std::is_same_v<DecayType, std::string>) return H5::PredType::C_S1;                 // NOLINT
  if constexpr (std::is_same_v<DecayType, char>) return H5::PredType::C_S1;                        // NOLINT
  // clang-format on

  throw std::logic_error("getH5_type(): Unable to map C++ type to a H5T");
}

H5::StrType METADATA_STR_TYPE() {
  H5::StrType t(H5::PredType::C_S1, H5T_VARIABLE);
  t.setStrpad(H5T_STR_NULLTERM);
  t.setCset(H5T_CSET_UTF8);
  return t;
}

attr_types getCpp_type(const H5::IntType &h5_type) {
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  attr_types type;
  assert(h5_type.getSign() != H5T_SGN_ERROR);  // NOLINT
  const bool is_signed = h5_type.getSign() == H5T_SGN_NONE;
  const auto size = h5_type.getSize();
  assert(size >= 1 && size <= 8);  // NOLINT
  if (is_signed) {
    switch (size) {
      case (sizeof(int8_t)):
        type = static_cast<int8_t>(0);
        return type;
      case (sizeof(int16_t)):
        type = static_cast<int16_t>(0);
        return type;
      case (sizeof(int32_t)):
        type = static_cast<int32_t>(0);
        return type;
      case (sizeof(int64_t)):
        type = static_cast<int64_t>(0);
        return type;
      default:
        utils::throw_with_trace(std::logic_error("Unreachable code"));
    }
  } else {
    switch (size) {
      case (sizeof(uint8_t)):
        type = static_cast<uint8_t>(0);
        return type;
      case (sizeof(uint16_t)):
        type = static_cast<uint16_t>(0);
        return type;
      case (sizeof(uint32_t)):
        type = static_cast<uint32_t>(0);
        return type;
      case (sizeof(uint64_t)):
        type = static_cast<uint64_t>(0);
        return type;
      default:
        utils::throw_with_trace(std::logic_error("Unreachable code"));
    }
  }
  DISABLE_WARNING_POP
  utils::throw_with_trace(std::logic_error("Unreachable code"));
}

attr_types getCpp_type(const H5::FloatType &h5_type) {
  attr_types type;
  const auto size = h5_type.getSize();
  switch (size) {
    case (sizeof(float)):
      type = static_cast<float>(0);
      return type;
    case (sizeof(double)):
      type = static_cast<double>(0);
      return type;
    case (sizeof(long double)):
      type = static_cast<long double>(0);
      return type;
    default:
      utils::throw_with_trace(std::logic_error("Unreachable code"));
  }
  utils::throw_with_trace(std::logic_error("Unreachable code"));
}

}  // namespace modle::hdf5
