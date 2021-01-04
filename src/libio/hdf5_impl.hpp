#pragma once

#include <H5Cpp.h>
#include <absl/time/clock.h>

#include <algorithm>
#include <experimental/type_traits>
#include <string_view>
#include <type_traits>
#include <vector>

#include "modle/contacts.hpp"
#include "modle/suppress_compiler_warnings.hpp"
#include "modle/utils.hpp"

namespace modle::hdf5 {

template <typename DataType>
H5::PredType getH5_type() {
  // Taken from github.com/DavidAce/h5pp:
  // https://github.com/DavidAce/h5pp/blob/master/h5pp/include/h5pp/details/h5ppUtils.h
  // clang-format off
  using DecayType = typename std::decay<DataType>::type;
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

  throw std::logic_error("getH5_type(): Unable to map C++ type to a H5T.");
}

std::string construct_error_stack() {
  std::string buff;
  // The '+' in front of the declaration is a trick that allows to use the lambda as C function
  // pointer
  auto err_handler = +[](int n, const H5E_error_t *err_desc, void *client_data) {
    auto *buff_ = reinterpret_cast<std::string *>(client_data);  // NOLINT
    absl::StrAppendFormat(
        buff_,
        "%d:\n  Major: %d\n  Minor: %d\n  Function %s raised an error while processing file "
        "'%s':\n    %s\n",
        n, err_desc->maj_num, err_desc->min_num, err_desc->func_name, err_desc->file_name,
        err_desc->desc != nullptr ? err_desc->desc : "error description not available");
    return 0;
  };

  H5::Exception::walkErrorStack(H5E_WALK_DOWNWARD,
                                reinterpret_cast<H5E_walk2_t>(err_handler),  // NOLINT
                                buff.data());
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
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error occurred while writing a collection of strings '{}' to dataset "
            "'{}' at offset {}:\n{}"),
        str, dataset.getObjName(), file_offset, construct_error_stack()));
  }
}

template <typename CS>
hsize_t write_strings(const CS &strings, const H5::DataSet &dataset, const H5::StrType &str_type,
                      hsize_t file_offset) {
  static_assert(std::is_convertible_v<decltype(*strings.begin()), H5std_string> &&
                    std::is_convertible_v<decltype(*strings.end()), H5std_string>,
                "CS should be collection whose elements have a type convertible to std::string.");
  assert(!strings.empty());
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
  H5::Exception::dontPrint();

  using N = std::remove_pointer_t<decltype(numbers.data())>;
  assert(!numbers.empty());
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
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "Failed to read a number from dataset '{}' at offset {} using a buffer size of 1:\n{}"),
        dataset.getObjName(), file_offset, construct_error_stack()));
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

bool has_attribute(const H5::Group &g, std::string_view attr_name) {
  absl::ConsumePrefix(&attr_name, "/");
  return g.attrExists(std::string{attr_name});
}

bool has_attribute(H5::H5File &f, std::string_view attr_name, std::string_view path) {
  auto g = f.openGroup(std::string{path});

  return has_attribute(g, attr_name);
}

template <typename T>
void read_attribute(H5::H5File &f, std::string_view attr_name, T &buff, std::string_view path) {
  auto g = f.openGroup(std::string{path});
  if (!has_attribute(g, attr_name)) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Unable to find an attribute named '{}' in path '/{}'"), attr_name, path));
  }

  auto attr = g.openAttribute(std::string{attr_name});
  if constexpr (std::is_constructible_v<H5std_string, T>) {
    buff.clear();
    attr.read(attr.getStrType(), buff);
  } else {
    attr.read(attr.getDataType(), &buff);
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

bool group_exists(H5::H5File &f, std::string_view name, std::string_view root_path) {
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

bool dataset_exists(H5::H5File &f, std::string_view name, std::string_view root_path) {
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

}  // namespace modle::hdf5
