// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "modle/hdf5.hpp"

// IWYU pragma: no_include <boost/core/checked_delete.hpp>
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <H5ArrayType.h>
// IWYU pragma: no_include <H5Attribute.h>
// IWYU pragma: no_include <H5DaccProp.h>
// IWYU pragma: no_include <H5DataSet.h>
// IWYU pragma: no_include <H5DataSpace.h>
// IWYU pragma: no_include <H5DataType.h>
// IWYU pragma: no_include <H5DcreatProp.h>
// IWYU pragma: no_include <H5Exception.h>
// IWYU pragma: no_include <H5File.h>
// IWYU pragma: no_include <H5FloatType.h>
// IWYU pragma: no_include <H5Group.h>
// IWYU pragma: no_include <H5IntType.h>
// IWYU pragma: no_include <H5PredType.h>
// IWYU pragma: no_include <H5Public.h>
// IWYU pragma: no_include "H5Public.h"
// IWYU pragma: no_include <H5Ppublic.h>
// IWYU pragma: no_include "H5Ppublic.h"
// IWYU pragma: no_include <H5public.h>
// IWYU pragma: no_include <H5StrType.h>
// IWYU pragma: no_include <H5Spublic.h>
// IWYU pragma: no_include <H5public.h>

#include <H5Cpp.h>       // IWYU pragma: keep
#include <fmt/format.h>  // for FMT_STRING, format, to_string

#include <cassert>      // for assert
#include <limits>       // for numeric_limits
#include <stdexcept>    // for runtime_error, logic_error
#include <string>       // for string
#include <string_view>  // for string_view
#include <type_traits>  // for decay_t, declval, remove_pointer_t
#include <vector>       // for vector

#include "modle/common/common.hpp"                      // for i32, i16, i64, i8
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNI...
#include "modle/common/utils.hpp"                       // for throw_with_trace, get_printable_ty...

namespace modle::hdf5 {

template <class S>
hsize_t write_str(const S &str_, const H5::DataSet &dataset, const H5::StrType &str_type,
                  hsize_t file_offset) {
  static_assert(std::is_constructible_v<S, H5std_string>,
                "S should be a type that can be used to construct a std::string.");
  const auto lck = internal::lock();
  H5::Exception::dontPrint();
  const hsize_t RANK{1};  // i.e. number of dimensions
  const hsize_t BUFF_SIZE{1};
  const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset

  const auto &str = [&]() {  // This is mostly useful to allow std::string_view to be passed in
    if constexpr (std::is_convertible_v<S, H5std_string>) {
      return str_;
    }
    return std::string{str_};
  }();

  if constexpr (utils::ndebug_not_defined()) {
    if (str.size() > str_type.getSize()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("The following error occurred while writing strings '{}' to dataset "
                     "'{}' at offset {}:\n string does not fit in the receiving dataset: "
                     "string length: {}; Max. string length: {}"),
          str, dataset.getObjName(), file_offset, str.size(), str_type.getSize()));
    }
  }

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
                               "'{}' at offset {}: {}"),
                    str, dataset.getObjName(), file_offset, construct_error_stack(e)));
  }
}

template <class CS>
hsize_t write_strings(const CS &strings, const H5::DataSet &dataset, const H5::StrType &str_type,
                      hsize_t file_offset, bool write_empty_strings) {
  static_assert(std::is_convertible_v<decltype(*strings.begin()), H5std_string> &&
                    std::is_convertible_v<decltype(*strings.end()), H5std_string>,
                "CS should be collection whose elements have a type convertible to std::string.");
  if (!write_empty_strings && strings.empty()) {
    return file_offset;
  }
  const auto lck = internal::lock();
  for (const auto &s : strings) {
    std::ignore = write_str(s, dataset, str_type, file_offset++);
  }
  return file_offset;
}

template <class N>
hsize_t write_number(N &num, const H5::DataSet &dataset, hsize_t file_offset) {
  static_assert(std::is_arithmetic_v<N>, "num should be a numeric type.");
  const auto lck = internal::lock();
  H5::Exception::dontPrint();

  const hsize_t RANK{1};
  const hsize_t MAXDIMS{H5S_UNLIMITED};
  const hsize_t BUFF_SIZE{1};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);

    auto file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(&num, getH5_type<N>(), mem_space, file_space);

    return file_size;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while writing number {} to dataset "
                               "'{}' at offset {}: {}"),
                    num, dataset.getObjName(), file_offset, construct_error_stack(e)));
  }
}

template <class CN>
hsize_t write_numbers(CN &numbers, const H5::DataSet &dataset, hsize_t file_offset) {
  static_assert(std::is_arithmetic_v<std::remove_pointer_t<decltype(std::declval<CN &>().data())>>,
                "numbers does not have a suitable ::data() member function.");
  static_assert(std::is_convertible_v<decltype(std::declval<CN &>().size()), usize>,
                "numbers does not have a suitable ::size() member function.");
  if (numbers.empty()) {
    return file_offset;
  }
  const auto lck = internal::lock();
  H5::Exception::dontPrint();

  using N = std::remove_pointer_t<decltype(numbers.data())>;
  const hsize_t RANK{1};
  const hsize_t MAXDIMS{H5S_UNLIMITED};
  const hsize_t BUFF_SIZE{numbers.size()};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);

    auto file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(numbers.data(), getH5_type<N>(), mem_space, file_space);

    return file_size;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error occurred while writing a collection of {} numbers to dataset "
            "'{}' at offset {}: {}"),
        numbers.size(), file_offset, dataset.getObjName(), construct_error_stack(e)));
  }
}

template <class N>
hsize_t read_number(const H5::DataSet &dataset, N &buff, hsize_t file_offset) {
  static_assert(std::is_integral<N>::value, "I should be a numeric type.");
  const auto lck = internal::lock();
  H5::Exception::dontPrint();
  const hsize_t DIMS = 1;
  const hsize_t RANK = 1;
  try {
    auto mem_space{H5::DataSpace(RANK, &DIMS)};

    const auto type = dataset.getDataType();
    auto file_space = dataset.getSpace();
    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &RANK, &file_offset);
    dataset.read(&buff, type, mem_space, file_space);

    return file_offset + 1;

  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Failed to read a number from dataset \"{}\" "
                                                    "at offset {} using a buffer size of 1: {}"),
                                         dataset.getObjName(), file_offset,
                                         construct_error_stack(e)));
  }
}

template <class CN>
hsize_t read_numbers(const H5::DataSet &dataset, CN &buff, hsize_t file_offset) {
  static_assert(std::is_arithmetic_v<std::decay_t<decltype(*std::declval<CN &>().begin())>>,
                "numbers does not have a suitable ::begin() member function.");
  static_assert(std::is_arithmetic_v<std::decay_t<decltype(*std::declval<CN &>().end())>>,
                "numbers does not have a suitable ::end() member function.");

  using N = std::remove_pointer_t<decltype(buff.data())>;
  const auto lck = internal::lock();
  H5::Exception::dontPrint();

  const hsize_t RANK{1};  // i.e. number of dimensions
  hsize_t BUFF_SIZE{buff.size()};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE)};
    auto file_space = dataset.getSpace();
    if (const auto n_items = static_cast<hsize_t>(file_space.getSimpleExtentNpoints());
        n_items < BUFF_SIZE || BUFF_SIZE == 0) {
      BUFF_SIZE = n_items;
      assert(BUFF_SIZE > 0);
      buff.resize(BUFF_SIZE);
    }

    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.read(buff.data(), getH5_type<N>(), mem_space, file_space);

    return file_offset + BUFF_SIZE;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to read numbers from dataset \"{}\" at offset {} using a "
                               "buffer of size {} ({} bytes): {}"),
                    dataset.getObjName(), file_offset, BUFF_SIZE, BUFF_SIZE * sizeof(N),
                    construct_error_stack(e)));
  }
}

template <class T>
void read_attribute(const boost::filesystem::path &path_to_file, std::string_view attr_name,
                    T &buff, std::string_view path) {
  const auto lck = internal::lock();
  auto f = open_file_for_reading(path_to_file);
  read_attribute(f, attr_name, buff, path);
}

template <class T>
void read_attribute(H5::H5File &f, std::string_view attr_name, T &buff, std::string_view path) {
  const auto lck = internal::lock();
  auto g = hdf5::open_group(f, std::string{path});
  if (!hdf5::has_attribute(g, attr_name)) {
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

  // Map HDF5 type to C++ type using absl::variant
  assert(h5_class == H5T_INTEGER || h5_class == H5T_FLOAT);
  auto vars = h5_class == H5T_INTEGER ? get_cpp_arithmetic_type(attr.getIntType())
                                      : get_cpp_arithmetic_type(attr.getFloatType());
  // Visit the appropriate variant and read the attribute
  absl::visit([&](auto &&var) { attr.read(attr.getDataType(), &var); }, vars);

  // Try to convert the variant to type T whenever it make sense
  if constexpr (std::is_constructible_v<H5std_string, T>) {
    absl::visit([&](auto &&v) { buff = fmt::to_string(v); }, vars);
    return;
  }

  using VT = std::decay_t<decltype(vars)>;
  if constexpr (std::is_integral_v<T> != std::is_integral_v<VT> &&
                std::is_constructible_v<H5std_string, VT>) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Unable to store attribute '{}' of type {} into a buffer of type {}"), attr_name,
        utils::get_printable_type_name<T>(), utils::get_printable_type_name<VT>()));
  }
  absl::visit(
      [&](auto &&v) {  // NOLINT unused-parameter
        if constexpr (std::is_arithmetic_v<T>) {
          DISABLE_WARNING_PUSH
          DISABLE_WARNING_IMPL_INT_TO_FLOAT
          DISABLE_WARNING_SIGN_COMPARE
#if __GNUC__ < 8
          // The following two warnings seem to be a false positive produced by GCC 7.5
          DISABLE_WARNING_CONVERSION
          DISABLE_WARNING_SIGN_CONVERSION
#endif
          if (v >= (std::numeric_limits<T>::min)() && v <= (std::numeric_limits<T>::max)()) {
            DISABLE_WARNING_POP
            buff = static_cast<T>(v);
          } else {
            throw std::runtime_error(fmt::format(
                FMT_STRING("Unable to store attribute '{}' with value {} in a buffer with "
                           "type {}. Reason: value does not fit in the receiver"),
                attr_name, v, utils::get_printable_type_name<T>()));
          }
        }
      },
      vars);
}

template <class T>
void read_attribute(const H5::DataSet &d, std::string_view attr_name, T &buff) {
  if (!has_attribute(d, attr_name)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find an attribute named '{}' in dataset '{}'"), attr_name,
                    d.getObjName()));
  }

  const auto lck = internal::lock();
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
    assert(buff_size > 1);
    buff.resize(buff_size);
    attr.read(attr.getArrayType(), buff.data());
  }
}

template <class T>
T read_attribute(const boost::filesystem::path &path_to_file, std::string_view attr_name,
                 std::string_view path) {
  T buff;
  read_attribute(path_to_file, attr_name, buff, path);
  return buff;
}

template <class T>
T read_attribute(H5::H5File &f, std::string_view attr_name, std::string_view path) {
  T buff;
  read_attribute(f, attr_name, buff, path);
  return buff;
}

template <class T>
T read_attribute(const H5::DataSet &d, std::string_view attr_name) {
  T buff;
  read_attribute(d, attr_name, buff);
  return buff;
}

template <class T>
inline void write_or_create_attribute(H5::H5File &f, std::string_view attr_name, T &buff,
                                      std::string_view path) {
  // TODO handle array attributes
  const auto lck = internal::lock();
  H5::DataSpace attr_space(H5S_SCALAR);
  if (!hdf5::has_group(f, path)) {
    throw std::runtime_error(fmt::format(FMT_STRING("Unable to find group '{}'"), path));
  }

  const auto attribute_exists = hdf5::has_attribute(f, std::string{attr_name}, std::string{path});

  auto g = f.openGroup(std::string{path});
  if constexpr (std::is_convertible_v<T, H5std_string>) {
    if (attribute_exists) {
      auto attr = g.openAttribute(std::string{attr_name});
      attr.write(attr.getDataType(), buff);
    } else {
      auto attr = g.createAttribute(std::string{attr_name}, METADATA_STR_TYPE(), attr_space);
      attr.write(attr.getDataType(), buff);
    }
  } else {
    if (attribute_exists) {
      auto attr = g.openAttribute(std::string{attr_name});
      attr.write(attr.getDataType(), &buff);
    } else {
      auto attr = g.createAttribute(std::string{attr_name}, getH5_type<T>(), attr_space);
      attr.write(attr.getDataType(), &buff);
    }
  }
}

template <class T>
bool check_dataset_type(const H5::DataSet &dataset, T type, bool throw_on_failure) {
  static_assert(std::is_base_of_v<H5::DataType, T>,
                "type should have a type that is derived from H5::DataType (such as PredType, "
                "IntType or StrType).");
  const auto lck = internal::lock();
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

// Internal functions

[[nodiscard]] inline H5::StrType METADATA_STR_TYPE() {
  H5::StrType t(H5::PredType::C_S1, H5T_VARIABLE);
  t.setStrpad(H5T_STR_NULLTERM);
  t.setCset(H5T_CSET_UTF8);
  return t;
}

template <class H5Type>
[[nodiscard]] attr_types get_cpp_arithmetic_type(const H5Type &h5_type) {
  static_assert(std::is_base_of_v<H5Type, H5::IntType> || std::is_base_of_v<H5Type, H5::FloatType>);

  attr_types type;
  const auto size = h5_type.getSize();
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  if constexpr (std::is_same_v<H5Type, H5::IntType>) {
    assert(h5_type.getSign() != H5T_SGN_ERROR);
    const bool is_signed = h5_type.getSign() == H5T_SGN_NONE;
    assert(size >= 1 && size <= 8);
    if (is_signed) {
      switch (size) {
        case (sizeof(i8)):
          type = static_cast<i8>(0);
          return type;
        case (sizeof(i16)):
          type = static_cast<i16>(0);
          return type;
        case (sizeof(i32)):
          type = static_cast<i32>(0);
          return type;
        case (sizeof(i64)):
          type = static_cast<i64>(0);
          return type;
        default:
          MODLE_UNREACHABLE_CODE;
      }
    }
    switch (size) {
      case (sizeof(u8)):
        type = static_cast<u8>(0);
        return type;
      case (sizeof(u16)):
        type = static_cast<u16>(0);
        return type;
      case (sizeof(u32)):
        type = static_cast<u32>(0);
        return type;
      case (sizeof(u64)):
        type = static_cast<u64>(0);
        return type;
      default:
        MODLE_UNREACHABLE_CODE;
    }
  } else {
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
        MODLE_UNREACHABLE_CODE;
    }
  }
  MODLE_UNREACHABLE_CODE;
  DISABLE_WARNING_POP
}

template <class DataType>
inline H5::PredType getH5_type() {
  // Taken from github.com/DavidAce/h5pp:
  // https://github.com/DavidAce/h5pp/blob/master/h5pp/include/h5pp/details/h5ppUtils.h
  // clang-format off
  using DecayType = typename std::decay_t<DataType>;
  if constexpr (std::is_same_v<DecayType, short>) return H5::PredType::NATIVE_SHORT;
  if constexpr (std::is_same_v<DecayType, int>) return H5::PredType::NATIVE_INT;
  if constexpr (std::is_same_v<DecayType, long>) return H5::PredType::NATIVE_LONG;
  if constexpr (std::is_same_v<DecayType, long long>) return H5::PredType::NATIVE_LLONG;
  if constexpr (std::is_same_v<DecayType, unsigned short>) return H5::PredType::NATIVE_USHORT;
  if constexpr (std::is_same_v<DecayType, unsigned int>) return H5::PredType::NATIVE_UINT;
  if constexpr (std::is_same_v<DecayType, unsigned long>) return H5::PredType::NATIVE_ULONG;
  if constexpr (std::is_same_v<DecayType, unsigned long long>) return H5::PredType::NATIVE_ULLONG;
  if constexpr (std::is_same_v<DecayType, double>) return H5::PredType::NATIVE_DOUBLE;
  if constexpr (std::is_same_v<DecayType, long double>) return H5::PredType::NATIVE_LDOUBLE;
  if constexpr (std::is_same_v<DecayType, float>) return H5::PredType::NATIVE_FLOAT;
  if constexpr (std::is_same_v<DecayType, bool>) return H5::PredType::NATIVE_UINT8;
  if constexpr (std::is_same_v<DecayType, std::string>) return H5::PredType::C_S1;
  if constexpr (std::is_same_v<DecayType, char>) return H5::PredType::C_S1;
  // clang-format on

  throw std::logic_error("getH5_type(): Unable to map C++ type to a H5T");
}

template <class DType>
H5::DataSet create_dataset(H5::H5File &f, const std::string &name, const DType &type,
                           const H5::DSetCreatPropList &create_prop,
                           const H5::DSetAccPropList &access_prop) {
  const hsize_t RANK{1};                 // i.e. number of dimensions
  const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  const hsize_t BUFF_SIZE{1};            // Dummy buffer size
  const auto lck = internal::lock();
  auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};
  return f.createDataSet(name, type, mem_space, create_prop, access_prop);
}

template <class T1, class T2>
H5::DSetCreatPropList generate_creat_prop_list(const hsize_t chunk_size, const u8 compression_lvl,
                                               [[maybe_unused]] const T1 type,
                                               [[maybe_unused]] const T2 fill_value) {
  const auto lck = internal::lock();
  H5::DSetCreatPropList prop{};
  prop.setChunk(1, &chunk_size);
  prop.setDeflate(compression_lvl);
  if constexpr (!std::is_constructible_v<H5std_string, T2>) {
    prop.setFillValue(type, &fill_value);
  }

  return prop;
}

template <class T>
H5::DSetAccPropList generate_acc_prop_list([[maybe_unused]] T type, hsize_t chunk_size,
                                           hsize_t cache_size, const double rdcc_w0,
                                           const double multiplier) {
  // https://support.hdfgroup.org/HDF5/doc/Advanced/Chunking/index.html
  // The w0 parameter affects how the library decides which chunk to evict when it needs room in
  // the cache. If w0 is set to 0, then the library will always evict the least recently used
  // chunk in cache. If w0 is set to 1, the library will always evict the least recently used
  // chunk which has been fully read or written, and if none have been fully read or written, it
  // will evict the least recently used chunk. If w0 is between 0 and 1, the behavior will be a
  // blend of the two.

  // Therefore, if the application will access the same data more than once, w0 should be set
  // closer to 0, and if the application does not, w0 should be set closer to 1.
  const auto nslots = static_cast<usize>(multiplier * static_cast<double>(cache_size / chunk_size));
  const auto lck = internal::lock();
  H5::DSetAccPropList prop{};
  prop.setChunkCache(nslots, cache_size, rdcc_w0);
  return prop;
}
}  // namespace modle::hdf5
