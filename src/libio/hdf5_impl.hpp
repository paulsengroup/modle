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
hsize_t write_str(const S &str, H5::DataSet &dataset, const H5::StrType &str_type,
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
    dataset.write(std::string{str}, str_type, mem_space, file_space);

    return file_size;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error occurred while writing a collectionof strings '{}' to dataset "
            "'{}' at offset {}:\n{}"),
        str, dataset.getObjName(), file_offset, construct_error_stack()));
  }
}

template <typename CS>
hsize_t write_strings(const CS &strings, H5::DataSet &dataset, const H5::StrType &str_type,
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
hsize_t write_number(N &num, H5::DataSet &dataset, hsize_t file_offset) {
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
hsize_t write_numbers(CN &numbers, H5::DataSet &dataset, hsize_t file_offset) {
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
hsize_t read_number(H5::DataSet &dataset, N &buff, hsize_t file_offset) {
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
hsize_t read_numbers(H5::DataSet &dataset, CN &buff, hsize_t file_offset) {
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

hsize_t read_str(H5::DataSet &dataset, std::string &buff, hsize_t file_offset) {
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

hsize_t read_strings(H5::DataSet &dataset, std::vector<std::string> &buff, hsize_t file_offset) {
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

std::string read_str(H5::DataSet &dataset, hsize_t file_offset) {
  std::string buff;
  (void)read_str(dataset, buff, file_offset);
  return buff;
}
std::vector<std::string> read_strings(H5::DataSet &dataset, hsize_t file_offset) {
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

bool has_attribute(H5::Group &g, std::string_view attr_name) {
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
bool check_dataset_type(H5::DataSet &dataset, T type, bool throw_on_failure) {
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

/*
template <typename S>
hsize_t write_str(const S &str, H5::H5File &f, std::string_view dataset_name,
                  hsize_t MAX_STR_LENGTH, hsize_t file_offset, hsize_t CHUNK_DIMS,
                  uint8_t COMPRESSION_LEVEL) {
  static_assert(std::is_convertible<S, H5std_string>::value,
                "S should be a type convertible to std::string.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr hsize_t RANK{1};     // i.e. number of dimensions
  constexpr hsize_t BUFF_SIZE{1};
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  const H5std_string FILL_VAL_STR{'\0'};

  try {
    auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    // Create the datatype as follows
    H5::StrType STD_STR(H5::PredType::C_S1, MAX_STR_LENGTH);
    STD_STR.setStrpad(H5T_STR_NULLPAD);  // This seems to be the most robust way of terminating
    // strings. Using NULLPAD sometimes causes some garbage to
    // appear at end of string
    STD_STR.setCset(H5T_CSET_ASCII);

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(STD_STR, &FILL_VAL_STR);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, STD_STR, mem_space, cparms);

    H5::DataSpace file_space;  // File space

    // This is probably not very efficient, but it should be fine, given that we don't expect to
    // write that many strings
    hsize_t file_size{file_offset + 1};
    dataset.extend(&file_size);
    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(std::string{str}, STD_STR, mem_space, file_space);

    return file_size;
  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::GroupIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSetIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSpaceIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}


template <typename S>
hsize_t write_vect_of_str(std::vector<S> &data, H5::H5File &f, std::string_view dataset_name,
                          hsize_t file_offset, hsize_t CHUNK_DIMS, uint8_t COMPRESSION_LEVEL) {
  static_assert(std::is_convertible<S, H5std_string>::value,
                "S should be a type convertible to std::string.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr hsize_t RANK{1};     // i.e. number of dimensions
  constexpr hsize_t BUFF_SIZE{1};
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  const H5std_string FILL_VAL_STR{'\0'};

  try {
    auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    const hsize_t MAX_STR_LENGTH =
        std::max_element(data.begin(), data.end(),
                         [](const auto &s1, const auto &s2) { return s1.size() < s2.size(); })
            ->size() +
        1;
    // Create the datatype as follows
    H5::StrType STD_STR(H5::PredType::C_S1, MAX_STR_LENGTH);
    STD_STR.setStrpad(H5T_STR_NULLTERM);  // This seems to be the most robust way of terminating
                                          // strings. Using NULLPAD sometimes causes some
                                          // garbage to appear at end of string
    STD_STR.setCset(H5T_CSET_ASCII);

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(STD_STR, &FILL_VAL_STR);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, STD_STR, mem_space, cparms);

    H5::DataSpace file_space;  // File space

    // This is probably not very efficient, but it should be fine, given that we don't expect to
    // write that many strings
    hsize_t file_size{file_offset + data.size()};
    dataset.extend(&file_size);
    for (const auto &s : data) {
      file_space = dataset.getSpace();
      file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
      dataset.write(s.c_str(), STD_STR, mem_space, file_space);
      ++file_offset;
    }
    return file_offset;
  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::GroupIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSetIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSpaceIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}


template <typename I>
hsize_t write_int(I num, H5::H5File &f, std::string_view dataset_name, hsize_t file_offset,
                  hsize_t CHUNK_DIMS, uint8_t COMPRESSION_LEVEL) {
  static_assert(std::is_integral<I>::value, "I should be an integer type.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr hsize_t RANK{1};     // i.e. number of dimensions
  const hsize_t BUFF_SIZE{1};
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  constexpr I FILL_VAL{0};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(getH5_type<I>(), &FILL_VAL);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, getH5_type<I>(), mem_space,
cparms);

    H5::DataSpace file_space;

    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);
    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(&num, getH5_type<I>(), mem_space, file_space);
    return file_size;

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::GroupIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSetIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSpaceIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

template <typename I>
hsize_t write_numbers(std::vector<I> &data, H5::DataSet &dataset, H5::DataSpace &mem_space,
                      hsize_t file_offset) {
  static_assert(std::is_integral<I>::value, "I should be an integer type.");

  const hsize_t BUFF_SIZE{data.size()};
  if (mem_space.getSimpleExtentNdims() != BUFF_SIZE) {
    constexpr hsize_t RANK{1};
    constexpr hsize_t MAXDIMS{H5S_UNLIMITED};
    mem_space = H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS);
  }

  try {
    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    const auto *cbuff = reinterpret_cast<I(*)[BUFF_SIZE]>(data.data());  // NOLINT
    DISABLE_WARNING_POP
    auto file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(cbuff, getH5_type<I>(), mem_space, file_space);
    return file_size;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

template <typename I>
hsize_t write_vect_of_int(std::vector<I> &data, H5::H5File &f, std::string_view dataset_name,
                          hsize_t file_offset, hsize_t CHUNK_DIMS, uint8_t COMPRESSION_LEVEL) {
  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr hsize_t RANK{1};     // i.e. number of dimensions
  const hsize_t BUFF_SIZE{data.size()};
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  constexpr I FILL_VAL{0};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(getH5_type<I>(), &FILL_VAL);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, getH5_type<I>(), mem_space,
cparms);

    return write_numbers(data, dataset, mem_space, file_offset);

  } catch (const H5::Exception &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}


template <typename I>
hsize_t write_numbers(std::vector<I> &data, H5::H5File &f, std::string_view dataset_name,
                          hsize_t file_offset, hsize_t CHUNK_DIMS, uint8_t COMPRESSION_LEVEL) {
  static_assert(std::is_integral<I>::value, "I should be an integer type.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr hsize_t RANK{1};     // i.e. number of dimensions
  const hsize_t BUFF_SIZE{data.size()};
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  constexpr I FILL_VAL{0};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(getH5_type<I>(), &FILL_VAL);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, getH5_type<I>(), mem_space,
cparms);

    H5::DataSpace file_space;

    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    const auto *cbuff = reinterpret_cast<I(*)[BUFF_SIZE]>(data.data());  // NOLINT
    DISABLE_WARNING_POP
    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(cbuff, getH5_type<I>(), mem_space, file_space);
    return file_size;

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::GroupIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSetIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSpaceIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}


template <typename I>
hsize_t read_int(H5::H5File &f, std::string_view dataset_name, I &BUFF, hsize_t file_offset) {
  static_assert(std::is_integral<I>::value, "I should be an integer type.");
  try {
    H5::Exception::dontPrint();

    const H5std_string d{dataset_name};
    constexpr hsize_t DIMS = 1;
    constexpr hsize_t RANK = 1;

    auto dataset = f.openDataSet(d);
    const auto DTYPE = dataset.getDataType();
    auto file_space = dataset.getSpace();

    auto mem_space{H5::DataSpace(RANK, &DIMS)};

    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &RANK, &file_offset);
    dataset.read(&BUFF, DTYPE, mem_space, file_space);

    return file_offset + 1;

  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to read data into an integer variable (dataset = {}; file offset = "
                   "{}, buff size = 1): Function '{}' failed with error '{}'"),
        dataset_name, file_offset, e.getFuncName(), e.getDetailMsg()));
  }
}

template <typename I>
hsize_t read_vect_of_int(H5::H5File &f, std::string_view dataset_name, std::vector<I> &BUFF,
                         hsize_t file_offset) {
  static_assert(std::is_integral<I>::value, "I should be an integer type.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr hsize_t RANK{1};           // i.e. number of dimensions
  hsize_t BUFF_SIZE{BUFF.size()};

  try {
    auto dataset = f.openDataSet(d);
    auto file_space = dataset.getSpace();

    if (const auto n_items = static_cast<hsize_t>(file_space.getSimpleExtentNpoints());
        n_items < BUFF_SIZE || BUFF_SIZE == 0) {
      BUFF_SIZE = n_items;
      BUFF.resize(BUFF_SIZE);
    }

    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE)};

    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    auto *cbuff = reinterpret_cast<I(*)[BUFF_SIZE]>(BUFF.data());  // NOLINT
    DISABLE_WARNING_POP

    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.read(cbuff, getH5_type<I>(), mem_space, file_space);

    return file_offset + BUFF_SIZE;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to read data into a vector of integers (dataset = {}; file offset = "
                   "{}; buff size = {}): Function '{}' failed with error '{}'"),
        dataset_name, file_offset, BUFF_SIZE, e.getFuncName(), e.getDetailMsg()));
  }
}

template <typename I>
void write_modle_cmatrix_to_cooler(const ContactMatrix<I> &cmatrix, std::string_view chr_name,
                                   uint64_t chr_start, uint64_t chr_end, uint64_t chr_length,
                                   std::string_view output_file, bool force_overwrite) {
  const std::vector<ContactMatrix<I>> cmatrices{cmatrix};
  write_modle_cmatrices_to_cooler(cmatrices, std::vector < std::string_view{chr_name},
                                std::vector<uint64_t>{chr_start},
std::vector<uint64_t>{chr_end}, std::vector<uint64_t>{chr_length}, output_file,
force_overwrite);
}


template <typename I>
void write_modle_cmatrices_to_cooler(const std::vector<ContactMatrix<I>> &cmatrices,
                                     const std::vector<std::string> &chr_names,
                                     const std::vector<uint64_t> &chr_starts,
                                     const std::vector<uint64_t> &chr_ends,
                                     const std::vector<uint64_t> &chr_sizes, uint64_t bin_size,
                                     std::string_view output_file, bool force_overwrite) {
  std::vector<const ContactMatrix<I> *> v(cmatrices.size());
  std::transform(cmatrices.begin(), cmatrices.end(), v.begin(), [](const auto &m) { return &m; });
  write_modle_cmatrices_to_cooler(v, chr_names, chr_starts, chr_ends, chr_sizes, bin_size,
                                  output_file, force_overwrite);
}

template <typename I>
void write_modle_cmatrices_to_cooler(const std::vector<const ContactMatrix<I> *> &cmatrices,
                                     const std::vector<std::string> &chr_names,
                                     const std::vector<uint64_t> &chr_starts,
                                     const std::vector<uint64_t> &chr_ends,
                                     const std::vector<uint64_t> &chr_sizes, uint64_t bin_size,
                                     std::string_view output_file, bool force_overwrite) {
  static_assert(std::is_integral<I>::value, "I should be an integral type.");

  fmt::print(stderr, FMT_STRING("Writing {} contact matrices for to file '{}'...\n"),
             cmatrices.size(), output_file);

  const auto n_chromosomes = cmatrices.size();
  if (n_chromosomes != chr_names.size() || n_chromosomes != chr_starts.size() ||
      n_chromosomes != chr_ends.size() || n_chromosomes != chr_sizes.size()) {
    utils::throw_with_trace(std::runtime_error(fmt::format(
        FMT_STRING("Message for the DEVS: The vectors passed to function "
                   "cooler::write_modle_cmatrices_to_cooler() should all have "
                   "the same size!\n  cmatrices.size()={}\n  chr_names.size()={}\n  "
                   "chr_starts={}\n  chr_ends={}\n  chr_sizes={}\n"),
        cmatrices.size(), chr_names.size(), chr_starts.size(), chr_ends.size(),
chr_sizes.size())));
  }

  auto t0 = absl::Now();

  if (force_overwrite && std::filesystem::exists(output_file)) {
    std::filesystem::remove(output_file);
  }

  const auto max_chrom_name_length =
      std::max_element(chr_names.begin(), chr_names.end(),
                       [](const auto &s1, const auto &s2) { return s1.size() < s2.size(); })
          ->size() +
      1;

  std::filesystem::create_directories(std::filesystem::path(output_file).parent_path());

  auto f = cooler::init_cooler_file(output_file, force_overwrite);
  auto datasets = cooler::init_cooler_datasets(f, max_chrom_name_length);

  std::size_t chr_name_foffset = 0;
  std::size_t chr_length_foffset = 0;

  constexpr std::size_t BUFF_SIZE = 1024 * 1024 / sizeof(int64_t);  // 1 MB
  std::vector<int64_t> bin_pos_buff;
  std::vector<int32_t> bin_chrom_buff;
  std::vector<int64_t> pixel_b1_idx_buff;
  std::vector<int64_t> pixel_b2_idx_buff;
  std::vector<int64_t> pixel_count_buff;
  std::vector<int64_t> idx_bin1_offset_buff;
  std::vector<int64_t> idx_chrom_offset_buff;

  bin_pos_buff.reserve(BUFF_SIZE);
  bin_chrom_buff.reserve(BUFF_SIZE);
  pixel_b1_idx_buff.reserve(BUFF_SIZE);
  pixel_b2_idx_buff.reserve(BUFF_SIZE);
  pixel_count_buff.reserve(BUFF_SIZE);
  idx_bin1_offset_buff.reserve(BUFF_SIZE);
  idx_chrom_offset_buff.reserve(cmatrices.size());

  hsize_t pixel_b1_idx_h5df_offset = 0;
  hsize_t pixel_b2_idx_h5df_offset = 0;
  hsize_t pixel_count_h5df_offset = 0;
  hsize_t idx_bin1_offset_h5df_offset = 0;

  auto write_pixels_to_file = [&]() {
    pixel_b1_idx_h5df_offset =
        cooler::write_vect_of_int(pixel_b1_idx_buff, f, "pixels/bin1_id",
pixel_b1_idx_h5df_offset); pixel_b2_idx_h5df_offset = cooler::write_vect_of_int(pixel_b2_idx_buff,
f, "pixels/bin2_id", pixel_b2_idx_h5df_offset); pixel_count_h5df_offset =
        cooler::write_vect_of_int(pixel_count_buff, f, "pixels/count", pixel_count_h5df_offset);

    pixel_b1_idx_buff.clear();
    pixel_b2_idx_buff.clear();
    pixel_count_buff.clear();
  };

  int64_t nnz = 0;
  int64_t offset = 0;
  int64_t nbins = 0;

  cooler::write_metadata(f, static_cast<int32_t>(bin_size));

  for (auto chr_idx = 0UL; chr_idx < n_chromosomes; ++chr_idx) {
    const auto &cmatrix = cmatrices[chr_idx];
    const auto &chr_name = chr_names[chr_idx];
    const auto &chr_total_len = chr_sizes[chr_idx];
    const auto &chr_start = chr_starts[chr_idx];
    const auto &chr_end = chr_ends[chr_idx];
    const auto chr_simulated_len = chr_end - chr_start;
    chr_name_foffset =
        cooler::write_str(chr_name, f, "chroms/name", max_chrom_name_length, chr_name_foffset);
    chr_length_foffset = cooler::write_int(chr_total_len, f, "chroms/length", chr_length_foffset);
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_CONVERSION
    DISABLE_WARNING_SIGN_CONVERSION
    fmt::print(stderr, FMT_STRING("Writing contacts for '{}' ({:.2f} Mbp)..."), chr_name,
               chr_simulated_len / 1.0e6);  // NOLINT
    const auto t1 = absl::Now();
    idx_chrom_offset_buff.push_back(nbins);
    nbins = cooler::write_bins(f, chr_idx++, chr_simulated_len, bin_size,  // NOLINT
                               bin_chrom_buff, bin_pos_buff, nbins);
    DISABLE_WARNING_POP

    for (auto i = 0UL; i < chr_start / bin_size; ++i) {
      idx_bin1_offset_buff.push_back(nnz);
      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {
        idx_bin1_offset_h5df_offset = cooler::write_vect_of_int(
            idx_bin1_offset_buff, f, "indexes/bin1_offset", idx_bin1_offset_h5df_offset);
        idx_bin1_offset_buff.clear();
      }
    }

    offset += static_cast<int64_t>(chr_start / bin_size);
    for (auto i = 0UL; i < cmatrix->n_cols(); ++i) {
      idx_bin1_offset_buff.push_back(nnz);
      for (auto j = i; j < i + cmatrix->n_rows() && j < cmatrix->n_cols(); ++j) {
        if (const auto m = cmatrix->get(i, j); m != 0) {
          pixel_b1_idx_buff.push_back(offset + static_cast<int64_t>(i));
          pixel_b2_idx_buff.push_back(offset + static_cast<int64_t>(j));
#ifndef NDEBUG
          if (offset + static_cast<int64_t>(i) > offset + static_cast<int64_t>(j)) {
            utils::throw_with_trace(std::runtime_error(
                fmt::format(FMT_STRING("b1 > b2: b1={}; b2={}; offset={}; m={}\n"),
                            pixel_b1_idx_buff.back(), pixel_b2_idx_buff.back(), offset, m)));
          }
#endif
          pixel_count_buff.push_back(m);
          ++nnz;
          if (pixel_b1_idx_buff.size() == BUFF_SIZE) {
            write_pixels_to_file();
          }
        }
      }
      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {
        idx_bin1_offset_h5df_offset = cooler::write_vect_of_int(
            idx_bin1_offset_buff, f, "indexes/bin1_offset", idx_bin1_offset_h5df_offset);
        idx_bin1_offset_buff.clear();
      }
    }
    for (auto i = 0UL; i < (chr_total_len - chr_end) / bin_size; ++i) {
      idx_bin1_offset_buff.push_back(nnz);
      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {
        idx_bin1_offset_h5df_offset = cooler::write_vect_of_int(
            idx_bin1_offset_buff, f, "indexes/bin1_offset", idx_bin1_offset_h5df_offset);
        idx_bin1_offset_buff.clear();
      }
    }

    offset = nbins;
    fmt::print(stderr, FMT_STRING(" DONE in {}!\n"), absl::FormatDuration(absl::Now() - t1));
  }

  if (!pixel_b1_idx_buff.empty()) {
    write_pixels_to_file();
  }

  idx_chrom_offset_buff.push_back(nbins);
  idx_bin1_offset_buff.push_back(nnz);

  cooler::write_vect_of_int(idx_chrom_offset_buff, f, "indexes/chrom_offset");
  cooler::write_vect_of_int(idx_bin1_offset_buff, f, "indexes/bin1_offset",
                            idx_bin1_offset_h5df_offset);

  H5::DataSpace att_space(H5S_SCALAR);
  auto att = f.createAttribute("nchroms", H5::PredType::NATIVE_INT64, att_space);
  att.write(H5::PredType::NATIVE_INT, &n_chromosomes);
  att = f.createAttribute("nbins", H5::PredType::NATIVE_INT64, att_space);
  att.write(H5::PredType::NATIVE_INT, &nbins);
  att = f.createAttribute("nnz", H5::PredType::NATIVE_INT64, att_space);
  att.write(H5::PredType::NATIVE_INT, &nnz);

  fmt::print(stderr, "DONE! Saved {} pixels in {}.\n", nnz, absl::FormatDuration(absl::Now() -
t0));
}

hsize_t write_vect_of_enums(std::vector<int32_t> &data, const H5::EnumType &ENUM, H5::H5File &f,
                            std::string_view dataset_name, hsize_t file_offset, hsize_t
CHUNK_DIMS, uint8_t COMPRESSION_LEVEL) { const H5std_string d{dataset_name};  // Probably
unnecessary constexpr hsize_t RANK{1};           // i.e. number of dimensions const hsize_t
BUFF_SIZE{data.size()}; constexpr hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset = f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, ENUM, mem_space,
cparms);

    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);
    auto file_space = dataset.getSpace();
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    const auto *cbuff = reinterpret_cast<int32_t(*)[BUFF_SIZE]>(data.data());  // NOLINT
    DISABLE_WARNING_POP
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(cbuff, ENUM, mem_space, file_space);
    return file_size;

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::GroupIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSetIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSpaceIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

H5::H5File init_cooler_file(std::string_view path_to_output_file, bool force_overwrite) {
  H5::H5File f(std::string{path_to_output_file}, force_overwrite ? H5F_ACC_TRUNC : H5F_ACC_CREAT);

  for (const auto &g : {"chroms", "bins", "pixels", "indexes"}) {
    f.createGroup(g);
  }
  return f;
}

H5::EnumType init_enum_from_strs(const std::vector<std::string> &data, int32_t offset) {
  H5::EnumType ENUM{getH5_type<int32_t>()};
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_COMPARE
  DISABLE_WARNING_SIGN_CONVERSION
  for (int32_t i = offset; i < data.size(); ++i) {
    ENUM.insert(data[i], &i);
  }
  DISABLE_WARNING_POP
  return ENUM;
}

std::vector<H5::DataSet> init_cooler_datasets(H5::H5File &f, hsize_t MAX_STR_LENGTH,
                                              hsize_t CHUNK_DIMS, uint8_t COMPRESSION_LEVEL) {
  constexpr hsize_t RANK{1};                 // i.e. number of dimensions
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  constexpr hsize_t BUFF_SIZE{1};            // Dummy buffer size
  constexpr int64_t FILL_VAL_INT64{0};
  constexpr int32_t FILL_VAL_INT32{0};
  const H5std_string FILL_VAL_STR{'\0'};
  std::vector<H5::DataSet> datasets;
  datasets.reserve(10);  // NOLINT

  try {
    H5::StrType STD_STR(H5::PredType::C_S1, MAX_STR_LENGTH);
    STD_STR.setStrpad(H5T_STR_NULLPAD);
    STD_STR.setCset(H5T_CSET_ASCII);
    auto &INT64 = H5::PredType::NATIVE_INT64;
    auto &INT32 = H5::PredType::NATIVE_INT;

    H5::DSetCreatPropList cparms_str{};
    cparms_str.setChunk(RANK, &CHUNK_DIMS);
    cparms_str.setDeflate(COMPRESSION_LEVEL);

    H5::DSetCreatPropList cparms_int64 = cparms_str;
    H5::DSetCreatPropList cparms_int32 = cparms_str;
    cparms_str.setFillValue(STD_STR, &FILL_VAL_STR);
    cparms_int64.setFillValue(INT64, &FILL_VAL_INT64);
    cparms_int64.setFillValue(INT32, &FILL_VAL_INT32);

    // rdcc_w0 == 1 means that HDf5lib will evict from cache the least recently
    // used chunk which has been fully read or written);
    H5::DSetAccPropList aparms_int{};
    aparms_int.setChunkCache(100 * 10, 10 * ONE_MB, 1.0);

    auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};
    // Do not change the order of these pushbacks
    datasets.push_back(
        f.createDataSet("chroms/length", INT64, mem_space, cparms_int64, aparms_int));
    datasets.push_back(f.createDataSet("chroms/name", STD_STR, mem_space, cparms_str));

    datasets.push_back(f.createDataSet("bins/chrom", INT64, mem_space, cparms_int64, aparms_int));
    datasets.push_back(f.createDataSet("bins/start", INT64, mem_space, cparms_int64, aparms_int));
    datasets.push_back(f.createDataSet("bins/end", INT64, mem_space, cparms_int64, aparms_int));

    datasets.push_back(
        f.createDataSet("pixels/bin1_id", INT64, mem_space, cparms_int64, aparms_int));
    datasets.push_back(
        f.createDataSet("pixels/bin2_id", INT64, mem_space, cparms_int64, aparms_int));
    datasets.push_back(f.createDataSet("pixels/count", INT32, mem_space, cparms_int32,
aparms_int));

    datasets.push_back(
        f.createDataSet("indexes/bin1_offset", INT64, mem_space, cparms_int64, aparms_int));
    datasets.push_back(
        f.createDataSet("indexes/chrom_offset", INT64, mem_space, cparms_int64, aparms_int));

    return datasets;

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

std::vector<H5::DataSet> open_cooler_datasets(H5::H5File &f, hsize_t CACHE_INT64_NBYTES,
                                              hsize_t CACHE_INT32_NBYTES, hsize_t
CACHE_STR_NBYTES, double CACHE_RDCC_W0) { std::vector<H5::DataSet> datasets; datasets.reserve(10);
// NOLINT

  try {
    H5::DSetAccPropList aparms_str{};
    aparms_str.setChunkCache(100 * (CACHE_STR_NBYTES / (ONE_MB / 2)), CACHE_STR_NBYTES,
                             CACHE_RDCC_W0);
    H5::DSetAccPropList aparms_int64{};
    aparms_int64.setChunkCache(100 * 64, 64 * ONE_MB, CACHE_RDCC_W0);
    H5::DSetAccPropList aparms_int32{};
    aparms_int32.setChunkCache(100 * 32, 32 * ONE_MB, CACHE_RDCC_W0);

    // Do not change the order of these pushbacks
    datasets.push_back(f.openDataSet("chroms/length", aparms_int64));
    datasets.push_back(f.openDataSet("chroms/name", aparms_str));

    datasets.push_back(f.openDataSet("bins/chrom", aparms_int64));
    datasets.push_back(f.openDataSet("bins/start", aparms_int64));
    datasets.push_back(f.openDataSet("bins/end", aparms_int64));

    datasets.push_back(f.openDataSet("pixels/bin1_id", aparms_int64));
    datasets.push_back(f.openDataSet("pixels/bin2_id", aparms_int64));
    datasets.push_back(f.openDataSet("pixels/count", aparms_int32));

    datasets.push_back(f.openDataSet("indexes/bin1_offset", aparms_int64));
    datasets.push_back(f.openDataSet("indexes/chrom_offset", aparms_int64));

    return datasets;

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

hsize_t write_bins(H5::H5File &f, int32_t chrom, int64_t length, int64_t bin_size,
                   std::vector<int32_t> &buff32, std::vector<int64_t> &buff64, hsize_t
file_offset, hsize_t BUFF_SIZE) { buff32.resize(BUFF_SIZE); buff64.resize(BUFF_SIZE);

  std::fill(buff32.begin(), buff32.end(), chrom);

  int64_t start = 0;
  int64_t end = std::min(length, bin_size);
  hsize_t bins_processed = 0;
  const int64_t nbins = (length / bin_size) + 1;
  const int64_t nchunks = (nbins / static_cast<int64_t>(BUFF_SIZE)) + 1;

  for (auto i = 0; i < nchunks; ++i) {
    const auto chunk_size = std::min(BUFF_SIZE, static_cast<hsize_t>(nbins) - bins_processed);
    if (chunk_size != BUFF_SIZE) {
      buff32.resize(chunk_size);
      buff64.resize(chunk_size);
    }
    write_vect_of_int(buff32, f, "bins/chrom", file_offset);

    std::generate(buff64.begin(), buff64.end(), [&]() { return (start += bin_size) - bin_size; });
    write_vect_of_int(buff64, f, "bins/start", file_offset);

    std::generate(buff64.begin(), buff64.end(), [&]() { return (end += bin_size) - bin_size; });
    if (chunk_size != BUFF_SIZE) {
      buff64.back() = length;
    }
    file_offset = write_vect_of_int(buff64, f, "bins/end", file_offset);
    bins_processed += chunk_size;
  }
  assert(buff64.back() == length);

  return file_offset;
}

void write_metadata(H5::H5File &f, int32_t bin_size, std::string_view assembly_name) {
  H5::StrType STR_TYPE(H5::PredType::C_S1, H5T_VARIABLE);
  STR_TYPE.setCset(H5T_CSET_UTF8);
  H5::DataSpace att_space(H5S_SCALAR);
  int32_t int_buff{};
  std::string str_buff{};

  auto att = f.createAttribute("format", STR_TYPE, att_space);
  str_buff = "HDF5::Cooler";
  att.write(STR_TYPE, str_buff);

  att = f.createAttribute("format-version", H5::PredType::NATIVE_INT64, att_space);
  int_buff = 3;
  att.write(H5::PredType::NATIVE_INT, &int_buff);

  att = f.createAttribute("bin-type", STR_TYPE, att_space);
  str_buff = "fixed";
  att.write(STR_TYPE, str_buff);

  att = f.createAttribute("bin-size", H5::PredType::NATIVE_INT64, att_space);
  att.write(H5::PredType::NATIVE_INT, &bin_size);

  att = f.createAttribute("storage-mode", STR_TYPE, att_space);
  str_buff = "symmetric-upper";
  att.write(STR_TYPE, str_buff);

  if (!assembly_name.empty()) {
    str_buff = std::string{assembly_name};
    att = f.createAttribute("assembly-name", STR_TYPE, att_space);
    att.write(STR_TYPE, &str_buff);
  }

  att = f.createAttribute("generated-by", STR_TYPE, att_space);
  str_buff = fmt::format("ModLE-v{}", "0.0.1");  // TODO make ModLE ver a tunable
  att.write(STR_TYPE, str_buff);

  att = f.createAttribute("creation-date", STR_TYPE, att_space);
  str_buff = absl::FormatTime(absl::Now(), absl::UTCTimeZone());
  att.write(STR_TYPE, str_buff);
}

std::size_t read_chr_offset_idx(H5::H5File &f, std::string_view chr_name) {
  try {
    H5::Exception::dontPrint();

    constexpr hsize_t DIMS = 1;
    constexpr hsize_t RANK = 1;
    std::string BUFF;

    auto dataset = f.openDataSet("chroms/name");
    auto STR_TYPE = dataset.getDataType();
    auto file_space = dataset.getSpace();

    const auto nchroms = static_cast<hsize_t>(file_space.getSimpleExtentNpoints());

    auto mem_space{H5::DataSpace(RANK, &DIMS)};

    for (hsize_t i = 0U; i < nchroms; ++i) {
      file_space = dataset.getSpace();
      file_space.selectHyperslab(H5S_SELECT_SET, &RANK, &i);
      dataset.read(BUFF, STR_TYPE, mem_space, file_space);
      if (BUFF == chr_name) {
        return i;
      }
    }

    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find a chromosome named '{}'"), chr_name));
  } catch (const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while looking for chromosome "
                               "'{}': Function '{}' failed with error '{}'"),
                    chr_name, e.getFuncName(), e.getDetailMsg()));
  }
}

std::vector<int64_t> read_chr_offset_idx(H5::H5File &f) {
  std::vector<int64_t> idx;
  read_vect_of_int(f, "indexes/chrom_offset", idx, 0);
  return idx;
}

std::vector<int64_t> read_bin1_offset_idx(H5::H5File &f) {
  std::vector<int64_t> idx;
  read_vect_of_int(f, "indexes/bin1_offset", idx, 0);
  return idx;
}

std::vector<int64_t> read_bin1_offset_idx(H5::H5File &f, std::string_view chr_name) {
  const auto chr_idx = read_chr_offset_idx(f, chr_name);
  return read_bin1_offset_idx(f, chr_idx);
}

std::vector<int64_t> read_bin1_offset_idx(H5::H5File &f, std::size_t chr_idx) {
  const auto chr_offset_idx = read_chr_offset_idx(f);
  assert(chr_idx < chr_offset_idx.size());  // NOLINT
  const auto chr_start_bin = static_cast<std::size_t>(chr_offset_idx[chr_idx]);
  const auto chr_end_bin = static_cast<std::size_t>(chr_offset_idx[chr_idx + 1] - 1);
  assert(chr_end_bin >= chr_start_bin);  // NOLINT

  std::vector<int64_t> bin1_offset_idx(chr_end_bin - chr_start_bin);
  read_vect_of_int(f, "indexes/bin1_offset", bin1_offset_idx, chr_start_bin);

  return bin1_offset_idx;
}

std::pair<int64_t, int64_t> read_chrom_pixels_boundaries(H5::H5File &f, std::string_view chr_name)
{ const auto chr_idx = read_chr_offset_idx(f, chr_name); return read_chrom_pixels_boundaries(f,
chr_idx);
}

std::pair<int64_t, int64_t> read_chrom_pixels_boundaries(H5::H5File &f, std::size_t chr_idx) {
  const auto chr_offset_idx = read_chr_offset_idx(f);
  assert(chr_idx < chr_offset_idx.size());  // NOLINT
  const auto chr_start_bin = static_cast<std::size_t>(chr_offset_idx[chr_idx]);
  const auto chr_end_bin = static_cast<std::size_t>(chr_offset_idx[chr_idx + 1]);
  assert(chr_end_bin >= chr_start_bin);  // NOLINT

  std::pair<int64_t, int64_t> pixel_boundaries{};
  read_int(f, "indexes/bin1_offset", pixel_boundaries.first, chr_start_bin);
  read_int(f, "indexes/bin1_offset", pixel_boundaries.second, chr_end_bin);

  return pixel_boundaries;
}

ContactMatrix<uint32_t> cooler_to_cmatrix(H5::H5File &f, int64_t bin_offset,
                                          const std::vector<int64_t> &bin1_offset_idx,
                                          std::size_t diagonal_width, std::size_t bin_size) {
  // TODO: Add check for Cooler bin size
  const auto nrows = diagonal_width / bin_size;
  return cooler_to_cmatrix(f, bin_offset, bin1_offset_idx, nrows);
}

ContactMatrix<uint32_t> cooler_to_cmatrix(H5::H5File &f, int64_t bin_offset,
                                          const std::vector<int64_t> &bin1_offset_idx,
                                          std::size_t nrows) {
  ContactMatrix<uint32_t> cmatrix(nrows, bin1_offset_idx.size());
  std::vector<int64_t> bin1_BUFF(nrows);
  std::vector<int64_t> bin2_BUFF(nrows);
  std::vector<int64_t> count_BUFF(nrows);

  for (auto i = 1UL; i < bin1_offset_idx.size(); ++i) {
    const auto file_offset = static_cast<hsize_t>(bin1_offset_idx[i - 1]);
    const auto buff_size =  // Figure out why buff size is always small
        std::min(static_cast<std::size_t>(bin1_offset_idx[i] - bin1_offset_idx[i - 1]), nrows);
    if (buff_size == 0) {
      continue;
    }
    // fmt::print(stderr, "i={}; buff_size={}\n", i, buff_size);

    bin1_BUFF.resize(buff_size);
    bin2_BUFF.resize(buff_size);
    count_BUFF.resize(buff_size);

    read_vect_of_int(f, "pixels/bin1_id", bin1_BUFF, file_offset);
    read_vect_of_int(f, "pixels/bin2_id", bin2_BUFF, file_offset);
    read_vect_of_int(f, "pixels/count", count_BUFF, file_offset);

    assert(bin1_BUFF.size() == buff_size);   // NOLINT
    assert(bin2_BUFF.size() == buff_size);   // NOLINT
    assert(count_BUFF.size() == buff_size);  // NOLINT

    for (auto j = 0UL; j < buff_size; ++j) {
      assert(count_BUFF[j] != 0);  // NOLINT
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      DISABLE_WARNING_SIGN_COMPARE
      if (const auto bin2_id = bin2_BUFF[j] - bin_offset;
          bin2_id >= i + nrows - 1 || bin2_id >= bin1_offset_idx.size()) {
        break;
      }
      // fmt::print(stderr, "m[{}][{}]={} (nrows={}; ncols={})\n", bin2_BUFF[j] - bin_offset,
      //           bin1_BUFF[j] - bin_offset, count_BUFF[j], cmatrix.n_rows(),
      //           cmatrix.n_cols());
      cmatrix.set(bin2_BUFF[j] - bin_offset, bin1_BUFF[j] - bin_offset, count_BUFF[j]);

      DISABLE_WARNING_POP
    }
  }
  return cmatrix;
}

ContactMatrix<uint32_t> cooler_to_cmatrix(std::string_view path_to_file, std::string_view
chr_name, std::size_t diagonal_width, std::size_t bin_size, bool try_common_chr_prefixes) {
  // TODO: Add check for Cooler bin size
  const auto nrows =
      (diagonal_width / bin_size) + static_cast<std::size_t>(diagonal_width % bin_size != 0);
  return cooler_to_cmatrix(path_to_file, chr_name, nrows, try_common_chr_prefixes);
}

ContactMatrix<uint32_t> cooler_to_cmatrix(std::string_view path_to_file, std::string_view
chr_name, std::size_t nrows, bool try_common_chr_prefixes) {
  // TODO: Add check for Cooler bin size
  auto search_chr_offset_idx = [&](H5::H5File &f) {
    if (!try_common_chr_prefixes) {
      return std::make_pair(chr_name, cooler::read_chr_offset_idx(f, chr_name));
    }
    // Here's the issue: for a given genome assembly (say hg19), some tools name chromosome as
    // chrN, where N is the chromosome number, while others only use the number N. The purpose
    // of this lambda is to try to guess few reasonably common prefixes, look them up in the
    // Cooler file, then return the chromosome name variant that produced a hit together with
    // the chr index
    std::vector<std::string_view> chr_names{chr_name};
    for (const auto &prefix : {"chr", "CHR", "Chr"}) {
      if (absl::StartsWith(chr_name, prefix)) {
        chr_names.emplace_back(absl::StripPrefix(chr_name, prefix));
      } else {
        chr_names.emplace_back(absl::StrCat(prefix, chr_name));
      }
    }
    // TODO: Consider reading the full array of chr names and then look up all the combinations.
    // This would allows us to avoid using exceptions, but will increase memory usage when
    // processing highly fragmented assemblies
    chr_name = "";  // Just to make sure we are not using this variable instead of actual_chr_name
    for (const auto name : chr_names) {
      try {
        return std::make_pair(name, cooler::read_chr_offset_idx(f, name));
      } catch (const std::runtime_error &e) {
        if (!absl::StartsWith(e.what(), "Unable to find a chromosome named")) {
          throw;
        }
      }
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find a chromosome named '{}'. The following chromosome "
                               "name variants were searched: '{}'"),
                    chr_name, absl::StrJoin(chr_names, "', '")));
  };

  auto f = open_for_reading(path_to_file);
  const auto [actual_chr_name, chr_idx] = search_chr_offset_idx(f);
  const auto bin1_offset_idx = cooler::read_bin1_offset_idx(f, actual_chr_name);
  const auto bin_offset = cooler::read_chr_offset_idx(f)[chr_idx];

  return cooler_to_cmatrix(f, bin_offset, bin1_offset_idx, nrows);
}

std::string read_attribute_str(H5::H5File &f, std::string_view attr_name, std::string_view path) {
  std::string buff;
  read_attribute(f, attr_name, buff, path);
  return buff;
}

std::string read_attribute_str(std::string_view path_to_file, std::string_view attr_name,
                               std::string_view path) {
  std::string buff;
  auto f = open_for_reading(path_to_file);
  read_attribute(f, attr_name, buff, path);
  return buff;
}

int64_t read_attribute_int(H5::H5File &f, std::string_view attr_name, std::string_view path) {
  int64_t buff;
  read_attribute(f, attr_name, buff, path);
  return buff;
}

int64_t read_attribute_int(std::string_view path_to_file, std::string_view attr_name,
                           std::string_view path) {
  int64_t buff;
  auto f = open_for_reading(path_to_file);
  read_attribute(f, attr_name, buff, path);
  return buff;
}
 */

}  // namespace modle::hdf5