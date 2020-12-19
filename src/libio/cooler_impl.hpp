#pragma once

#include <H5Cpp.h>
#include <absl/time/clock.h>

#include <algorithm>
#include <range/v3/action/remove_if.hpp>
#include <range/v3/action/transform.hpp>
#include <range/v3/algorithm/transform.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/generate.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <string_view>
#include <type_traits>
#include <vector>

#include "modle/contacts.hpp"
#include "modle/suppress_compiler_warnings.hpp"

namespace modle::cooler {

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

template <typename S>
hsize_t write_str(const S &str, H5::H5File &f, std::string_view dataset_name,
                  hsize_t MAX_STR_LENGTH, hsize_t file_offset, hsize_t CHUNK_DIMS,
                  uint8_t COMPRESSION_LEVEL) {
  static_assert(std::is_convertible<S, H5std_string>::value,
                "S should be a type convertible to std::string.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr const hsize_t RANK{1};     // i.e. number of dimensions
  constexpr const hsize_t BUFF_SIZE{1};
  constexpr const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
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
    dataset.write(str.c_str(), STD_STR, mem_space, file_space);

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
  constexpr const hsize_t RANK{1};     // i.e. number of dimensions
  constexpr const hsize_t BUFF_SIZE{1};
  constexpr const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
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
  constexpr const hsize_t RANK{1};     // i.e. number of dimensions
  const hsize_t BUFF_SIZE{1};
  constexpr const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  constexpr const I FILL_VAL{0};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(getH5_type<I>(), &FILL_VAL);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, getH5_type<I>(), mem_space, cparms);

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
hsize_t write_vect_of_int(std::vector<I> &data, H5::H5File &f, std::string_view dataset_name,
                          hsize_t file_offset, hsize_t CHUNK_DIMS, uint8_t COMPRESSION_LEVEL) {
  static_assert(std::is_integral<I>::value, "I should be an integer type.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr const hsize_t RANK{1};     // i.e. number of dimensions
  const hsize_t BUFF_SIZE{data.size()};
  constexpr const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  constexpr const I FILL_VAL{0};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(getH5_type<I>(), &FILL_VAL);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, getH5_type<I>(), mem_space, cparms);

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

}  // namespace modle::cooler