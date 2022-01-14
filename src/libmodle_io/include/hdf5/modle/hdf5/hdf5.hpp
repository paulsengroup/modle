// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

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
// IWYU pragma: no_include <H5public.h>
// IWYU pragma: no_include <H5StrType.h>

#include <H5Cpp.h>               // IWYU pragma: keep
#include <absl/types/variant.h>  // for variant

#include <boost/filesystem/path.hpp>  // for path
#include <mutex>                      // recursive_mutex, unique_lock
#include <string>                     // for string
#include <string_view>                // for string_view
#include <vector>                     // for vector

#include "modle/common/common.hpp"  // for i64

namespace modle::hdf5 {

namespace internal {
[[nodiscard]] std::scoped_lock<std::recursive_mutex> lock();
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
inline static std::recursive_mutex global_mtx;

}  // namespace internal

template <class DataType>
inline H5::PredType getH5_type();

using attr_types = absl::variant<u8, u16, u32, u64, i8, i16, i32, i64, float, double, long double>;
template <class H5Type>
[[nodiscard]] attr_types get_cpp_arithmetic_type(const H5Type &h5_type);

[[nodiscard]] H5::StrType METADATA_STR_TYPE();

[[nodiscard]] std::string construct_error_stack(std::string_view function_name = "",
                                                std::string_view detail_msg = "");

[[nodiscard]] std::string construct_error_stack(const H5::Exception &e);

template <class S>
[[nodiscard]] inline hsize_t write_str(const S &str, const H5::DataSet &dataset,
                                       const H5::StrType &str_type, hsize_t file_offset);
template <class CS>
[[nodiscard]] inline hsize_t write_strings(const CS &strings, const H5::DataSet &dataset,
                                           const H5::StrType &str_type, hsize_t file_offset,
                                           bool write_empty_strings = false);

template <class N>
[[nodiscard]] inline hsize_t write_number(N &num, const H5::DataSet &dataset, hsize_t file_offset);
template <class CN>
[[nodiscard]] inline hsize_t write_numbers(CN &numbers, const H5::DataSet &dataset,
                                           hsize_t file_offset);

template <class N>
[[nodiscard]] inline hsize_t read_number(const H5::DataSet &dataset, N &buff, hsize_t file_offset);
template <class CN>
[[nodiscard]] inline hsize_t read_numbers(const H5::DataSet &dataset, CN &buff,
                                          hsize_t file_offset);

[[nodiscard]] hsize_t read_str(const H5::DataSet &dataset, std::string &buff, hsize_t file_offset);
[[nodiscard]] hsize_t read_strings(const H5::DataSet &dataset, std::vector<std::string> &buff,
                                   hsize_t file_offset);
[[nodiscard]] std::string read_str(const H5::DataSet &dataset, hsize_t file_offset);
[[nodiscard]] std::vector<std::string> read_strings(const H5::DataSet &dataset,
                                                    hsize_t file_offset);

template <class T>
inline void read_attribute(H5::H5File &f, std::string_view attr_name, T &buff,
                           std::string_view path = "/");
template <class T>
inline void read_attribute(const H5::DataSet &d, std::string_view attr_name, T &buff);
template <class T>
inline void read_attribute(const boost::filesystem::path &path_to_file, std::string_view attr_name,
                           T &buff, std::string_view path = "/");

template <class T>
[[nodiscard]] inline T read_attribute(H5::H5File &f, std::string_view attr_name,
                                      std::string_view path = "/");
template <class T>
[[nodiscard]] inline T read_attribute(const H5::DataSet &d, std::string_view attr_name);
template <class T>
[[nodiscard]] inline T read_attribute(const boost::filesystem::path &path_to_file,
                                      std::string_view attr_name, std::string_view path = "/");

[[nodiscard]] std::string read_attribute_str(H5::H5File &f, std::string_view attr_name,
                                             std::string_view path = "/");
[[nodiscard]] std::string read_attribute_str(const boost::filesystem::path &path_to_file,
                                             std::string_view attr_name,
                                             std::string_view path = "/");

[[nodiscard]] i64 read_attribute_int(H5::H5File &f, std::string_view attr_name,
                                     std::string_view path = "/");
[[nodiscard]] i64 read_attribute_int(const boost::filesystem::path &path_to_file,
                                     std::string_view attr_name, std::string_view path = "/");

[[nodiscard]] bool has_attribute(const H5::Group &g, std::string_view attr_name);
[[nodiscard]] bool has_attribute(const H5::DataSet &d, std::string_view attr_name);
[[nodiscard]] bool has_attribute(H5::H5File &f, std::string_view attr_name,
                                 std::string_view path = "/");

template <class T>
inline void write_or_create_attribute(H5::H5File &f, std::string_view attr_name, T &buff,
                                      std::string_view path = "/");
template <class DType>
[[nodiscard]] inline H5::DataSet create_dataset(H5::H5File &f, const std::string &name,
                                                const DType &type,
                                                const H5::DSetCreatPropList &create_prop,
                                                const H5::DSetAccPropList &access_prop);

[[nodiscard]] H5::DataSet open_dataset(H5::H5File &f, const std::string &name,
                                       const H5::DSetAccPropList &access_prop);

[[nodiscard]] H5::Group open_group(H5::H5File &f, const std::string &name);
[[nodiscard]] H5::Group create_group(H5::H5File &f, const std::string &name);
[[nodiscard]] H5::Group open_or_create_group(H5::H5File &f, const std::string &name);

[[nodiscard]] bool has_group(H5::H5File &f, std::string_view name,
                             std::string_view root_path = "/");
[[nodiscard]] bool has_dataset(H5::H5File &f, std::string_view name,
                               std::string_view root_path = "/");

template <class T>
[[nodiscard]] bool check_dataset_type(const H5::DataSet &dataset, T type,
                                      bool throw_on_failure = true);

[[nodiscard]] H5::H5File open_file_for_reading(const boost::filesystem::path &path_to_file);
[[nodiscard]] H5::H5File open_file_for_writing(const boost::filesystem::path &path_to_file);

template <class T1, class T2>
[[nodiscard]] inline H5::DSetCreatPropList generate_creat_prop_list(hsize_t chunk_size,
                                                                    u8 compression_lvl, T1 type,
                                                                    T2 fill_value);
template <class T>
[[nodiscard]] inline H5::DSetAccPropList generate_acc_prop_list(T type, hsize_t chunk_size,
                                                                hsize_t cache_size, double rdcc_w0,
                                                                double multiplier = 100.0);

[[nodiscard]] std::string get_file_name(H5::H5File &f);

}  // namespace modle::hdf5

#include "../../../../hdf5_impl.hpp"  // IWYU pragma: export
