#pragma once

#include <H5Cpp.h>

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "modle/contacts.hpp"

namespace modle::hdf5 {

[[nodiscard]] inline std::string construct_error_stack();

template <typename S>
[[nodiscard]] inline hsize_t write_str(const S &str, const H5::DataSet &dataset,
                                       const H5::StrType &str_type, hsize_t file_offset);
template <typename CS>
[[nodiscard]] inline hsize_t write_strings(const CS &strings, const H5::DataSet &dataset,
                                           const H5::StrType &str_type, hsize_t file_offset);

template <typename N>
[[nodiscard]] inline hsize_t write_number(N &num, const H5::DataSet &dataset, hsize_t file_offset);
template <typename CN>
[[nodiscard]] inline hsize_t write_numbers(CN &numbers, const H5::DataSet &dataset,
                                           hsize_t file_offset);

template <typename N>
[[nodiscard]] inline hsize_t read_number(const H5::DataSet &dataset, N &buff, hsize_t file_offset);
template <typename CN>
[[nodiscard]] inline hsize_t read_numbers(const H5::DataSet &dataset, CN &buff,
                                          hsize_t file_offset);

[[nodiscard]] inline hsize_t read_str(const H5::DataSet &dataset, std::string &buff,
                                      hsize_t file_offset);
[[nodiscard]] inline hsize_t read_strings(const H5::DataSet &dataset,
                                          std::vector<std::string> &buff, hsize_t file_offset);
[[nodiscard]] inline std::string read_str(const H5::DataSet &dataset, hsize_t file_offset);
[[nodiscard]] inline std::vector<std::string> read_strings(const H5::DataSet &dataset,
                                                           hsize_t file_offset);

template <typename T>
inline void read_attribute(H5::H5File &f, std::string_view attr_name, T &buff,
                           std::string_view path = "/");
template <typename T>
inline void read_attribute(const H5::DataSet &d, std::string_view attr_name, T &buff);
template <typename T>
inline void read_attribute(std::string_view path_to_file, std::string_view attr_name, T &buff,
                           std::string_view path = "/");

[[nodiscard]] inline std::string read_attribute_str(H5::H5File &f, std::string_view attr_name,
                                                    std::string_view path = "/");
[[nodiscard]] inline std::string read_attribute_str(std::string_view path_to_file,
                                                    std::string_view attr_name,
                                                    std::string_view path = "/");

[[nodiscard]] inline int64_t read_attribute_int(H5::H5File &f, std::string_view attr_name,
                                                std::string_view path = "/");
[[nodiscard]] inline int64_t read_attribute_int(std::string_view path_to_file,
                                                std::string_view attr_name,
                                                std::string_view path = "/");

[[nodiscard]] inline bool has_attribute(const H5::Group &g, std::string_view attr_name);
[[nodiscard]] inline bool has_attribute(const H5::DataSet &d, std::string_view attr_name);
[[nodiscard]] inline bool has_attribute(H5::H5File &f, std::string_view attr_name,
                                        std::string_view path = "/");

template <typename T>
inline void write_or_create_attribute(H5::H5File &f, std::string_view attr_name, T &buff,
                                      std::string_view path = "/");

[[nodiscard]] inline bool has_group(H5::H5File &f, std::string_view name,
                                    std::string_view root_path = "/");
[[nodiscard]] inline bool has_dataset(H5::H5File &f, std::string_view name,
                                      std::string_view root_path = "/");

template <typename T>
[[nodiscard]] inline bool check_dataset_type(const H5::DataSet &dataset, T type,
                                             bool throw_on_failure = true);

[[nodiscard]] inline H5::H5File open_file_for_reading(std::string_view path_to_file);

}  // namespace modle::hdf5

#include "../../hdf5_impl.hpp"
