// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cpp-sort/comparators/natural_less.h>
#include <cpp-sort/sorters/insertion_sorter.h>
#include <fmt/format.h>   // for compile_string_to_view, FMT_STRING, formatbu...
#include <fmt/ostream.h>  // for compile_string_to_view, FMT_STRING, formatbu...

#include <algorithm>                         // for transform
#include <boost/filesystem/file_status.hpp>  // for file_type, regular_file, directory_file, fil...
#include <boost/filesystem/operations.hpp>   // for status
#include <boost/filesystem/path.hpp>         // for operator<<, path
#include <cctype>                            // for isalpha
#include <cmath>                             // for trunc
#include <exception>                         // for exception
#include <range/v3/view/any_view.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/zip.hpp>
#include <string>       // for string, basic_string
#include <string_view>  // for string_view

#include "modle/common/numeric_utils.hpp"  // for parse_numeric_or_throw

namespace modle::utils {

// Try to convert str representations like "1.0" or "1.000000" to "1"
std::string str_float_to_str_int(const std::string &s) {
  try {
    auto n = parse_numeric_or_throw<double>(s);
    if (std::trunc(n) == n) {
      return fmt::format(FMT_STRING("{:.0f}"), n);
    }
  } catch (const std::exception &e) {  // Let the callee deal with invalid numbers
    return s;
  }
  return s;
}

template <char replacement>
std::string replace_non_alpha_char(const std::string &s) {
  return replace_non_alpha_char(s, replacement);
}

std::string replace_non_alpha_char(const std::string &s, const char replacement) {
  auto ss = s;
  std::transform(ss.begin(), ss.end(), ss.begin(),
                 [&](const auto c) { return c ? std::isalpha(c) : replacement; });
  return ss;
}

template <class Collection>
std::string format_collection_to_english_list(const Collection &collection,
                                              const std::string_view sep,
                                              std::string_view last_sep) {
  if (last_sep.empty()) {
    last_sep = sep;
  }

  const auto size = std::distance(std::begin(collection), std::end(collection));
  if (size <= 2) {
    return fmt::format(FMT_STRING("{}"), fmt::join(collection, last_sep));
  }

  const auto first = std::begin(collection);
  const auto second_to_last = std::end(collection) - 1;
  return fmt::format(FMT_STRING("{}{}{}"), fmt::join(first, second_to_last, sep), last_sep,
                     *second_to_last);
}

std::string detect_path_collision(const boost::filesystem::path &p, bool force_overwrite,
                                  boost::filesystem::file_type expected_type) {
  std::string error_msg;
  detect_path_collision(p, error_msg, force_overwrite, expected_type);
  return error_msg;
}

bool detect_path_collision(const boost::filesystem::path &p, std::string &error_msg,
                           bool force_overwrite, boost::filesystem::file_type expected_type) {
  const auto path_type = boost::filesystem::status(p).type();
  if (force_overwrite && path_type == boost::filesystem::file_not_found) {
    return true;
  }

  if (expected_type != path_type) {
    switch (path_type) {
      case boost::filesystem::regular_file:
        error_msg +=
            fmt::format(FMT_STRING("Path {} already exists and is actually a file. Please remove "
                                   "the file and try again"),
                        p);
        return false;
      case boost::filesystem::directory_file:
        error_msg += fmt::format(
            FMT_STRING("Path {} already exists and is actually a directory. Please remove "
                       "the directory and try again"),
            p);
        return false;
      default:  // For the time being we only handle regular files and folders
        return true;
    }
  }

  if (force_overwrite) {
    return true;
  }

  if (path_type == boost::filesystem::regular_file) {
    error_msg += fmt::format(FMT_STRING("File {} already exists. Pass --force to overwrite"), p);
    return false;
  }
  return true;
}

template <class EnumT, class StringT>
CliEnumMappings<EnumT, StringT>::CliEnumMappings(const std::initializer_list<value_type> mappings,
                                                 const bool sort_by_key)
    : _mappings(mappings) {
  if (sort_by_key) {
    cppsort::insertion_sort(_mappings.begin(), _mappings.end(), [](const auto &a, const auto &b) {
      return cppsort::natural_less(a.first, b.first);
    });
  }
}

template <class EnumT, class StringT>
CliEnumMappings<EnumT, StringT>::CliEnumMappings(const std::initializer_list<StringT> labels,
                                                 const std::initializer_list<EnumT> enums,
                                                 const bool sort_by_key)
    : CliEnumMappings(ranges::views::zip(labels, enums), sort_by_key) {
  assert(labels.size() == enums.size());
}

template <class EnumT, class StringT>
auto CliEnumMappings<EnumT, StringT>::begin() const -> const_iterator {
  return this->_mappings.begin();
}

template <class EnumT, class StringT>
auto CliEnumMappings<EnumT, StringT>::end() const -> const_iterator {
  return this->_mappings.end();
}

template <class EnumT, class StringT>
auto CliEnumMappings<EnumT, StringT>::cbegin() const -> const_iterator {
  return this->begin();
}

template <class EnumT, class StringT>
auto CliEnumMappings<EnumT, StringT>::cend() const -> const_iterator {
  return this->end();
}

template <class EnumT, class StringT>
auto CliEnumMappings<EnumT, StringT>::find(const EnumT key) const -> const_iterator {
  return std::find_if(this->_mappings.begin(), this->_mappings.end(),
                      [&](const auto &v) { return v.second == key; });
}

template <class EnumT, class StringT>
auto CliEnumMappings<EnumT, StringT>::find(const StringT &key) const -> const_iterator {
  return std::find_if(this->_mappings.begin(), this->_mappings.end(),
                      [&](const auto &v) { return v.first == key; });
}

template <class EnumT, class StringT>
const StringT &CliEnumMappings<EnumT, StringT>::at(const EnumT key) const {
  auto match = this->find(key);
  if (match == this->_mappings.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("Invalid key {}"), int(key)));
  }
  return match->first;
}

template <class EnumT, class StringT>
EnumT CliEnumMappings<EnumT, StringT>::at(const StringT &key) const {
  auto match = this->find(key);
  if (match == this->_mappings.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("Invalid key {}"), key));
  }
  return match->second;
}

template <class EnumT, class StringT>
auto CliEnumMappings<EnumT, StringT>::keys_view() const
    -> decltype(ranges::views::keys(this->_mappings)) {
  return this->_mappings | ranges::views::keys;
}

template <class EnumT, class StringT>
auto CliEnumMappings<EnumT, StringT>::values_view() const
    -> decltype(ranges::views::values(this->_mappings)) {
  auto foo = this->_mappings | ranges::views::values;
  return this->_mappings | ranges::views::values;
}

}  // namespace modle::utils

// IWYU pragma: private, include "modle/utils.hpp"
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <cstdlib>
