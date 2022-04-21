// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/strings/str_replace.h>
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

std::string trim_trailing_zeros_from_decimal_digits(std::string &&s) {
  try {
    auto n = parse_numeric_or_throw<double>(s);
    if (std::trunc(n) == n) {
      return fmt::format(FMT_STRING("{:.0f}"), n);
    }

  } catch ([[maybe_unused]] const std::exception &e) {
    // Let the caller deal with invalid numbers
  }
  return s;
}

template <char replacement>
std::string replace_non_alpha_chars(std::string &&s) {
  std::replace_if(
      s.begin(), s.end(), [](const auto c) { return !std::isalpha(c); }, replacement);
  return s;
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

namespace cli {
std::string Formatter::make_option_opts(const CLI::Option *opt) const {
  if (!opt->get_option_text().empty()) {
    return opt->get_option_text();
  }

  auto str_contains = [](const auto s, const auto query) {
    return s.find(query) != decltype(s)::npos;
  };

  std::string out;
  if (opt->get_type_size() != 0) {
    // Format default values so that the help string reads like: --my-option=17.0
    if (!opt->get_default_str().empty()) {
      if (absl::StartsWith(opt->get_type_name(), "FLOAT")) {
        auto s = opt->get_default_str();
        if (s.find('.') == std::string::npos) {
          s += ".0";
        }
        out += fmt::format(FMT_STRING("={}"), s);
      } else {
        out += fmt::format(FMT_STRING("={}"), opt->get_default_str());
      }
    }

    // Format param domain using open/closed interval notation
    if (const auto &t = opt->get_type_name(); str_contains(t, " in ")) {
      const auto p1 = t.find("[", t.find(" in "));
      const auto p2 = t.find("]", t.find(" in "));
      if (p1 != std::string::npos && p2 != std::string::npos && p2 > p1) {
        out += " " + absl::StrReplaceAll(t.substr(p1, p2), {{" - ", ", "}});
      }
    } else if (str_contains(t, "POSITIVE")) {
      out += " (0, inf]";
    } else if (str_contains(t, "NONNEGATIVE") || str_contains(t, "UINT")) {
      out += " [0, inf]";
    }

    if (opt->get_expected_max() == CLI::detail::expected_max_vector_size) {
      out += " ...";
    } else if (opt->get_expected_min() > 1) {
      out += fmt::format(FMT_STRING(" x {}"), opt->get_expected());
    }

    if (opt->get_required()) {
      out += " REQUIRED";
    }
  }
  if (!opt->get_envname().empty()) {
    out += fmt::format(FMT_STRING(" ({}: {})"), get_label("env"), opt->get_envname());
  }
  if (!opt->get_needs().empty()) {
    out += fmt::format(FMT_STRING(" {}:"), get_label("needs"));
    for (const auto *op : opt->get_needs()) {
      out += fmt::format(FMT_STRING(" {}"), op->get_name());
    }
  }
  if (!opt->get_excludes().empty()) {
    out += fmt::format(FMT_STRING(" {}:"), get_label("excludes"));
    for (const auto *op : opt->get_excludes()) {
      out += fmt::format(FMT_STRING(" {}"), op->get_name());
    }
  }

  return out;
}

IsFinite::IsFinite(bool nan_ok) {
  description("ISFINITE");

  func_ = [nan_ok](const std::string &input) -> std::string {
    try {
      auto n = parse_numeric_or_throw<double>(input);
      if (std::isfinite(n) || (nan_ok && !std::isnan(n))) {
        return "";
      }
      return fmt::format(FMT_STRING("Value {} is not a finite number"), n);
    } catch ([[maybe_unused]] const std::exception &e) {
      return fmt::format(FMT_STRING("Value {} could not be converted"), input);
    }
  };
}

}  // namespace cli

}  // namespace modle::utils

// IWYU pragma: private, include "modle/utils.hpp"
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <cstdlib>
