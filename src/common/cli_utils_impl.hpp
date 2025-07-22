// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/strings/str_replace.h>
#include <cpp-sort/comparators/natural_less.h>
#include <cpp-sort/sorters/insertion_sorter.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/std.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <filesystem>
#include <string>
#include <string_view>

#include "modle/common/fmt_helpers.hpp"
#include "modle/common/numeric_utils.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"

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
    return fmt::format("{}", fmt::join(collection, last_sep));
  }

  const auto first = std::begin(collection);
  const auto second_to_last = std::end(collection) - 1;
  return fmt::format("{}{}{}", fmt::join(first, second_to_last, sep), last_sep, *second_to_last);
}

std::string trim_trailing_zeros_from_decimal_digits(std::string &s) {
  try {
    auto n = parse_numeric_or_throw<double>(s);
    if (std::trunc(n) == n) {
      return fmt::format("{:.0f}", n);
    }

  } catch ([[maybe_unused]] const std::exception &e) {
    // Let the caller deal with invalid numbers
  }
  return s;
}

template <char replacement>
std::string replace_non_alpha_chars(std::string &s) {
  std::replace_if(s.begin(), s.end(), [](const auto c) { return !std::isalpha(c); }, replacement);
  return s;
}

std::string detect_path_collision(const std::filesystem::path &p, bool force_overwrite,
                                  std::filesystem::file_type expected_type) {
  std::string error_msg;
  detect_path_collision(p, error_msg, force_overwrite, expected_type);
  return error_msg;
}

bool detect_path_collision(const std::filesystem::path &p, std::string &error_msg,
                           bool force_overwrite, std::filesystem::file_type expected_type) {
  using file_type = std::filesystem::file_type;
  const auto path_type = std::filesystem::status(p).type();
  if (force_overwrite && path_type == file_type::not_found) {
    return true;
  }

  if (expected_type != path_type) {
    switch (path_type) {
      case file_type::regular:
        error_msg += fmt::format(
            "path {} already exists and is actually a file. Please remove the file and try again",
            p);
        return false;
      case file_type::directory:
        error_msg += fmt::format(
            "path {} already exists and is actually a directory. Please remove the directory and "
            "try again",
            p);
        return false;
      default:  // For the time being we only handle regular files and folders
        return true;
    }
  }

  if (force_overwrite) {
    return true;
  }

  if (path_type == file_type::regular) {
    error_msg += fmt::format("file {} already exists. Pass --force to overwrite", p);
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

// TODO refactor
namespace internal {
template <class EnumT, class StringT>
[[nodiscard]] inline std::vector<std::pair<StringT, EnumT>> zip_enum_labels(
    const std::initializer_list<StringT> labels, const std::initializer_list<EnumT> enums) {
  assert(labels.size() == enums.size());
  std::vector<std::pair<StringT, EnumT>> result;
  result.reserve(labels.size());
  for (usize i = 0; i < labels.size(); ++i) {
    result.emplace_back(labels[i], enums[i]);
  }

  return result;
}

}  // namespace internal

template <class EnumT, class StringT>
CliEnumMappings<EnumT, StringT>::CliEnumMappings(const std::initializer_list<StringT> labels,
                                                 const std::initializer_list<EnumT> enums,
                                                 const bool sort_by_key)
    : CliEnumMappings(internal::zip_enum_labels(labels, enums), sort_by_key) {
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
    throw std::out_of_range(fmt::format("invalid key {}", int(key)));
  }
  return match->first;
}

template <class EnumT, class StringT>
EnumT CliEnumMappings<EnumT, StringT>::at(const StringT &key) const {
  auto match = this->find(key);
  if (match == this->_mappings.end()) {
    throw std::out_of_range(fmt::format("invalid key {}", key));
  }
  return match->second;
}

template <class EnumT, class StringT>
std::vector<std::string> CliEnumMappings<EnumT, StringT>::keys() const {
  std::vector<std::string> buff{_mappings.size()};
  std::transform(_mappings.begin(), _mappings.end(), buff.begin(),
                 [](const auto &kv) { return std::string{kv.first}; });
  return buff;
}

template <class EnumT, class StringT>
std::vector<std::string> CliEnumMappings<EnumT, StringT>::values() const {
  std::vector<std::string> buff{_mappings.size()};
  std::transform(_mappings.begin(), _mappings.end(), buff.begin(),
                 [](const auto &kv) { return std::string{kv.second}; });
  return buff;
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
        out += fmt::format("={}", s);
      } else {
        out += fmt::format("={}", opt->get_default_str());
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
      out += " (0, inf)";
    } else if (str_contains(t, "NONNEGATIVE") || str_contains(t, "UINT")) {
      out += " [0, inf)";
    }

    if (opt->get_expected_max() == CLI::detail::expected_max_vector_size) {
      out += " ...";
    } else if (opt->get_expected_min() > 1) {
      out += fmt::format(" x {}", opt->get_expected());
    }

    if (opt->get_required()) {
      out += " REQUIRED";
    }
  }
  if (!opt->get_envname().empty()) {
    out += fmt::format(" ({}: {})", get_label("env"), opt->get_envname());
  }
  if (!opt->get_needs().empty()) {
    out += fmt::format(" {}:", get_label("needs"));
    for (const auto *op : opt->get_needs()) {
      out += fmt::format(" {}", op->get_name());
    }
  }
  if (!opt->get_excludes().empty()) {
    out += fmt::format(" {}:", get_label("excludes"));
    for (const auto *op : opt->get_excludes()) {
      out += fmt::format(" {}", op->get_name());
    }
  }

  return out;
}

IsFiniteValidator::IsFiniteValidator(bool nan_ok) {
  description("ISFINITE");

  func_ = [nan_ok](const std::string &input) -> std::string {
    try {
      auto n = parse_numeric_or_throw<double>(input);
      if (std::isfinite(n) || (nan_ok && !std::isnan(n))) {
        return "";
      }
      return fmt::format("value {} is not a finite number", n);
    } catch ([[maybe_unused]] const std::exception &e) {
      return fmt::format("value {} could not be converted", input);
    }
  };
}

TrimTrailingZerosFromDecimalDigitValidator::TrimTrailingZerosFromDecimalDigitValidator()
    : CLI::Transformer({}) {
  description("Trim trailing zeros from the decimal portion of FP numbers");

  func_ = [](std::string &input) {
    std::ignore = trim_trailing_zeros_from_decimal_digits(input);
    return std::string{};
  };
}

AsGenomicDistanceTransformer::AsGenomicDistanceTransformer() : CLI::CheckedTransformer({}) {
  description("Convert common multiple of genomic distances to bp");

  func_ = [](std::string &input) -> std::string {
    const auto &mappings = genomic_distance_unit_multiplier_map;
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_CONVERSION
    std::string_view s{input};
    DISABLE_WARNING_POP
    // Look for the first non-alpha char
    auto unit_rbegin = std::find_if(s.rbegin(), s.rend(),
                                    [&](const auto &c) { return !std::isalpha(c, std::locale()); });

    // input is empty or there are no digits preceding the unit
    if (unit_rbegin == s.rend()) {
      throw CLI::ValidationError(fmt::format("value {} could not be converted", input));
    }

    // input does not have a unit
    if (unit_rbegin == s.rbegin()) {
      try {  // Make sure input can be parsed to a positive integer
        std::ignore = utils::parse_numeric_or_throw<bp_t>(input);
      } catch ([[maybe_unused]] const std::exception &e) {
        throw CLI::ValidationError(fmt::format("unable to convert {} to a number", input));
      }
      return "";  // input validation was successful
    }

    const auto i = s.size() - static_cast<usize>(std::distance(s.rbegin(), unit_rbegin));
    auto num_str = input.substr(0, i);
    auto unit = input.substr(i);

    // Look-up the appropriate multiplier
    auto it = mappings.find(absl::AsciiStrToLower(unit));
    if (it == mappings.end()) {
      throw CLI::ValidationError(
          fmt::format("{} unit not recognized.\n"
                      "Valid units:\n - {}",
                      unit, fmt::join(mappings.keys(), "\n - ")));
    }

    const auto &multiplier = *it.second;
    try {  // Make sure input can be parsed to a positive integer
      const auto m =
          static_cast<double>(multiplier) * utils::parse_numeric_or_throw<double>(num_str);

      if (std::trunc(m) != m) {  // Ensure double can be represented as a whole number
        throw CLI::ValidationError(fmt::format(
            "Unable to convert {} to a number of base-pairs ({} is not an integral number)", s, m));
      }
      input = fmt::format("{:.0f}", m);
      return "";  // input validation and transformation were successful
    } catch (const CLI::ValidationError &e) {
      throw;
    } catch ([[maybe_unused]] const std::exception &e) {
      throw CLI::ValidationError(fmt::format("unable to convert {} to a number", num_str));
    }
  };
}

}  // namespace cli

}  // namespace modle::utils

// IWYU pragma: private, include "modle/utils.hpp"
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include <cstdlib>
