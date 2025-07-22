// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/container/fixed_array.h>

#include <CLI/Formatter.hpp>
#include <filesystem>
#include <initializer_list>
#include <string>
#include <string_view>

#include "modle/common/common.hpp"
#include "modle/common/const_map.hpp"

namespace modle::utils {

// Try to convert str representations like "1.0" or "1.000000" to "1"
[[nodiscard]] inline std::string trim_trailing_zeros_from_decimal_digits(std::string& s);

template <char replacement = '_'>
[[nodiscard]] inline std::string replace_non_alpha_char(std::string& s);

template <class Collection>
[[nodiscard]] inline std::string format_collection_to_english_list(const Collection& collection,
                                                                   std::string_view sep = ", ",
                                                                   std::string_view last_sep = "");

template <class N = double>
[[nodiscard]] inline bool is_finite_number(const std::string& s);

// Returns false in case a collision was detected
inline bool detect_path_collision(
    const std::filesystem::path& p, std::string& error_msg, bool force_overwrite = false,
    std::filesystem::file_type expected_type = std::filesystem::file_type::regular);
[[nodiscard]] inline std::string detect_path_collision(
    const std::filesystem::path& p, bool force_overwrite = false,
    std::filesystem::file_type expected_type = std::filesystem::file_type::regular);

// This class stores pairs of string labels and enums and provides the in a way that is compatible
// with CLI11
template <class EnumT, class StringT = std::string>
class CliEnumMappings {
  static_assert(std::is_convertible_v<StringT, std::string>);

 private:
  absl::FixedArray<std::pair<StringT, EnumT>> _mappings;

 public:
  using key_type = StringT;
  using mapped_type = EnumT;
  using value_type = std::pair<key_type, mapped_type>;
  using reference = typename decltype(_mappings)::const_reference;
  using const_reference = typename decltype(_mappings)::const_reference;
  using pointer = typename decltype(_mappings)::const_pointer;
  using const_pointer = typename decltype(_mappings)::const_pointer;
  using size_type = typename decltype(_mappings)::size_type;
  using iterator = typename decltype(_mappings)::const_iterator;
  using const_iterator = typename decltype(_mappings)::const_iterator;

  inline CliEnumMappings() = default;
  inline CliEnumMappings(std::initializer_list<value_type> mappings, bool sort_by_key = true);
  inline CliEnumMappings(std::initializer_list<StringT> labels, std::initializer_list<EnumT> enums,
                         bool sort_by_key = true);

  [[nodiscard]] inline const_iterator begin() const;
  [[nodiscard]] inline const_iterator end() const;
  [[nodiscard]] inline const_iterator cbegin() const;
  [[nodiscard]] inline const_iterator cend() const;

  [[nodiscard]] inline const_iterator find(EnumT key) const;
  [[nodiscard]] inline const_iterator find(const StringT& key) const;

  [[nodiscard]] inline const StringT& at(EnumT key) const;
  [[nodiscard]] inline EnumT at(const StringT& key) const;

  [[nodiscard]] inline std::vector<std::string> keys() const;
  [[nodiscard]] inline std::vector<std::string> values() const;
};

namespace cli {
class Formatter : public CLI::Formatter {
  [[nodiscard]] inline std::string make_option_opts(const CLI::Option* opt) const override;
};

static constexpr ConstMap<std::string_view, bp_t, 10> genomic_distance_unit_multiplier_map{
    {"bp", bp_t(1)},
    {"k", bp_t(1'000)},
    {"kb", bp_t(1'000)},
    {"kbp", bp_t(1'000)},
    {"m", bp_t(1'000'000)},
    {"mb", bp_t(1'000'000)},
    {"mbp", bp_t(1'000'000)},
    {"g", bp_t(1'000'000'000)},
    {"gb", bp_t(1'000'000'000)},
    {"gbp", bp_t(1'000'000'000)}};

struct IsFiniteValidator : public CLI::Validator {
  inline explicit IsFiniteValidator(bool nan_ok = false);
};

struct TrimTrailingZerosFromDecimalDigitValidator : public CLI::Transformer {
  inline TrimTrailingZerosFromDecimalDigitValidator();
};

struct AsGenomicDistanceTransformer : public CLI::CheckedTransformer {
  inline AsGenomicDistanceTransformer();
};

inline const auto AsGenomicDistance = AsGenomicDistanceTransformer();  // NOLINT(cert-err58-cpp)
inline const auto IsFinite = IsFiniteValidator();                      // NOLINT(cert-err58-cpp)
// NOLINTNEXTLINE(cert-err58-cpp)
inline const auto TrimTrailingZerosFromDecimalDigit = TrimTrailingZerosFromDecimalDigitValidator();

}  // namespace cli

}  // namespace modle::utils

#include "../../../cli_utils_impl.hpp"  // IWYU pragma: export
