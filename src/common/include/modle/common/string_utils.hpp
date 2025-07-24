// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cctype>
#include <string>
#include <string_view>

namespace modle {

[[nodiscard]] inline std::string to_lower(std::string_view s) {
  std::string buff(s.size(), '\0');
  std::ranges::transform(s, buff.begin(), [&](char c) { return std::tolower(c); });
  return buff;
}

inline void lower_inplace(std::string& s) {
  std::ranges::transform(s, s.begin(), [&](char c) { return std::tolower(c); });
}

[[nodiscard]] inline std::string_view strip_trailing_whitespace(std::string_view s) {
  const auto it = std::find_if_not(s.rbegin(), s.rend(), [](char c) { return std::isspace(c); });
  return s.substr(0, static_cast<size_t>(s.rend() - it));
}

[[nodiscard]] inline bool str_contains(std::string_view s, std::string_view query) noexcept {
  return s.find(query) != std::string_view::npos;
}

[[nodiscard]] inline std::string_view str_strip_suffix(std::string_view s,
                                                       std::string_view suffix) {
  if (s.ends_with(suffix)) {
    s.remove_suffix(suffix.size());
  }
  return s;
}

}  // namespace modle
