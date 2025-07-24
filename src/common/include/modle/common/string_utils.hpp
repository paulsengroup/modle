// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

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

[[nodiscard]] inline std::vector<std::string_view> str_split(std::string_view s, char delim) {
  if (s.empty()) {
    return {};
  }

  std::vector<std::string_view> toks;

  std::size_t i0 = 0;
  std::size_t i1 = s.find(delim, 0);

  while (i1 != std::string_view::npos) {
    if (i1 > i0) {
      toks.emplace_back(s.substr(i0, i1 - i0));
    }
    i0 = i1 + 1;
    i1 = s.find(delim, i0);
  }

  if (i0 < s.size()) {
    // deal with the last token
    toks.emplace_back(s.substr(i0));
  } else if (i0 == s.size() && toks.empty()) {
    // string only contains delimiters
    toks.resize(s.size());
    std::ranges::generate(toks, []() -> std::string_view { return {}; });
  }

  return toks;
}

[[nodiscard]] inline std::vector<std::string_view> str_split(std::string_view s,
                                                             std::string_view delimiters) {
  if (delimiters.empty()) {
    return {s};
  }

  if (s.empty()) {
    return {};
  }

  std::vector<std::string_view> toks;

  std::size_t i0 = 0;
  std::size_t i1 = s.find_first_of(delimiters, 0);

  while (i1 != std::string_view::npos) {
    if (i1 > i0) {
      toks.emplace_back(s.substr(i0, i1 - i0));
    }
    i0 = i1 + 1;
    i1 = s.find_first_of(delimiters, i0);
  }

  if (i0 < s.size()) {
    // deal with the last token
    toks.emplace_back(s.substr(i0));
  } else if (i0 == s.size() && toks.empty()) {
    // string only contains delimiters
    toks.resize(s.size());
    std::ranges::generate(toks, []() -> std::string_view { return {}; });
  }

  return toks;
}

[[nodiscard]] inline std::string str_replace(std::string_view s, std::string_view query,
                                             std::string_view replace) {
  if (s.empty() || query.empty()) {
    return std::string{s};
  }

  std::size_t count = 0;
  auto i = std::string_view{s}.find(query);
  while (i != std::string_view::npos) {
    ++count;
    i = std::string_view{s}.find(query, i + query.size());
  }

  if (count == 0) {
    return std::string{s};
  }

  const auto ssize = static_cast<std::int64_t>(s.size());
  const auto scount = static_cast<std::int64_t>(count);
  const auto ssize_query = static_cast<std::int64_t>(query.size());
  const auto ssize_rep = static_cast<std::int64_t>(replace.size());

  const auto new_size = ssize + (scount * (ssize_rep - ssize_query));
  assert(new_size > 0);
  std::string buff;
  buff.reserve(static_cast<std::size_t>(new_size));

  std::size_t i0 = 0;
  std::size_t i1 = std::string_view{s}.find(query);
  while (i1 != std::string_view::npos) {
    buff.append(s.substr(i0, i1 - i0));
    buff.append(replace);
    i1 = std::string_view{s}.find(query, i0 + query.size());
  }

  if (i0 < s.size()) {
    // deal with the last token
    buff.append(s.substr(i0));
  }

  return buff;
}

}  // namespace modle
