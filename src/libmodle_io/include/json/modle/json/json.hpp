// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "modle/toml/toml.hpp"

namespace modle {

[[nodiscard]] inline nlohmann::json reformat_json_nulls(nlohmann::json attributes) {
  std::vector<std::string> null_fields{};
  for (const auto& field : attributes.items()) {
    if (field.value() == "null") {
      null_fields.emplace_back(field.key());
    }
  }

  for (const auto& k : null_fields) {
    attributes[k] = nullptr;
  }

  return attributes;
}

[[nodiscard]] inline nlohmann::json toml_to_json(const toml::table& t) {
  std::stringstream buff;

  // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
  buff << toml::json_formatter(t);

  auto j = reformat_json_nulls(nlohmann::json::parse(buff.str()));
  if (const auto metadata = t.find("metadata");
      metadata != t.end() && metadata->second.is_string()) {
    try {
      j["metadata"] =
          reformat_json_nulls(nlohmann::json::parse(metadata->second.ref<std::string>()));
      // NOLINTNEXTLINE
    } catch (...) {
    }
  }

  return j;
}

}  // namespace modle
