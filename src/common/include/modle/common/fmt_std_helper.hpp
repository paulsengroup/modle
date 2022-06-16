// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <filesystem>

template <>
struct fmt::formatter<std::filesystem::path> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw fmt::format_error("invalid format");
    }
    return ctx.end();
  }

  template <typename FormatContext>
  auto format(const std::filesystem::path& p, FormatContext& ctx) const -> decltype(ctx.out()) {
    return fmt::format_to(ctx.out(), FMT_STRING("\"{}\""), p.string());
  }
};
