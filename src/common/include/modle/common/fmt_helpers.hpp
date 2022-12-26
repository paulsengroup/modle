// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include "modle/common/pixel.hpp"

template <>
struct fmt::formatter<modle::PixelCoordinates> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw fmt::format_error("invalid format");
    }
    return ctx.end();
  }

  template <class FormatContext>
  auto format(const modle::PixelCoordinates& c, FormatContext& ctx) const -> decltype(ctx.out()) {
    return fmt::format_to(ctx.out(), FMT_STRING("[{}, {}]"), c.row, c.col);
  }
};

template <class N>
struct fmt::formatter<modle::Pixel<N>> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw fmt::format_error("invalid format");
    }
    return ctx.end();
  }

  template <class FormatContext>
  auto format(const modle::Pixel<N>& p, FormatContext& ctx) const -> decltype(ctx.out()) {
    return fmt::format_to(ctx.out(), FMT_STRING("{}: {}"), p.coords, p.count);
  }
};
