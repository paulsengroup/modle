// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <string_view>

#include "modle/common/common.hpp"

namespace modle::config::version {

[[nodiscard]] u8f major() noexcept;
[[nodiscard]] u8f minor() noexcept;
[[nodiscard]] u8f patch() noexcept;

[[nodiscard]] std::string_view suffix() noexcept;

[[nodiscard]] std::string_view str() noexcept;
[[nodiscard]] std::string_view str_long(std::string_view prefix = "MoDLE") noexcept;

}  // namespace modle::config::version
