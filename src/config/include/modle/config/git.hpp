// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string_view>

namespace modle::config::git {

[[nodiscard]] bool state_available() noexcept;
[[nodiscard]] std::string_view head_sha1() noexcept;
[[nodiscard]] bool is_dirty() noexcept;
[[nodiscard]] std::string_view author_name() noexcept;
[[nodiscard]] std::string_view author_email() noexcept;
[[nodiscard]] std::string_view commit_date() noexcept;
[[nodiscard]] std::string_view commit_subject() noexcept;
[[nodiscard]] std::string_view commit_body() noexcept;
[[nodiscard]] std::string_view describe() noexcept;
[[nodiscard]] std::string_view branch() noexcept;
[[nodiscard]] std::string_view tag() noexcept;

}  // namespace modle::config::git
