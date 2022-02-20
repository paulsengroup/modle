// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string_view>   // for string_view
#include <system_error>  // for errc
#include <vector>        // for vector

#include "modle/common/common.hpp"  // for usize, u64

namespace modle::utils {

template <class N>
inline auto from_chars(const char* first, const char* last, N& value) noexcept;

template <class N>
inline void parse_numeric_or_throw(std::string_view tok, N& field);

template <class N>
inline N parse_numeric_or_throw(std::string_view tok);

template <class N>
inline void parse_numeric_or_throw(const std::vector<std::string_view>& toks, usize idx, N& field);

template <class N>
inline void parse_vect_of_numbers_or_throw(const std::vector<std::string_view>& toks, usize idx,
                                           std::vector<N>& fields, u64 expected_size);

namespace detail {
template <class N>
inline void throw_except_from_errc(std::string_view tok, usize idx, const N& field, const char* c,
                                   std::errc e);
}
}  // namespace modle::utils

#include "../../../numeric_utils_impl.hpp"  // IWYU pragma: export
// IWYU pragma: no_include <iterator>
// IWYU pragma: no_include <boost/container/detail/std_fwd.hpp>
