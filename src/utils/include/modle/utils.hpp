#pragma once
#include <cstdint>
#include <string_view>
#include <system_error>
#include <vector>

namespace modle::utils {

template <typename N>
static void parse_numeric_or_throw(std::string_view tok, N& field);
// This is a temporary workaround to deal with the fact that libstdc++ 10 does not come with the
// float/double overloads for std::from_chars
template <typename R>
static void parse_real_or_throw(std::string_view tok, R& field);

template <typename N>
static void parse_numeric_or_throw(const std::vector<std::string_view>& toks, std::size_t idx,
                                   N& field);
// This is a temporary workaround to deal with the fact that libstdc++ 10 does not come with the
// float/double overloads for std::from_chars
template <typename R>
static void parse_real_or_throw(const std::vector<std::string_view>& toks, std::size_t idx,
                                R& field);

// Because of the issue outlined for parse_real_or_throw(), this currently only works for integral
// numbers (which should be fine for our purposes)
template <typename N>
static void parse_vect_of_numbers_or_throw(const std::vector<std::string_view>& toks,
                                           std::size_t idx, std::vector<N>& field,
                                           uint64_t expected_size);

template <typename N>
static void throw_except_from_errc(std::string_view tok, std::size_t idx, const N& field,
                                   const char* c, std::errc e);

}  // namespace modle::utils

#include "../../utils_impl.hpp"