#pragma once
#include <boost/exception/all.hpp>
#include <boost/stacktrace.hpp>
#include <cstdint>
#include <string_view>
#include <system_error>
#include <vector>

namespace modle::utils {

template <typename N>
inline void parse_numeric_or_throw(std::string_view tok, N& field);

template <typename I>
inline void parse_int_or_throw(std::string_view tok, I& field);
// This is a temporary workaround to deal with the fact that libstdc++ 10 does not come with the
// float/double overloads for std::from_chars
template <typename R>
inline void parse_real_or_throw(std::string_view tok, R& field);

template <typename N>
inline void parse_numeric_or_throw(const std::vector<std::string_view>& toks, std::size_t idx,
                                   N& field);
template <typename I>
inline void parse_int_or_throw(const std::vector<std::string_view>& toks, std::size_t idx,
                               I& field);
// This is a temporary workaround to deal with the fact that libstdc++ 10 does not come with the
// float/double overloads for std::from_chars
template <typename R>
inline void parse_real_or_throw(const std::vector<std::string_view>& toks, std::size_t idx,
                                R& field);

// Because of the issue outlined for parse_real_or_throw(), this currently only works for integral
// numbers (which should be fine for our purposes)
template <typename N>
inline void parse_vect_of_numbers_or_throw(const std::vector<std::string_view>& toks,
                                           std::size_t idx, std::vector<N>& field,
                                           uint64_t expected_size);

template <typename N>
inline void throw_except_from_errc(std::string_view tok, std::size_t idx, const N& field,
                                   const char* c, std::errc e);

[[nodiscard]] inline std::string init_juicer_tools_argv(std::string_view path_to_juicer_tools,
                                                        uint64_t juicer_tools_mem);

[[nodiscard]] inline bool chr_equal_operator(std::string_view chr1, std::string_view chr2);
[[nodiscard]] inline bool chr_equal_operator(const std::pair<std::string_view, int64_t>& chr1,
                                             const std::pair<std::string_view, int64_t>& chr2);

[[nodiscard]] inline bool chr_less_than_operator(std::string_view chr1, std::string_view chr2);
[[nodiscard]] inline bool chr_less_than_operator(const std::pair<std::string_view, int64_t>& chr1,
                                                 const std::pair<std::string_view, int64_t>& chr2);

template <class E>
[[noreturn]] inline void throw_with_trace(const E& e);

// Typetraits stuff

// Usage: static_assert(std::experimental::is_detected_convertible_v<const char*,
// utils::has_data_member_func, T>, "error msg");
// This static_assert will fail for types that don't have a member function .data() that returns a
// const char*

template <class T>
using has_data_member_func = decltype(std::declval<T&>().data());

template <class T>
using has_size_member_func = decltype(std::declval<T&>().size());

template <typename T>
[[nodiscard]] inline std::string get_printable_type_name(const T& var);

// https://stackoverflow.com/a/56766138
template <typename T>
constexpr auto get_printable_type_name() noexcept;

}  // namespace modle::utils

#include "../../utils_impl.hpp"
