#include <charconv>
#include <type_traits>

#include "absl/strings/str_format.h"
#include "absl/strings/str_split.h"
#include "modle/utils.hpp"

namespace modle::utils {

template <typename N>
void parse_numeric_or_throw(std::string_view tok, N &field) {
  static_assert(std::is_integral<N>());
  auto [ptr, err] = std::from_chars(tok.data(), tok.end(), field);
  if (ptr != tok.end() && err != std::errc{}) {
    throw_except_from_errc(tok, -1, field, ptr, err);
  }
}

template <typename R>
void parse_real_or_throw(std::string_view tok, R &field) {
  static_assert(std::is_floating_point<R>());
  char *end = nullptr;
  R tmp = std::strtod(tok.data(), &end);
  if (tmp == HUGE_VALF)
    throw_except_from_errc(tok, -1, tmp, nullptr, std::errc::result_out_of_range);
  else if (tmp == 0 && end == tok.data())
    throw_except_from_errc(tok, -1, tmp, nullptr, std::errc::invalid_argument);

  field = tmp;
}

template <typename N>
void parse_numeric_or_throw(const std::vector<std::string_view> &toks, uint8_t idx, N &field) {
  static_assert(std::is_integral<N>());
  auto [ptr, err] = std::from_chars(toks[idx].data(), toks[idx].data() + toks[idx].size(), field);
  if (ptr != toks[idx].end() && err != std::errc{}) {
    throw_except_from_errc(toks[idx], idx, field, ptr, err);
  }
}

template <typename R>
void parse_real_or_throw(const std::vector<std::string_view> &toks, uint8_t idx, R &field) {
  static_assert(std::is_floating_point<R>());
  const std::string tok(toks[idx].begin(), toks[idx].end());
  char *end = nullptr;
  R tmp = std::strtod(tok.data(), &end);
  if (tmp == HUGE_VALF)
    throw_except_from_errc(toks[idx], idx, tmp, nullptr, std::errc::result_out_of_range);
  else if (tmp == 0 && end == tok.data())
    throw_except_from_errc(tok, idx, tmp, nullptr, std::errc::invalid_argument);

  field = tmp;
}

template <typename N>
void parse_vect_of_numbers_or_throw(const std::vector<std::string_view> &toks, uint8_t idx,
                                    std::vector<N> &field, uint64_t expected_size) {
  static_assert(std::is_arithmetic<N>());
  std::vector<std::string_view> ns = absl::StrSplit(toks[idx], ',');
  if (ns.size() != expected_size)
    throw std::runtime_error(
        absl::StrFormat("Expected %lu fields, got %lu.", expected_size, ns.size()));
  field = std::vector<N>(ns.size());
  for (auto i = 0UL; i < expected_size; ++i) parse_numeric_or_throw(ns, i, field[i]);
}

template <typename N>
void throw_except_from_errc(std::string_view tok, int32_t idx, const N &field, const char *c,
                            std::errc e) {
  static_assert(std::is_arithmetic<N>());
  std::string base_error;
  if (idx != -1)
    base_error = absl::StrFormat("Unable to convert field %lu ('%s') to a ", idx, tok);
  else
    base_error = absl::StrFormat("Unable to convert field '%s' to", tok);
  if (std::is_integral<N>()) {
    if (std::is_unsigned<N>())
      base_error += " a positive integral number";
    else
      base_error += "an integral number";
  } else
    base_error += "a real number";
  if (e == std::errc::invalid_argument) {
    if (c != nullptr)
      throw std::runtime_error(
          absl::StrFormat("%s. Reason: found an invalid character '%c'.", base_error, *c));
    throw std::runtime_error(
        absl::StrFormat("%s. Reason: found an invalid character.", base_error));
  }
  if (e == std::errc::result_out_of_range) {
    throw std::runtime_error(absl::StrFormat(
        "%s. Reason: number %s is outside the range of representable numbers [%s, %s].", base_error,
        tok, std::to_string(std::numeric_limits<N>::min()),
        std::to_string(std::numeric_limits<N>::max())));
  }
  throw std::logic_error(
      absl::StrFormat("%s. If you see this error, report it to the developers on "
                      "github.\nBED::throw_except_from_errc "
                      "called with an invalid std::errc. This should not be possible!",
                      base_error));
}
}  // namespace modle::utils
