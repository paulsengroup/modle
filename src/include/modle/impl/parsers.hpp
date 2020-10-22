#pragma once

#include <charconv>

#include "absl/strings/str_format.h"
#include "absl/strings/str_split.h"

namespace modle {
template <typename N, typename>
void BED::parse_numeric_or_throw(const std::vector<std::string_view> &toks, uint8_t idx, N &field) {
  auto [ptr, err] = std::from_chars(toks[idx].data(), toks[idx].data() + toks[idx].size(), field);
  if (ptr != toks[idx].data() + toks.size() && err != std::errc{}) {
    throw_except_from_errc(toks[idx], idx, err);
  }
}

template <typename R, typename>
void BED::parse_real_or_throw(const std::vector<std::string_view> &toks, uint8_t idx, R &field) {
  const std::string tok(toks[idx].begin(), toks[idx].end());
  char *end;
  field = std::strtod(tok.data(), &end);
  if (auto ptr = tok.data() + tok.size() ; ptr != end) {
    throw std::runtime_error(absl::StrFormat(
        "Failed to convert field #%lu to a floating-point number. Field #4: '%s'", idx, toks[idx]));
  }
}

template <typename N, typename>
void BED::parse_vect_of_numbers_or_throw(const std::vector<std::string_view> &toks, uint8_t idx,
                                         std::vector<N> &field, uint64_t expected_size) {
  std::vector<std::string_view> ns = absl::StrSplit(toks[idx], ',');
  if (ns.size() != expected_size)
    throw std::runtime_error(
        absl::StrFormat("Expected %lu fields, got %lu.", expected_size, ns.size()));
  field = std::vector<N>(ns.size());
  for (auto i = 0UL; i < expected_size; ++i) parse_numeric_or_throw(ns, i, field[i]);
}

}  // namespace modle