// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/chrom_sizes/chrom_sizes.hpp"  // for Parser

#include <absl/container/flat_hash_set.h>  // for flat_hash_set
#include <absl/strings/str_split.h>        // for StrSplit, Splitter
#include <fmt/compile.h>
#include <fmt/format.h>  // for format, FMT_COMPILE_STRING, FMT_STRING...

#include <cassert>      // for assert
#include <filesystem>   // for filesystem::path
#include <stdexcept>    // for runtime_error, invalid_argument, out_of_range
#include <string>       // for string, operator==
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "modle/bed/bed.hpp"  // for BED

namespace modle::chrom_sizes {

Parser::Parser(const std::filesystem::path& path_to_chrom_sizes) : _reader(path_to_chrom_sizes) {}

std::vector<bed::BED> Parser::parse_all(char sep) {
  std::string buff;
  absl::flat_hash_set<std::string> chrom_names{};
  std::vector<bed::BED> chrom_sizes;

  for (usize i = 1UL, id = 0; this->_reader.getline(buff); ++i) {
    if (buff.empty()) {
      continue;
    }

    const auto splitter = absl::StrSplit(buff, sep);
    const auto num_toks = std::distance(splitter.begin(), splitter.end());
    try {
      if (num_toks < 2) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("expected 2 or more tokens, got {}: \"{}\""), num_toks, buff));
      }
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_NULL_DEREF
      if (const auto chrom_name = *splitter.begin(); chrom_names.contains(chrom_name)) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("found multiple records for chrom \"{}\""), chrom_name));
      }
      DISABLE_WARNING_POP
      chrom_sizes.emplace_back(
          fmt::format(FMT_COMPILE("{}\t0\t{}"), *splitter.begin(), *std::next(splitter.begin())),
          id++, bed::BED::BED3);
    } catch (const std::runtime_error& e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("encountered a malformed record at line {} of file \"{}\": {}.\n "
                                 "Line that triggered the error:\n\"{}\""),
                      i, this->_reader.path_string(), e.what(), buff.data()));
    }
  }
  return chrom_sizes;
}

}  // namespace modle::chrom_sizes
