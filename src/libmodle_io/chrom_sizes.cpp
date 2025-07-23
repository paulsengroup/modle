// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/chrom_sizes/chrom_sizes.hpp"

#include <absl/strings/ascii.h>
#include <absl/strings/str_split.h>
#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <parallel_hashmap/phmap.h>

#include <cassert>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "modle/bed/bed.hpp"
#include "modle/common/utils.hpp"

namespace modle::chrom_sizes {

Parser::Parser(const std::filesystem::path& path_to_chrom_sizes) : _reader(path_to_chrom_sizes) {}

std::vector<bed::BED> Parser::parse_all(char sep) {
  std::string buff;
  phmap::flat_hash_set<std::string> chrom_names{};
  std::vector<bed::BED> chrom_sizes;

  for (usize i = 1UL, id = 0; _reader.getline(buff); ++i) {
    buff = absl::StripTrailingAsciiWhitespace(buff);
    if (buff.empty()) {
      continue;
    }

    const auto splitter = absl::StrSplit(buff, sep);
    const auto num_toks = std::distance(splitter.begin(), splitter.end());
    try {
      if (num_toks != 2) {
        throw std::runtime_error(
            fmt::format("expected exactly 2 fields, found {}: \"{}\"", num_toks, buff));
      }
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_NULL_DEREF
      const auto chrom_name = utils::strip_quote_pairs(*splitter.begin());
      const auto chrom_size = *std::next(splitter.begin());
      DISABLE_WARNING_POP
      if (chrom_names.contains(chrom_name)) {
        throw std::runtime_error(
            fmt::format("found multiple records for chrom \"{}\"", chrom_name));
      }

      if (chrom_size == "0") {
        throw std::runtime_error(fmt::format("chrom \"{}\" has a length of 0bp", chrom_name));
      }
      chrom_sizes.emplace_back(
          fmt::format(FMT_COMPILE("{}\t0\t{}"), chrom_name, *std::next(splitter.begin())), id++,
          bed::BED::BED3);
    } catch (const std::runtime_error& e) {
      throw std::runtime_error(
          fmt::format("encountered a malformed record at line {} of file {}: {}.\n "
                      "Line that triggered the error:\n\"{}\"",
                      i, _reader.path(), e.what(), buff.data()));
    }
  }
  return chrom_sizes;
}

}  // namespace modle::chrom_sizes
