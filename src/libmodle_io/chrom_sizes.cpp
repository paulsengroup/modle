// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/chrom_sizes/chrom_sizes.hpp"  // for Parser

#include <absl/container/flat_hash_set.h>  // for flat_hash_set
#include <absl/strings/str_join.h>         // for StrJoin
#include <absl/strings/str_split.h>        // for StrSplit, Splitter
#include <fmt/format.h>                    // for format, FMT_COMPILE_STRING, FMT_STRING...

#include <boost/filesystem/path.hpp>  // for filesystem::path
#include <cassert>                    // for assert
#include <stdexcept>                  // for runtime_error, invalid_argument, out_of_range
#include <string>                     // for string, operator==
#include <string_view>                // for string_view
#include <vector>                     // for vector

#include "modle/bed/bed.hpp"  // for BED

namespace modle::chrom_sizes {

Parser::Parser(const boost::filesystem::path& path_to_chrom_sizes) : _reader(path_to_chrom_sizes) {}

std::vector<bed::BED> Parser::parse_all(char sep) {
  std::string buff;
  std::vector<std::string_view> tokens;
  absl::flat_hash_set<std::string> chrom_names{};
  std::vector<bed::BED> chrom_sizes;

  for (usize i = 1UL, id = 0; this->_reader.getline(buff); ++i) {
    if (buff.empty()) {
      continue;
    }

    tokens = absl::StrSplit(buff, sep);
    assert(!tokens.empty());
    try {
      if (tokens.size() < 2) {
        throw std::runtime_error(fmt::format(FMT_STRING("Expected 2 or more tokens, got {}: \"{}\""),
                                             tokens.size(), buff));
      }
      tokens.insert(tokens.begin() + 1, "0");
      if (const auto& chrom_name = tokens.front(); chrom_names.contains(chrom_name)) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Found multiple records for chrom \"{}\""), chrom_name));
      }
      chrom_sizes.emplace_back(absl::StrJoin(tokens, "\t"), id++, bed::BED::BED3);
    } catch (const std::runtime_error& e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Encountered a malformed record at line {} of file \"{}\": {}.\n "
                                 "Line that triggered the error:\n\"{}\""),
                      i, this->_reader.path_string(), e.what(), buff.data()));
    }
  }
  return chrom_sizes;
}

}  // namespace modle::chrom_sizes
