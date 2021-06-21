#include "modle/chrom_sizes.hpp"  // for Parser

#include <absl/container/flat_hash_set.h>  // for flat_hash_set
#include <absl/strings/str_join.h>         // for StrJoin
#include <absl/strings/str_split.h>        // for StrSplit, Splitter
#include <fmt/format.h>                    // for format, FMT_COMPILE_STRING, FMT_STRING...
#include <fmt/ostream.h>

#include <cassert>      // for assert
#include <cerrno>       // for errno
#include <cstddef>      // for size_t
#include <filesystem>   // for filesystem::path
#include <fstream>      // for ifstream, basic_ios
#include <stdexcept>    // for runtime_error, invalid_argument, out_of_range
#include <string>       // for string, basic_string, getline, operator==
#include <string_view>  // for string_view
#include <utility>      // for move
#include <vector>       // for vector

#include "modle/bed.hpp"  // for BED

namespace modle::chrom_sizes {

Parser::Parser(std::filesystem::path path_to_chrom_sizes) : _path(std::move(path_to_chrom_sizes)) {
  this->_fp.open(this->_path);
  if (!this->_fp) {
    throw fmt::system_error(errno, "Unable to open file {} for reading", this->_path);
  }
}

std::vector<bed::BED> Parser::parse_all(char sep) {
  std::string buff;
  std::vector<std::string_view> tokens;
  absl::flat_hash_set<std::string> chrom_names{};
  std::vector<bed::BED> chrom_sizes;

  for (auto i = 1UL, id = 0UL; this->_fp.good(); ++i) {
    if (std::getline(this->_fp, buff); buff.empty()) {
      continue;
    }

    if (!this->_fp) {
      throw fmt::system_error(errno, "IO error while reading file {}", this->_path);
    }
    tokens = absl::StrSplit(buff, sep);
    assert(!tokens.empty());  // NOLINT; This should only happen at EOF, which is handled elsewhere
    try {
      if (tokens.size() < 2) {
        throw std::runtime_error(fmt::format(FMT_STRING("Expected 2 or more tokens, got {}: '{}'"),
                                             tokens.size(), buff));
      }
      tokens.insert(tokens.begin() + 1, "0");
      if (const auto& chrom_name = tokens.front(); chrom_names.contains(chrom_name)) {
        throw std::runtime_error(fmt::format("Found multiple records for chrom '{}'", chrom_name));
      }
      chrom_sizes.emplace_back(absl::StrJoin(tokens, "\t"), id++, bed::BED::BED3);
    } catch (const std::runtime_error& e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Encountered a malformed record at line {} of file {}: {}.\n "
                                 "Line that triggered the error:\n'{}'"),
                      i, this->_path, e.what(), buff.data()));
    }
  }
  if (!this->_fp && !this->_fp.eof()) {
    throw fmt::system_error(errno, "IO error while reading file {}", this->_path);
  }
  return chrom_sizes;
}

}  // namespace modle::chrom_sizes
