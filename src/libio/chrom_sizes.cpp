#include "modle/chrom_sizes.hpp"

#include <absl/container/flat_hash_set.h>  // for flat_hash_set, BitMask
#include <absl/strings/str_join.h>         // for StrJoin
#include <absl/strings/str_split.h>        // for StrSplit, Splitter
#include <fmt/format.h>                    // for format, FMT_COMPILE_STRING, FMT_STRING...

#include <algorithm>    // for copy, max
#include <cassert>      // for assert
#include <cerrno>       // for errno
#include <cstddef>      // for size_t
#include <cstdint>      // for uint64_t
#include <fstream>      // for ifstream, basic_ios
#include <stdexcept>    // for runtime_error, invalid_argument, out_of_range
#include <string>       // for string, basic_string, stoull, getline, operator==
#include <string_view>  // for string_view
#include <utility>      // for move
#include <vector>       // for vector

namespace modle::chrom_sizes {

ChromSize::ChromSize(std::string_view chrom_name, uint64_t chrom_start, uint64_t chrom_end)
    : name(chrom_name), start(chrom_start), end(chrom_end) {}

ChromSize::ChromSize(std::string_view chrom_name, uint64_t chrom_length)
    : ChromSize(chrom_name, 0UL, chrom_length) {}

ChromSize::ChromSize(std::vector<std::string>& toks) {
  assert(toks.size() > 1);
  this->name = std::move(toks[0]);
  size_t idx = 0;  // Used to print useful errors
  try {
    if (toks.size() == 2) {
      this->start = 0ULL;
      this->end = std::stoull(toks[++idx]);
    } else {
      this->start = std::stoull(toks[++idx]);
      this->end = std::stoull(toks[++idx]);
      if (this->start > this->end) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("Sequence '{}' has a starting coord. that is larger than its end coord.: "
                       "start={} > end={}"),
            this->name, this->start, this->end));
      }
    }

  } catch (const std::invalid_argument& e) {
    if (toks.size() == 2) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Sequence '{}' has an invalid simulated_length of '{}': {}"),
                      this->name, toks[idx], e.what()));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Sequence '{}' has an invalid {} coord. of '{}': {}"), this->name,
                    idx == 1 ? "start" : "end", toks[idx], e.what()));

  } catch (const std::out_of_range& e) {
    if (toks.size() == 2) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Sequence '{}' has an out of bound simulated_length of '{}': {}"),
                      this->name, toks[idx], e.what()));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Sequence '{}' has an out of bound {} coord. of '{}': {}"),
                    this->name, idx == 1 ? "start" : "end", toks[idx], e.what()));
  }
}

bool ChromSize::operator==(const ChromSize& other) const {
  return this->name == other.name && this->start == other.start && this->end == other.end;
}

bool ChromSize::operator<(const ChromSize& other) const {
  return this->end - this->start < other.end - other.start;
}

Parser::Parser(std::string path_to_chrom_sizes) : _path(std::move(path_to_chrom_sizes)) {}
Parser::Parser(std::string_view path_to_chrom_sizes) : Parser(std::string{path_to_chrom_sizes}) {}

std::vector<ChromSize> Parser::parse_all(char sep) {
  std::string buff;
  std::vector<std::string> tokens;
  std::vector<ChromSize> chrom_sizes;
  this->_f = std::ifstream(this->_path);
  for (auto i = 1UL; this->_f.good(); ++i) {
    if (std::getline(this->_f, buff); buff.empty()) {
      continue;
    }

    if (!this->_f) {
      throw fmt::system_error(errno, "IO error while reading file '{}'", this->_path);
    }
    tokens = absl::StrSplit(buff, sep);
    assert(!tokens.empty());  // This should only happen at EOF, which is handled elsewhere
    try {
      if (tokens.size() < 2) {
        throw std::runtime_error(fmt::format(FMT_STRING("Expected 2 or more tokens, got {}: '{}'"),
                                             tokens.size(), buff));
      }
      ChromSize record(tokens);
      if (this->_chrs.contains(record)) {
        throw std::runtime_error(fmt::format("Found multiple records for chrom '{}'", record.name));
      }
      this->_chrs.emplace(record);
      chrom_sizes.emplace_back(std::move(record));

    } catch (const std::runtime_error& e) {
      this->_errors.push_back(
          fmt::format(FMT_STRING("Encountered a malformed record at line {} of file '{}': {}.\n "
                                 "Line that triggered the error:\n'{}'"),
                      i, this->_path, e.what(), buff.data()));
    }
  }
  if (!this->_f && !this->_f.eof()) {
    throw fmt::system_error(errno, "IO error while reading file '{}'", this->_path);
  }
  if (!this->_errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error(s) occurred while parsing file '{}':\n - {}"),
                    this->_path, absl::StrJoin(this->_errors, "\n - ")));
  }
  return chrom_sizes;
}

}  // namespace modle::chrom_sizes
