#include "modle/chr_sizes.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "fmt/format.h"

namespace modle::chr_sizes {
Parser::Parser(std::string path_to_chr_sizes) : _path(std::move(path_to_chr_sizes)) {}
Parser::Parser(std::string_view path_to_chr_sizes) : _path(std::string(path_to_chr_sizes)) {}

std::vector<ChrSize> Parser::parse(char sep) {
  std::string buff;
  std::vector<std::string> tokens;
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
        throw std::runtime_error(
            fmt::format("Expected 2 or more tokens, got {}: '{}'", tokens.size(), buff));
      }
      ChrSize record(tokens);
      if (this->_chrs.contains(record)) {
        throw std::runtime_error(fmt::format("Found multiple records for chr '{}'", record.name));
      }
      this->_chrs.emplace(std::move(record));

    } catch (const std::runtime_error& e) {
      this->_errors.push_back(
          fmt::format("Encountered a malformed record at line {} of file '{}': {}.\n "
                      "Line that triggered the error:\n'{}'",
                      i, this->_path, e.what(), buff.data()));
    }
  }
  if (!this->_f && !this->_f.eof()) {
    throw fmt::system_error(errno, "IO error while reading file '{}'", this->_path);
  }
  if (!this->_errors.empty()) {
    throw std::runtime_error(
        fmt::format("The following error(s) occurred while parsing file '{}':\n - {}", this->_path,
                    absl::StrJoin(this->_errors, "\n - ")));
  }
  return {this->_chrs.begin(), this->_chrs.end()};
}

ChrSize::ChrSize(std::vector<std::string>& toks) {
  assert(toks.size() > 1);
  this->name = std::move(toks[0]);
  std::size_t idx = 0;  // Used to print useful errors
  try {
    if (toks.size() == 2) {
      this->start = 0ULL;
      this->end = std::stoull(toks[++idx]);
    } else {
      this->start = std::stoull(toks[++idx]);
      this->end = std::stoull(toks[++idx]);
      if (this->start > this->end) {
        throw std::runtime_error(
            fmt::format("Sequence '{}' has a starting coord. that is larger than its end coord.: "
                        "start={} > end={}",
                        this->name, this->start, this->end));
      }
    }

  } catch (const std::invalid_argument& e) {
    if (toks.size() == 2) {
      throw std::runtime_error(fmt::format("Sequence '{}' has an invalid length of '{}': {}",
                                           this->name, toks[idx], e.what()));
    }
    throw std::runtime_error(fmt::format("Sequence '{}' has an invalid {} coord. of '{}': {}",
                                         this->name, idx == 1 ? "start" : "end", toks[idx],
                                         e.what()));

  } catch (const std::out_of_range& e) {
    if (toks.size() == 2) {
      throw std::runtime_error(fmt::format("Sequence '{}' has an out of bound length of '{}': {}",
                                           this->name, toks[idx], e.what()));
    }
    throw std::runtime_error(fmt::format("Sequence '{}' has an out of bound {} coord. of '{}': {}",
                                         this->name, idx == 1 ? "start" : "end", toks[idx],
                                         e.what()));
  }
}

ChrSize::ChrSize(std::string_view chr_name, uint64_t length)
    : name(chr_name), start(0), end(length) {}
ChrSize::ChrSize(std::string_view chr_name, uint64_t start, uint64_t end)
    : name(chr_name), start(start), end(end) {}

bool ChrSize::operator==(const ChrSize& other) const {
  return this->name == other.name && this->start == other.start && this->end == other.end;
}

bool ChrSize::operator<(const ChrSize& other) const {
  return this->end - this->start < other.end - other.start;
}
}  // namespace modle::chr_sizes