#pragma once
#include <cstdint>
#include <fstream>
#include <string_view>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"

namespace modle::chr_sizes {

struct ChrSize {
  ChrSize(std::string_view chr_name, uint64_t length);
  ChrSize(std::string_view chr_name, uint64_t start, uint64_t end);
  explicit ChrSize(std::vector<std::string>& toks);

  [[nodiscard]] bool operator==(const ChrSize& other) const;
  [[nodiscard]] bool operator<(const ChrSize& other) const;

  std::string name;
  uint64_t start;
  uint64_t end;

  template <typename H>
  friend H AbslHashValue(H h, const ChrSize& c) {
    return H::combine(std::move(h), c.name);
  }
};

class Parser {
 public:
  explicit Parser(std::string path_to_chr_sizes);
  explicit Parser(std::string_view path_to_chr_sizes);
  std::vector<ChrSize> parse(char sep = '\t');

 private:
  std::string _path;
  std::ifstream _f{};
  absl::flat_hash_set<ChrSize> _chrs{};
  std::vector<std::string> _errors{};
};
}  // namespace modle::chr_sizes