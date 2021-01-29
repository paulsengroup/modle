#pragma once

#include <absl/container/flat_hash_set.h>

#include <cstdint>  // for uint64_t
#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include <utility>  // for move
#include <vector>   // for vector

namespace modle::chr_sizes {

struct ChrSize {
  inline ChrSize(std::string_view chr_name, uint64_t chr_start, uint64_t chr_end);
  inline ChrSize(std::string_view chr_name, uint64_t chr_length);
  inline explicit ChrSize(std::vector<std::string>& toks);

  [[nodiscard]] inline bool operator==(const ChrSize& other) const;
  [[nodiscard]] inline bool operator<(const ChrSize& other) const;

  std::string name;
  uint64_t start;
  uint64_t end;

 private:
  template <typename H>
  [[maybe_unused]] inline friend H AbslHashValue(H h, const ChrSize& c) {
    return H::combine(std::move(h), c.name);
  }
};

class Parser {
 public:
  inline explicit Parser(std::string path_to_chr_sizes);
  inline explicit Parser(std::string_view path_to_chr_sizes);
  [[nodiscard]] inline std::vector<ChrSize> parse_all(char sep = '\t');

 private:
  std::string _path;
  std::ifstream _f{};
  absl::flat_hash_set<ChrSize> _chrs{};
  std::vector<std::string> _errors{};
};
}  // namespace modle::chr_sizes

#include "../../chr_sizes_impl.hpp"