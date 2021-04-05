#pragma once

#include <absl/container/flat_hash_set.h>  // for flat_hash_set

#include <cstdint>      // for uint64_t
#include <fstream>      // for ifstream
#include <string>       // for string, basic_string
#include <string_view>  // for string_view
#include <utility>      // for move
#include <vector>       // for vector

namespace modle::chrom_sizes {

struct ChromSize {
  inline ChromSize(std::string_view chrom_name, uint64_t chrom_start, uint64_t chrom_end);
  inline ChromSize(std::string_view chrom_name, uint64_t chrom_length);
  inline explicit ChromSize(std::vector<std::string>& toks);

  [[nodiscard]] inline bool operator==(const ChromSize& other) const;
  [[nodiscard]] inline bool operator<(const ChromSize& other) const;

  std::string name;
  uint64_t start;
  uint64_t end;

 private:
  template <typename H>
  [[maybe_unused]] inline friend H AbslHashValue(H h, const ChromSize& c) {
    return H::combine(std::move(h), c.name);
  }
};

class Parser {
 public:
  inline explicit Parser(std::string path_to_chrom_sizes);
  inline explicit Parser(std::string_view path_to_chrom_sizes);
  [[nodiscard]] inline std::vector<ChromSize> parse_all(char sep = '\t');

 private:
  std::string _path;
  std::ifstream _f{};
  absl::flat_hash_set<ChromSize> _chrs{};
  std::vector<std::string> _errors{};
};
}  // namespace modle::chrom_sizes

#include "../../chrom_sizes_impl.hpp"  // IWYU pragma: keep
