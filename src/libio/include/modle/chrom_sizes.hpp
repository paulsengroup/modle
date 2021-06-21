#pragma once

#include <filesystem>  // for filesystem::path
#include <fstream>     // for ifstream
#include <vector>      // for vector

namespace modle::bed {
struct BED;
}

namespace modle::chrom_sizes {

class Parser {
 public:
  explicit Parser(std::filesystem::path path_to_chrom_sizes);
  [[nodiscard]] std::vector<bed::BED> parse_all(char sep = '\t');

 private:
  std::filesystem::path _path;
  std::ifstream _fp{};
};
}  // namespace modle::chrom_sizes
