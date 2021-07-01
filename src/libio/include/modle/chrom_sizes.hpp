#pragma once

#include <fstream>  // for ifstream
#include <vector>   // for vector

#include <boost/filesystem/path.hpp>  // for filesystem::path

namespace modle::bed {
struct BED;
}

namespace modle::chrom_sizes {

class Parser {
 public:
  explicit Parser(boost::filesystem::path path_to_chrom_sizes);
  [[nodiscard]] std::vector<bed::BED> parse_all(char sep = '\t');

 private:
  boost::filesystem::path _path;
  std::ifstream _fp{};
};
}  // namespace modle::chrom_sizes
