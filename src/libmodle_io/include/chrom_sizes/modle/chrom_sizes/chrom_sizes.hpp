// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <boost/filesystem/path.hpp>  // for filesystem::path
#include <vector>                     // for vector

#include "modle/compressed_io/compressed_io.hpp"  // for Reader

namespace modle::bed {
struct BED;
}

namespace modle::chrom_sizes {

class Parser {
 public:
  explicit Parser(const boost::filesystem::path& path_to_chrom_sizes);
  [[nodiscard]] std::vector<bed::BED> parse_all(char sep = '\t');

 private:
  compressed_io::Reader _reader{};
};
}  // namespace modle::chrom_sizes
