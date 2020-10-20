#pragma once
#include <array>
#include <fstream>
#include <string>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"

namespace modle {
struct SimpleBED {
  [[nodiscard]] bool empty() const;
  [[nodiscard]] std::string to_string() const;
  enum strand { plus = '+', minus = '-', unknown = '.' };

  std::string chr;
  uint32_t chr_start;
  uint32_t chr_end;
  double score{0};
  strand strnd;
};

class SimpleBEDParser {
 public:
  explicit SimpleBEDParser(std::string_view path_to_bed);
  std::vector<SimpleBED> parse_all();
  SimpleBED parse_next();

 private:
  std::string _path_to_bed;
  std::ifstream _bed;
  std::array<char, 2048> _buff{0};
  uint32_t _line{1};
};

struct ChrSize {
  ChrSize(std::string_view chr_name, uint64_t length);
  explicit ChrSize(std::vector<std::string>& toks);

  [[nodiscard]] bool operator==(const ChrSize&other) const;
  [[nodiscard]] bool operator<(const ChrSize&other) const;

  std::string name;
  uint64_t size;

  template <typename H>
  friend H AbslHashValue(H h, const ChrSize& c) {
    return H::combine(std::move(h), c.name);
  }
};

class ChrSizeParser {
 public:
  explicit ChrSizeParser(std::string path_to_chr_sizes);
  std::vector<ChrSize> parse(char sep = '\t');

 private:
  std::string _path;
  std::ifstream _f{};
  absl::flat_hash_set<ChrSize> _chrs{};
  std::vector<std::string> _errors{};
};

}  // namespace modle
