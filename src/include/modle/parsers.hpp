#pragma once
#include <array>
#include <fstream>
#include <string>
#include <vector>

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

}  // namespace modle
