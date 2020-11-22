#pragma once
#include <cstdint>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/hash/hash.h"

namespace modle::bed {

using namespace std::literals::string_view_literals;
static const absl::flat_hash_map<std::string_view, char> bed_strand_encoding{
    {"+"sv, '+'},       {"plus"sv, '+'},    {"fwd"sv, '+'},     {"Fwd"sv, '+'},
    {"forward"sv, '+'}, {"Forward"sv, '+'}, {"FWD"sv, '+'},     {"FORWARD"sv, '+'},
    {"-"sv, '-'},       {"minus"sv, '-'},   {"rev"sv, '-'},     {"Rev"sv, '-'},
    {"reverse"sv, '-'}, {"Reverse"sv, '-'}, {"REV"sv, '-'},     {"REVERSE"sv, '-'},
    {"."sv, '.'},       {""sv, '.'},        {"none"sv, '.'},    {"None"sv, '.'},
    {"NONE"sv, '.'},    {"unknown"sv, '.'}, {"Unknown"sv, '.'}, {"unk"sv, '.'},
    {"Unk"sv, '.'},     {"UNK"sv, '.'}};

struct RGB {
  uint8_t r;
  uint8_t g;
  uint8_t b;

  [[nodiscard]] std::string to_string() const;
  [[nodiscard]] bool empty() const;
};

struct BED {
  enum Standard { BED3 = 3U, BED4 = 4U, BED5 = 5U, BED6 = 6U, BED9 = 9U, BED12 = 12U, none = 13U };
  enum Fields {
    BED_CHROM = 1,
    BED_CHROM_START = 2,
    BED_CHROM_END = 3,
    BED_NAME = 4,
    BED_SCORE = 5,
    BED_STRAND = 6,
    BED_THICK_START = 7,
    BED_THICK_END = 8,
    BED_ITEM_RGB = 9,
    BED_BLOCK_COUNT = 10,
    BED_BLOCK_SIZES = 11,
    BED_BLOCK_STARTS = 12
  };
  enum FieldsIdx {
    BED_CHROM_IDX = 0,
    BED_CHROM_START_IDX = 1,
    BED_CHROM_END_IDX = 2,
    BED_NAME_IDX = 3,
    BED_SCORE_IDX = 4,
    BED_STRAND_IDX = 5,
    BED_THICK_START_IDX = 6,
    BED_THICK_END_IDX = 7,
    BED_ITEM_RGB_IDX = 8,
    BED_BLOCK_COUNT_IDX = 9,
    BED_BLOCK_SIZES_IDX = 10,
    BED_BLOCK_STARTS_IDX = 11
  };
  BED() = default;
  explicit BED(std::string_view record, Standard bed_standard = none);
  std::string chrom{};
  uint64_t chrom_start{};
  uint64_t chrom_end{};
  std::string name{};
  double score{};
  char strand{'.'};  // TODO: change this to an enum or something. Also deal with
  // non-standard encodings (Forward, Reverse etc.)
  uint64_t thick_start{};
  uint64_t thick_end{};
  RGB rgb{};
  uint64_t block_count{};
  std::vector<uint64_t> block_sizes{};
  std::vector<uint64_t> block_starts{};

  [[nodiscard]] bool operator==(const BED& other) const;
  [[nodiscard]] bool operator<(const BED& other) const;
  [[nodiscard]] uint8_t size() const;
  [[nodiscard]] std::string to_string() const;
  [[nodiscard]] bool empty() const;
  template <typename H>
  friend H AbslHashValue(H h, const BED& c) {
    return H::combine(std::move(h), c.chrom, c.chrom_start, c.chrom_end);
  }

 private:
  uint8_t _size;
  static void parse_rgb_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                 RGB& field);
  static void parse_strand_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                    char& field);
};

class Parser {
 public:
  // For now we always skip the header
  explicit Parser(std::string path_to_bed, BED::Standard bed_standard = BED::Standard::none);
  explicit Parser(std::string_view path_to_bed, BED::Standard bed_standard = BED::Standard::none);
  std::vector<BED> parse_all(bool throw_on_duplicates = true);
  void reset();

 private:
  std::string _path_to_bed;
  std::ifstream _fp;
  std::string _buff{};
  bool _skip_header;
  BED::Standard _standard;
  // uint8_t _ncols;
};
}  // namespace modle::bed