#pragma once
#include <array>
#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"
#include "modle/contacts.hpp"

namespace modle {

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

class BEDParser {
 public:
  // For now we always skip the header
  explicit BEDParser(std::string path_to_bed, BED::Standard bed_standard = BED::Standard::none);
  explicit BEDParser(std::string_view path_to_bed,
                     BED::Standard bed_standard = BED::Standard::none);
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

struct ChrSize {
  ChrSize(std::string_view chr_name, uint64_t length);
  explicit ChrSize(std::vector<std::string>& toks);

  [[nodiscard]] bool operator==(const ChrSize& other) const;
  [[nodiscard]] bool operator<(const ChrSize& other) const;

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
  explicit ChrSizeParser(std::string_view path_to_chr_sizes);
  std::vector<ChrSize> parse(char sep = '\t');

 private:
  std::string _path;
  std::ifstream _f{};
  absl::flat_hash_set<ChrSize> _chrs{};
  std::vector<std::string> _errors{};
};

struct Contact {
  Contact() = default;
  template <typename I, typename R>
  Contact(I b1, I b2, R contacts);
  uint64_t bin1;
  uint64_t bin2;
  double contacts;
};

class ContactsParser {
 public:
  explicit ContactsParser(std::string contact_file);
  ContactMatrix<uint32_t> parse_into_contact_matrix(uint64_t width, std::string_view sep = "\t");

 private:
  std::string _path;
  std::ifstream _f{};

  struct MatrixProperties {
    uint64_t min_bin1{UINT64_MAX};
    uint64_t max_bin1{0};
    uint64_t min_bin2{UINT64_MAX};
    uint64_t max_bin2{0};
    uint64_t bin_size{0};
  };

  [[nodiscard]] MatrixProperties get_matrix_properties(std::string_view sep);
  bool parse_next(std::string& buff, Contact& c, std::string_view sep);
};

}  // namespace modle