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
  template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
  static void parse_numeric_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                     N& field);
  uint8_t _size;

  // This is a temporary workaround to deal with the fact that libstdc++ 10 does not come with the
  // float/double overloads for std::from_chars
  template <typename R,
            typename = typename std::enable_if<std::is_floating_point<R>::value, R>::type>
  static void parse_real_or_throw(const std::vector<std::string_view>& toks, uint8_t idx, R& field);
  static void parse_strand_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                    char& field);
  static void parse_rgb_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                 RGB& field);

  // Because of the issue outlined for parse_real_or_throw(), this currently only works for integral
  // numbers (which should be fine for our purposes)
  template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
  static void parse_vect_of_numbers_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                             std::vector<N>& field, uint64_t expected_size);
  static void throw_except_from_errc(std::string_view tok, uint8_t idx, std::errc e);
};

class BEDParser {
 public:
  // For now we always skip the header
  explicit BEDParser(std::string path_to_bed, BED::Standard bed_standard = BED::Standard::none);
  explicit BEDParser(std::string_view path_to_bed, BED::Standard bed_standard = BED::Standard::none);
  std::vector<BED> parse_all(bool throw_on_duplicates = true);
  void reset();

 private:
  std::string _path_to_bed;
  std::ifstream _fp;
  std::string _buff{};
  bool _skip_header;
  BED::Standard _standard;
  uint8_t _ncols;
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

}  // namespace modle

#include "modle/impl/parsers.hpp"