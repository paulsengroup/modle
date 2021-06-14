#pragma once

#include <absl/container/flat_hash_map.h>  // for flat_hash_map

#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for uint64_t, uint_fast8_t
#include <fstream>      // for ifstream
#include <limits>       // for numeric_limits
#include <memory>       // for unique_ptr
#include <string>       // for string, basic_string
#include <string_view>  // for operator""sv, basic_string_view, string_view
#include <utility>      // IWYU pragma: keep for move
#include <vector>       // for vector

namespace modle::bed {

struct RGB {
  uint8_t r;
  uint8_t g;
  uint8_t b;

  [[nodiscard]] std::string to_string() const noexcept;
};

struct BED {
  enum Dialect : uint_fast8_t {
    BED3 = 3U,
    BED4 = 4U,
    BED5 = 5U,
    BED6 = 6U,
    BED9 = 9U,
    BED12 = 12U,
    autodetect = 253U,
    none = 254U
  };

  [[nodiscard]] static Dialect str_to_dialect(std::string_view s);
  [[nodiscard]] static std::string dialect_to_str(Dialect d);

  enum FieldsIdx : uint_fast8_t {
    BED_CHROM_IDX = 0U,
    BED_CHROM_START_IDX = 1U,
    BED_CHROM_END_IDX = 2U,
    BED_NAME_IDX = 3U,
    BED_SCORE_IDX = 4U,
    BED_STRAND_IDX = 5U,
    BED_THICK_START_IDX = 6U,
    BED_THICK_END_IDX = 7U,
    BED_ITEM_RGB_IDX = 8U,
    BED_BLOCK_COUNT_IDX = 9U,
    BED_BLOCK_SIZES_IDX = 10U,
    BED_BLOCK_STARTS_IDX = 11U
  };

  enum Fields : uint_fast8_t {
    BED_CHROM = BED_CHROM_IDX + 1,
    BED_CHROM_START = BED_CHROM_START_IDX + 1,
    BED_CHROM_END = BED_CHROM_END_IDX + 1,
    BED_NAME = BED_NAME_IDX + 1,
    BED_SCORE = BED_SCORE_IDX + 1,
    BED_STRAND = BED_STRAND_IDX + 1,
    BED_THICK_START = BED_THICK_START_IDX + 1,
    BED_THICK_END = BED_THICK_END_IDX + 1,
    BED_ITEM_RGB = BED_ITEM_RGB_IDX + 1,
    BED_BLOCK_COUNT = BED_BLOCK_COUNT_IDX + 1,
    BED_BLOCK_SIZES = BED_BLOCK_SIZES_IDX + 1,
    BED_BLOCK_STARTS = BED_BLOCK_STARTS_IDX + 1
  };

  BED() = default;
  explicit BED(std::string_view record, Dialect bed_standard = autodetect,
               bool enforce_std_compliance = true);
  BED(const BED& other);
  BED(BED&& other) = default;
  ~BED() = default;

  BED& operator=(const BED& other);
  BED& operator=(BED&& other) = default;

  std::string chrom{};
  uint64_t chrom_start{};
  uint64_t chrom_end{};
  std::string name{};
  double score{};
  char strand{'.'};
  uint64_t thick_start{std::numeric_limits<uint64_t>::max()};
  uint64_t thick_end{std::numeric_limits<uint64_t>::max()};
  std::unique_ptr<RGB> rgb{nullptr};
  uint64_t block_count{std::numeric_limits<uint64_t>::max()};
  std::vector<uint64_t> block_sizes{};
  std::vector<uint64_t> block_starts{};
  std::string extra_tokens{};

  [[nodiscard]] bool operator==(const BED& other) const noexcept;
  [[nodiscard]] bool operator<(const BED& other) const noexcept;
  [[nodiscard]] size_t size() const noexcept;
  [[nodiscard]] Dialect get_standard() const noexcept;
  [[nodiscard]] std::string to_string() const noexcept;
  [[nodiscard]] bool empty() const;
  template <typename H>
  inline friend H AbslHashValue(H h, const BED& c) {
    return H::combine(std::move(h), c.chrom, c.chrom_start, c.chrom_end);
  }

 private:
  Dialect _standard;
  static void parse_rgb_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                 RGB& field);
  [[nodiscard]] static RGB parse_rgb_or_throw(const std::vector<std::string_view>& toks,
                                              uint8_t idx);
  static void parse_strand_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                    char& field);
  [[nodiscard]] static Dialect detect_standard(const std::vector<std::string_view>& toks);
  static void validate_record(const std::vector<std::string_view>& toks, Dialect standard);

  // parse_* methods that return a bool will return true if there are no more tokens to be processed
  void parse_chrom(const std::vector<std::string_view>& toks);
  void parse_chrom_start(const std::vector<std::string_view>& toks);
  bool parse_chrom_end(const std::vector<std::string_view>& toks);
  bool parse_name(const std::vector<std::string_view>& toks);
  bool parse_score(const std::vector<std::string_view>& toks, bool validate);
  bool parse_strand(const std::vector<std::string_view>& toks);
  bool parse_thick_start(const std::vector<std::string_view>& toks, bool validate);
  bool parse_thick_end(const std::vector<std::string_view>& toks, bool validate);
  bool parse_item_rgb(const std::vector<std::string_view>& toks);
  bool parse_block_count(const std::vector<std::string_view>& toks);
  bool parse_block_sizes(const std::vector<std::string_view>& toks);
  bool parse_block_starts(const std::vector<std::string_view>& toks);
  void parse_extra_tokens(const std::vector<std::string_view>& toks);
};

class Parser {
 public:
  // For now we always skip the header
  explicit Parser(std::string path_to_bed, BED::Dialect bed_standard = BED::Dialect::autodetect,
                  bool enforce_std_compliance = true);
  explicit Parser(std::string_view path_to_bed,
                  BED::Dialect bed_standard = BED::Dialect::autodetect,
                  bool enforce_std_compliance = true);
  [[nodiscard]] std::vector<BED> parse_n(size_t nrecords, bool throw_on_duplicates = true);
  [[nodiscard]] std::string validate(size_t nrecords = 100,  // NOLINT
                                     bool throw_on_duplicates = true);
  [[nodiscard]] std::vector<BED> parse_all(bool throw_on_duplicates = true);
  void reset();

 private:
  std::string _path_to_bed;
  std::ifstream _fp;
  std::string _buff{};
  bool _skip_header;
  BED::Dialect _standard;
  bool _enforce_std_compliance;
};

using namespace std::literals::string_view_literals;
static const absl::flat_hash_map<std::string_view, char> bed_strand_encoding{
    {"+"sv, '+'},       {"plus"sv, '+'},    {"fwd"sv, '+'},     {"Fwd"sv, '+'},
    {"forward"sv, '+'}, {"Forward"sv, '+'}, {"FWD"sv, '+'},     {"FORWARD"sv, '+'},
    {"-"sv, '-'},       {"minus"sv, '-'},   {"rev"sv, '-'},     {"Rev"sv, '-'},
    {"reverse"sv, '-'}, {"Reverse"sv, '-'}, {"REV"sv, '-'},     {"REVERSE"sv, '-'},
    {"."sv, '.'},       {""sv, '.'},        {"none"sv, '.'},    {"None"sv, '.'},
    {"NONE"sv, '.'},    {"unknown"sv, '.'}, {"Unknown"sv, '.'}, {"unk"sv, '.'},
    {"Unk"sv, '.'},     {"UNK"sv, '.'}};

static const std::vector<std::string_view> bed_dialects{"BED3"sv, "BED4"sv, "BED5"sv,
                                                        "BED6"sv, "BED9"sv, "BED12"sv};
static const absl::flat_hash_map<std::string_view, BED::Dialect> str_to_bed_dialect_mappings{
    {"BED3"sv, BED::Dialect::BED3}, {"BED4"sv, BED::Dialect::BED4},
    {"BED5"sv, BED::Dialect::BED5}, {"BED6"sv, BED::Dialect::BED6},
    {"BED9"sv, BED::Dialect::BED9}, {"BED12"sv, BED::Dialect::BED12}};
static const absl::flat_hash_map<BED::Dialect, std::string_view> bed_dialect_to_str_mappings{
    {BED::Dialect::BED3, "BED3"sv}, {BED::Dialect::BED4, "BED4"sv},
    {BED::Dialect::BED5, "BED5"sv}, {BED::Dialect::BED6, "BED6"sv},
    {BED::Dialect::BED9, "BED9"sv}, {BED::Dialect::BED12, "BED12"sv}};

}  // namespace modle::bed
