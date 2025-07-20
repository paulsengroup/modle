// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>  // for btree_map
#include <absl/types/span.h>           // for Span
#include <fmt/format.h>                // for format_parse_context, formatter
#include <xxhash.h>                    // for XXH3_state_t, XXH_INLINE_XXH3_state_t

#include <array>        // for array
#include <filesystem>   // for path
#include <limits>       // for numeric_limits
#include <memory>       // for unique_ptr
#include <string>       // for string
#include <string_view>  // for operator""sv, string_view, basic_string_view, stri...
#include <type_traits>  // for __strip_reference_wrapper<>::__type
#include <utility>      // for make_pair, pair
#include <vector>       // for vector

#include "modle/common/common.hpp"     // for bp_t, u64, u8
#include "modle/common/const_map.hpp"  // for ConstMap, ConstMap::ConstMap<Key, Value, Size>
#include "modle/compressed_io/compressed_io.hpp"  // for Reader
#include "modle/interval_tree.hpp"                // for IITree, IITree::IITree<I, T>

namespace modle::bed {

struct RGB {
  u8 r;
  u8 g;
  u8 b;

  bool operator==(const RGB& other) const noexcept;
  bool operator!=(const RGB& other) const noexcept;
  [[nodiscard]] std::string to_string() const;
};

struct BED {
  friend class Parser;
  enum Dialect : u8f {
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

  enum FieldsIdx : u8f {
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

  enum Fields : u8f {
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
  explicit BED(Dialect d);
  explicit BED(std::string_view record, usize id_ = null_id, Dialect bed_standard = autodetect,
               bool validate = true);
  BED(std::string_view chrom_, bp_t chrom_start_, bp_t chrom_end_);
  BED(const BED& other);
  BED(BED&& other) = default;
  ~BED() = default;

  BED& operator=(const BED& other);
  BED& operator=(BED&& other) = default;

  std::string chrom{};
  bp_t chrom_start{};
  bp_t chrom_end{};
  std::string name{};
  double score{};
  char strand{'.'};
  bp_t thick_start{(std::numeric_limits<bp_t>::max)()};
  bp_t thick_end{(std::numeric_limits<bp_t>::max)()};
  std::unique_ptr<RGB> rgb{nullptr};
  u64 block_count{(std::numeric_limits<u64>::max)()};
  std::vector<u64> block_sizes{};
  std::vector<u64> block_starts{};
  std::string extra_tokens{};

  [[nodiscard]] constexpr explicit operator bool() const noexcept;
  [[nodiscard]] bool operator==(const BED& other) const noexcept;
  [[nodiscard]] bool operator<(const BED& other) const noexcept;
  [[nodiscard]] constexpr usize id() const noexcept;
  [[nodiscard]] constexpr usize size() const noexcept;
  [[nodiscard]] usize num_fields() const noexcept;
  [[nodiscard]] constexpr Dialect get_standard() const noexcept;
  [[nodiscard]] bool empty() const;
  template <typename H>
  inline friend H AbslHashValue(H h, const BED& c) {
    return H::combine(std::move(h), c.chrom, c.chrom_start, c.chrom_end);
  }

 private:
  static constexpr auto null_id = (std::numeric_limits<usize>::max)();
  usize _id{null_id};
  Dialect _standard{none};
  static void parse_rgb_or_throw(const std::vector<std::string_view>& toks, u8 idx, RGB& field);
  [[nodiscard]] static RGB parse_rgb_or_throw(const std::vector<std::string_view>& toks, u8 idx);
  static void parse_strand_or_throw(const std::vector<std::string_view>& toks, u8 idx, char& field);
  [[nodiscard]] static Dialect detect_standard(std::string_view line);
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

  [[nodiscard]] u64 hash(XXH3_state_t* state, u64 seed = 17039577131913730910ULL) const;
};

template <typename K = std::string, typename I = bp_t>
class BED_tree {
  using IITree_t = IITree<I, BED>;
  using BED_tree_t = phmap::btree_map<K, IITree<I, BED>>;

  friend IITree_t;

 public:
  inline BED_tree() = default;
  inline explicit BED_tree(const std::filesystem::path& path_to_bed,
                           BED::Dialect dialect = BED::Dialect::autodetect);

  using node_type = typename BED_tree_t::node_type;
  using key_type = typename BED_tree_t::value_type::first_type;
  using value_type = typename BED_tree_t::value_type::second_type;
  using iterator = typename BED_tree_t::iterator;
  using const_iterator = typename BED_tree_t::const_iterator;

  inline std::pair<iterator, bool> insert(BED interval);
  inline std::pair<iterator, bool> insert(const K& chrom_name, value_type tree);
  inline std::pair<iterator, bool> insert(const K& chrom_name, I chrom_start, I chrom_end);
  inline std::pair<iterator, bool> emplace(BED&& interval);
  inline std::pair<iterator, bool> emplace(const K& chrom_name, value_type&& tree);

  inline void insert(absl::Span<const BED> intervals);
  inline void emplace(std::vector<BED>&& intervals);

  [[nodiscard]] inline const value_type& at(const K& chrom) const;

  inline void index();
  inline void index(const K& chrom_name);

  [[nodiscard]] inline bool contains(const K& chrom_name) const;
  [[nodiscard]] inline bool contains_overlap(const BED& interval) const;
  [[nodiscard]] inline bool contains_overlap(const K& chrom_name, u64 chrom_start,
                                             u64 chrom_end) const;

  [[nodiscard]] inline usize count_overlaps(const BED& interval) const;
  [[nodiscard]] inline usize count_overlaps(const K& chrom_name, u64 chrom_start,
                                            u64 chrom_end) const;

  [[nodiscard]] inline absl::Span<const BED> find_overlaps(const BED& interval) const;
  [[nodiscard]] inline absl::Span<const BED> find_overlaps(const K& chrom_name, u64 chrom_start,
                                                           u64 chrom_end) const;

  [[nodiscard]] inline bool empty() const;
  [[nodiscard]] inline usize size() const;
  [[nodiscard]] inline usize size(const K& chrom_name) const;

  inline void clear();
  inline void clear(const K& chrom_name);

  [[nodiscard]] inline iterator begin();
  [[nodiscard]] inline iterator end();
  [[nodiscard]] inline const_iterator begin() const;
  [[nodiscard]] inline const_iterator end() const;
  [[nodiscard]] inline const_iterator cbegin() const;
  [[nodiscard]] inline const_iterator cend() const;

  // TODO: add erase methods

 private:
  BED_tree_t _trees{};
};

class Parser {
 public:
  explicit Parser(const std::filesystem::path& path_to_bed,
                  BED::Dialect bed_standard = BED::Dialect::autodetect,
                  bool enforce_std_compliance = true);

  [[nodiscard]] BED parse_next();
  [[nodiscard]] std::vector<BED> parse_n(usize num_records);
  [[nodiscard]] BED_tree<> parse_n_in_interval_tree(usize num_records);
  [[nodiscard]] std::string validate(usize nrecords = 100);
  [[nodiscard]] std::vector<BED> parse_all();
  [[nodiscard]] BED_tree<> parse_all_in_interval_tree();
  void reset();

 private:
  compressed_io::Reader _reader{};
  BED::Dialect _dialect;
  bool _enforce_std_compliance;
  std::string _buff{};
  usize _num_records_parsed{0};
  usize _num_lines_read{0};

  usize skip_header();
};

using namespace std::literals::string_view_literals;

// clang-format off
static constexpr utils::ConstMap<std::string_view, char, 26> bed_strand_encoding{
    {"+"sv, '+'},       {"plus"sv, '+'},
    {"fwd"sv, '+'},     {"Fwd"sv, '+'},
    {"forward"sv, '+'}, {"Forward"sv, '+'},
    {"FWD"sv, '+'},     {"FORWARD"sv, '+'},
    {"-"sv, '-'},       {"minus"sv, '-'},
    {"rev"sv, '-'},     {"Rev"sv, '-'},
    {"reverse"sv, '-'}, {"Reverse"sv, '-'},
    {"REV"sv, '-'},     {"REVERSE"sv, '-'},
    {"."sv, '.'},       {""sv, '.'},
    {"none"sv, '.'},    {"None"sv, '.'},
    {"NONE"sv, '.'},    {"unknown"sv, '.'},
    {"Unknown"sv, '.'}, {"unk"sv, '.'},
    {"Unk"sv, '.'},     {"UNK"sv, '.'}};
// clang-format on

static constexpr std::array<std::string_view, 6> bed_dialects{"BED3"sv, "BED4"sv, "BED5"sv,
                                                              "BED6"sv, "BED9"sv, "BED12"sv};
static constexpr utils::ConstMap<std::string_view, BED::Dialect, 6> str_to_bed_dialect_mappings{
    {"BED3"sv, BED::Dialect::BED3}, {"BED4"sv, BED::Dialect::BED4},
    {"BED5"sv, BED::Dialect::BED5}, {"BED6"sv, BED::Dialect::BED6},
    {"BED9"sv, BED::Dialect::BED9}, {"BED12"sv, BED::Dialect::BED12}};

static constexpr utils::ConstMap<BED::Dialect, std::string_view, 6> bed_dialect_to_str_mappings{
    {BED::Dialect::BED3, "BED3"sv}, {BED::Dialect::BED4, "BED4"sv},
    {BED::Dialect::BED5, "BED5"sv}, {BED::Dialect::BED6, "BED6"sv},
    {BED::Dialect::BED9, "BED9"sv}, {BED::Dialect::BED12, "BED12"sv}};

}  // namespace modle::bed

template <>
struct fmt::formatter<modle::bed::RGB> {
  inline constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  inline auto format(const modle::bed::RGB& rgb, FormatContext& ctx) const -> decltype(ctx.out());
};

template <>
struct fmt::formatter<modle::bed::BED> {
  // Presentation can be any of the following:
  //  - auto
  //  - none
  //  - bed3
  //  - bed4
  //  - bed5
  //  - bed6
  //  - bed9
  //  - bed12

  inline constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  inline auto format(const modle::bed::BED& b, FormatContext& ctx) const -> decltype(ctx.out());

 private:
  modle::bed::BED::Dialect presentation{modle::bed::BED::none};
  constexpr static std::array<std::pair<std::string_view, modle::bed::BED::Dialect>, 8>
      presentation_mappings{std::make_pair("auto", modle::bed::BED::autodetect),
                            std::make_pair("none", modle::bed::BED::none),
                            std::make_pair("bed3", modle::bed::BED::BED3),
                            std::make_pair("bed4", modle::bed::BED::BED4),
                            std::make_pair("bed5", modle::bed::BED::BED5),
                            std::make_pair("bed6", modle::bed::BED::BED6),
                            std::make_pair("bed9", modle::bed::BED::BED9),
                            std::make_pair("bed12", modle::bed::BED::BED12)};
};

#include "../../../../bed_impl.hpp"  // IWYU pragma: export
