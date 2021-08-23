#pragma once

#include <absl/container/btree_map.h>      // for btree_map
#include <absl/container/flat_hash_map.h>  // for flat_hash_map
#include <absl/types/span.h>               // for Span
#include <xxh3.h>                          // for XXH3_state_t

#include <boost/filesystem/path.hpp>  // for path
#include <cstddef>                    // IWYU pragma: keep for size_t
#include <cstdint>                    // for uint64_t, uint_fast8_t
#include <fstream>                    // for ifstream
#include <limits>                     // for numeric_limits
#include <memory>                     // for unique_ptr
#include <string>                     // for string, basic_string
#include <string_view>                // for operator""sv, basic_string_view, string_view
#include <utility>                    // IWYU pragma: keep for move, pair
#include <vector>                     // for vector

#include "modle/chrom_sizes.hpp"    // for ChromSizes
#include "modle/common/common.hpp"  // for bp_t
#include "modle/compressed_io.hpp"  // for Reader
#include "modle/interval_tree.hpp"  // for IITree

namespace modle::bed {

struct RGB {
  uint8_t r;
  uint8_t g;
  uint8_t b;

  [[nodiscard]] std::string to_string() const noexcept;
};

struct BED {
  friend class Parser;
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
  explicit BED(Dialect d);
  explicit BED(std::string_view record, size_t id_ = null_id, Dialect bed_standard = autodetect,
               bool enforce_std_compliance = true);
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
  bp_t thick_start{(std::numeric_limits<uint64_t>::max)()};
  bp_t thick_end{(std::numeric_limits<uint64_t>::max)()};
  std::unique_ptr<RGB> rgb{nullptr};
  uint64_t block_count{(std::numeric_limits<uint64_t>::max)()};
  std::vector<uint64_t> block_sizes{};
  std::vector<uint64_t> block_starts{};
  std::string extra_tokens{};

  [[nodiscard]] bool operator==(const BED& other) const noexcept;
  [[nodiscard]] bool operator<(const BED& other) const noexcept;
  [[nodiscard]] size_t id() const noexcept;
  [[nodiscard]] size_t size() const noexcept;
  [[nodiscard]] size_t num_fields() const noexcept;
  [[nodiscard]] Dialect get_standard() const noexcept;
  [[nodiscard]] bool empty() const;
  template <typename H>
  inline friend H AbslHashValue(H h, const BED& c) {
    return H::combine(std::move(h), c.chrom, c.chrom_start, c.chrom_end);
  }

 private:
  static constexpr auto null_id = (std::numeric_limits<size_t>::max)();
  size_t _id{null_id};
  Dialect _standard{none};
  static void parse_rgb_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                 RGB& field);
  [[nodiscard]] static RGB parse_rgb_or_throw(const std::vector<std::string_view>& toks,
                                              uint8_t idx);
  static void parse_strand_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                    char& field);
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

  [[nodiscard]] uint64_t hash(XXH3_state_t* state,
                              uint64_t seed = 17039577131913730910ULL) const;  // NOLINT
};

template <typename K = std::string, typename I = bp_t>
class BED_tree {
  using IITree_t = IITree<I, BED>;
  using BED_tree_t = absl::btree_map<K, IITree<I, BED>>;

  friend IITree_t;

 public:
  inline BED_tree() = default;
  inline explicit BED_tree(const boost::filesystem::path& path_to_bed,
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
  [[nodiscard]] inline bool contains_overlap(const K& chrom_name, uint64_t chrom_start,
                                             uint64_t chrom_end) const;

  [[nodiscard]] inline size_t count_overlaps(const BED& interval) const;
  [[nodiscard]] inline size_t count_overlaps(const K& chrom_name, uint64_t chrom_start,
                                             uint64_t chrom_end) const;

  [[nodiscard]] inline absl::Span<const BED> find_overlaps(const BED& interval) const;
  [[nodiscard]] inline absl::Span<const BED> find_overlaps(const K& chrom_name,
                                                           uint64_t chrom_start,
                                                           uint64_t chrom_end) const;

  [[nodiscard]] inline bool empty() const;
  [[nodiscard]] inline size_t size() const;
  [[nodiscard]] inline size_t size(const K& chrom_name) const;

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
  explicit Parser(boost::filesystem::path path_to_bed,
                  BED::Dialect bed_standard = BED::Dialect::autodetect,
                  bool enforce_std_compliance = true);

  [[nodiscard]] BED parse_next();
  [[nodiscard]] std::vector<BED> parse_n(size_t num_records);
  [[nodiscard]] BED_tree<> parse_n_in_interval_tree(size_t num_records);
  [[nodiscard]] std::string validate(size_t nrecords = 100);  // NOLINT
  [[nodiscard]] std::vector<BED> parse_all();
  [[nodiscard]] BED_tree<> parse_all_in_interval_tree();
  void reset();

 private:
  compressed_io::Reader _reader{};
  BED::Dialect _dialect;
  bool _enforce_std_compliance;
  std::string _buff{};
  size_t _num_records_parsed{0};
  size_t _num_lines_read{0};

  size_t skip_header();
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

template <>
struct fmt::formatter<modle::bed::RGB> {
  inline constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  inline auto format(const modle::bed::RGB& rgb, FormatContext& ctx) -> decltype(ctx.out());
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
  inline auto format(const modle::bed::BED& b, FormatContext& ctx) -> decltype(ctx.out());

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

#include "../../bed_impl.hpp"  // IWYU pragma: export
