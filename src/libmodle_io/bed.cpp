// IWYU pragma: no_include <ext/alloc_traits.h>

#include "modle/bed.hpp"

#include <absl/container/flat_hash_map.h>  // for flat_hash_map, BitMask, operator!=
#include <absl/hash/hash.h>                // for Hash
#include <absl/strings/str_cat.h>          // for StrCat
#include <absl/strings/str_join.h>         // for StrJoin
#include <absl/strings/str_split.h>        // for StrSplit, Splitter, SplitIterator, operator!=
#include <fmt/format.h>                    // for format, FMT_STRING, system_error
#include <fmt/ostream.h>

#include <algorithm>                         // for transform, _max
#include <boost/filesystem/file_status.hpp>  // for filesystem::file_type
#include <boost/filesystem/operations.hpp>   // filesystem::status
#include <boost/filesystem/path.hpp>         // for filesystem::path
#include <cassert>                           // for assert
#include <cerrno>                            // for errno
#include <cstdint>                           // for uint64_t, uint8_t, uint16_t
#include <limits>                            // for numeric_limits
#include <stdexcept>                         // for runtime_error
#include <string>       // for string, allocator, char_traits, getline, opera...
#include <string_view>  // for string_view, operator==, basic_string_view
#include <utility>      // for move
#include <vector>       // for vector

#include "absl/strings/match.h"
#include "modle/common/suppress_compiler_warnings.hpp"
#include "modle/common/utils.hpp"   // for parse_numeric_or_throw, parse_vect_of_numbers_...
#include "modle/compressed_io.hpp"  // for Reader, getline
#include "modle/interval_tree.hpp"  // for IITree

namespace modle::bed {

std::string RGB::to_string() const noexcept { return fmt::to_string(*this); }

void BED::parse_strand_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                char& field) {
  const auto match = bed_strand_encoding.find(toks[idx]);
  if (match == bed_strand_encoding.end()) {
    throw std::runtime_error(fmt::format("Unrecognized strand '{}'", toks[idx]));
  }
  field = match->second;
}

void BED::parse_rgb_or_throw(const std::vector<std::string_view>& toks, uint8_t idx, RGB& field) {
  if (toks[idx] == "0") {
    field = RGB{0, 0, 0};
    return;
  }
  const std::vector<std::string_view> channels = absl::StrSplit(toks[idx], ',');
  if (channels.size() != 3) {
    throw std::runtime_error(
        fmt::format("RGB: expected 3 fields, got {}: '{}'", channels.size(), toks[idx]));
  }
  utils::parse_numeric_or_throw(channels, 0, field.r);
  utils::parse_numeric_or_throw(channels, 1, field.g);
  utils::parse_numeric_or_throw(channels, 2, field.b);
}

RGB BED::parse_rgb_or_throw(const std::vector<std::string_view>& toks, uint8_t idx) {
  RGB buff;  // NOLINT
  BED::parse_rgb_or_throw(toks, idx, buff);
  return buff;
}

BED::Dialect BED::detect_standard(std::string_view line) {
  return BED::detect_standard(absl::StrSplit(line, absl::ByAnyChar("\t ")));
}

BED::Dialect BED::detect_standard(const std::vector<std::string_view>& toks) {
  if (toks.size() < BED3) {
    throw std::runtime_error(
        fmt::format("Expected at least 3 fields, got {}.\nRecord that caused the error: '{}'",
                    toks.size(), absl::StrJoin(toks, "\t")));
  }

  switch (toks.size()) {
    case BED3:
      return BED3;
    case BED4:
      return BED4;
    case BED5:
      return BED5;
    case BED6:
      return BED6;
    case BED9:
      return BED9;
    case BED12:
      return BED12;
    default:
      return none;
  }
}

void BED::validate_record(const std::vector<std::string_view>& toks, const Dialect standard) {
  assert(toks.size() >= BED3);  // NOLINT
  if (standard == none) {
    return;
  }

  const auto& ntoks = toks.size();
  if (standard == autodetect && detect_standard(toks) == none) {
    throw std::runtime_error(fmt::format(
        "Invalid BED record detected: expected 3, 4, 5, 6 or 12 fields, got {}.\nRefer to "
        "https://bedtools.readthedocs.io/_end/latest/content/general-usage.html#bed-format for the "
        "BED format specification",
        ntoks));
  }

  if (detect_standard(toks) < standard) {
    throw std::runtime_error(fmt::format(
        "Invalid BED record detected: Expected BED record with at least {} fields, got {}",
        standard, ntoks));
  }
}

void BED::parse_chrom(const std::vector<std::string_view>& toks) {
  this->chrom = toks[BED_CHROM_IDX];
}

void BED::parse_chrom_start(const std::vector<std::string_view>& toks) {
  utils::parse_numeric_or_throw(toks, BED_CHROM_START_IDX, this->chrom_start);
}

bool BED::parse_chrom_end(const std::vector<std::string_view>& toks) {
  utils::parse_numeric_or_throw(toks, BED_CHROM_END_IDX, this->chrom_end);
  if (this->chrom_start > this->chrom_end) {
    throw std::runtime_error(
        fmt::format("Invalid BED record detected: chrom_start > chrom_end: chrom='{}'; "
                    "overlap_start={}; _end={}",
                    this->chrom, this->chrom_start, this->chrom_end));
  }
  return this->_standard == BED3;
}

bool BED::parse_name(const std::vector<std::string_view>& toks) {
  assert(this->_standard >= BED4);  // NOLINT
  this->name = toks[BED_NAME_IDX];
  return this->_standard == BED4;
}

bool BED::parse_score(const std::vector<std::string_view>& toks, bool validate) {
  assert(this->_standard >= BED5);  // NOLINT
  utils::parse_numeric_or_throw(toks, BED_SCORE_IDX, this->score);
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  if (this->_standard != none && validate && (this->score < 0 || this->score > 1000)) {
    throw std::runtime_error(fmt::format(
        "Invalid BED record detected: score field should be between 0.0 and 1000.0, is {}.",
        this->score));
  }
  return this->_standard == BED5;
}

bool BED::parse_strand(const std::vector<std::string_view>& toks) {
  assert(this->_standard >= BED6);  // NOLINT
  parse_strand_or_throw(toks, BED_STRAND_IDX, this->strand);
  return this->_standard == BED6;
}

bool BED::parse_thick_start(const std::vector<std::string_view>& toks, bool validate) {
  assert(this->_standard >= 7);  // NOLINT
  utils::parse_numeric_or_throw(toks, BED_THICK_START_IDX, this->thick_start);
  if (validate && this->thick_start < this->chrom_start) {
    throw std::runtime_error(
        fmt::format("Invalid BED record detected: thick_start < chrom_start: chrom='{}'; "
                    "overlap_start={}; thick_start={}",
                    this->chrom, this->chrom_start, this->thick_start));
  }
  return this->_standard == none && toks.size() == BED_THICK_START;
}

bool BED::parse_thick_end(const std::vector<std::string_view>& toks, bool validate) {
  assert(this->_standard >= 8);  // NOLINT
  utils::parse_numeric_or_throw(toks, BED_THICK_END_IDX, this->thick_end);
  if (validate && this->thick_end > this->chrom_end) {
    throw std::runtime_error(
        fmt::format("Invalid BED record detected: thick_end > chrom_end: chrom='{}'; "
                    "overlap_start={}; thick_start={}",
                    this->chrom, this->chrom_end, this->thick_end));
  }
  if (validate && this->thick_start > this->thick_end) {
    throw std::runtime_error(
        fmt::format("Invalid BED record detected: thick_start > thick_end: chrom='{}'; "
                    "thick_start={}; thick_end={}",
                    this->chrom, this->chrom_start, this->chrom_end));
  }
  return this->_standard == none && toks.size() == BED_THICK_END;
}

bool BED::parse_item_rgb(const std::vector<std::string_view>& toks) {
  assert(this->_standard >= BED9);  // NOLINT
  this->rgb = std::make_unique<RGB>(parse_rgb_or_throw(toks, BED_ITEM_RGB_IDX));
  return this->_standard == BED9;
}

bool BED::parse_block_count(const std::vector<std::string_view>& toks) {
  assert(this->_standard >= 10);  // NOLINT
  utils::parse_numeric_or_throw(toks, BED_BLOCK_COUNT_IDX, this->block_count);
  return this->_standard == none && toks.size() == BED_BLOCK_COUNT;
}

bool BED::parse_block_sizes(const std::vector<std::string_view>& toks) {
  assert(this->_standard >= 11);  // NOLINT
  utils::parse_vect_of_numbers_or_throw(toks, BED_BLOCK_SIZES_IDX, this->block_sizes,
                                        this->block_count);
  return this->_standard == none && toks.size() == BED_BLOCK_SIZES;
}
bool BED::parse_block_starts(const std::vector<std::string_view>& toks) {
  assert(this->_standard >= 12);  // NOLINT
  utils::parse_vect_of_numbers_or_throw(toks, BED_BLOCK_STARTS_IDX, this->block_starts,
                                        this->block_count);
  return this->_standard == BED12;
}
void BED::parse_extra_tokens(const std::vector<std::string_view>& toks) {
  assert(this->_standard == none);  // NOLINT
  assert(toks.size() >= BED12);     // NOLINT
  // Copy non-whitespace tokens
  this->extra_tokens = absl::StrJoin(toks.begin() + BED12, toks.end(), "\t");
}

BED::Dialect BED::str_to_dialect(std::string_view s) {
  assert(bed::str_to_bed_dialect_mappings.contains(s));  // NOLINT
  return bed::str_to_bed_dialect_mappings.at(s);
}

std::string BED::dialect_to_str(Dialect d) {
  assert(bed::bed_dialect_to_str_mappings.contains(d));  // NOLINT
  return std::string{bed::bed_dialect_to_str_mappings.at(d)};
}

BED::BED(std::string_view chrom_, bp_t chrom_start_, bp_t chrom_end_)
    : chrom(chrom_), chrom_start(chrom_start_), chrom_end(chrom_end_), _standard(BED3) {}

BED::BED(BED::Dialect d) : _standard(d) {}

BED::BED(std::string_view record, size_t id_, BED::Dialect bed_standard, bool validate) : _id(id_) {
  std::vector<std::string_view> toks;
  for (std::string_view tok : absl::StrSplit(record, absl::ByAnyChar("\t "))) {
    if (!tok.empty()) {
      toks.push_back(tok);
    }
  }
  if (bed_standard == autodetect) {
    if ((bed_standard = detect_standard(toks)) == none && validate) {
      throw std::runtime_error(
          fmt::format("Expected 3, 4, 5, 6 or 12 fields, got {}.\nRefer to "
                      "https://bedtools.readthedocs.io/_end/latest/content/"
                      "general-usage.html#bed-format for the BED format specification.\nRecord "
                      "that caused the error: '{}'",
                      toks.size(), absl::StrJoin(toks, "\t")));
    }
  }

  validate_record(toks, bed_standard);
  this->_standard = bed_standard;

  try {
    this->parse_chrom(toks);
    this->parse_chrom_start(toks);
    if (this->parse_chrom_end(toks)) {
      return;
    }
    if (this->parse_name(toks)) {
      return;
    }
    if (this->parse_score(toks, validate)) {
      return;
    }
    if (this->parse_strand(toks)) {
      return;
    }
    if (this->parse_thick_start(toks, validate)) {
      assert(bed_standard == none);  // NOLINT
      return;
    }
    if (this->parse_thick_end(toks, validate)) {
      assert(bed_standard == none);  // NOLINT
      return;
    }
    if (this->parse_item_rgb(toks)) {
      return;
    }
    if (this->parse_block_count(toks)) {
      assert(bed_standard == none);  // NOLINT
      return;
    }
    if (this->parse_block_sizes(toks)) {
      assert(bed_standard == none);  // NOLINT
      return;
    }
    if (this->parse_block_starts(toks)) {
      return;
    }
    this->parse_extra_tokens(toks);

  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        "An error occurred while parsing the following BED record '{}':\n  {}", record, e.what()));
  }
}

BED::BED(const BED& other)
    : chrom(other.chrom),
      chrom_start(other.chrom_start),
      chrom_end(other.chrom_end),
      name(other.name),
      score(other.score),
      strand(other.strand),
      thick_start(other.thick_start),
      thick_end(other.thick_end),
      rgb(other.rgb ? std::make_unique<RGB>(*other.rgb) : nullptr),
      block_count(other.block_count),
      block_sizes(other.block_sizes),
      block_starts(other.block_starts),
      extra_tokens(other.extra_tokens),
      _id(other._id),
      _standard(other._standard) {}

BED& BED::operator=(const BED& other) {
  if (&other == this) {
    return *this;
  }
  chrom = other.chrom;
  chrom_start = other.chrom_start;
  chrom_end = other.chrom_end;
  name = other.name;
  score = other.score;
  strand = other.strand;
  thick_start = other.thick_start;
  thick_end = other.thick_end;
  rgb = other.rgb ? std::make_unique<RGB>(*other.rgb) : nullptr;
  block_count = other.block_count;
  block_sizes = other.block_sizes;
  block_starts = other.block_starts;
  extra_tokens = other.extra_tokens;
  _id = other._id;
  _standard = other._standard;

  return *this;
}

bool BED::operator==(const BED& other) const noexcept {
  return this->chrom == other.chrom && this->chrom_start == other.chrom_start &&
         this->chrom_end == other.chrom_end;
}

bool BED::operator<(const BED& other) const noexcept {
  if (this->chrom != other.chrom) {
    return this->chrom < other.chrom;
  }
  if (this->chrom_start != other.chrom_start) {
    return this->chrom_start < other.chrom_start;
  }
  return this->chrom_end < other.chrom_end;
}

BED::Dialect BED::get_standard() const noexcept { return this->_standard; }

size_t BED::id() const noexcept { return this->_id; }

size_t BED::size() const noexcept { return this->chrom_end - this->chrom_start; }

size_t BED::num_fields() const noexcept {
  assert(this->_standard != autodetect);  // NOLINT
  if (this->_standard != none) {
    return static_cast<size_t>(this->_standard);
  }
  if (thick_end == -1ULL) {  // NOLINT
    return BED_THICK_START;
  }
  if (!rgb) {
    return BED_THICK_END;
  }
  assert(block_count != -1ULL);  // NOLINT
  if (block_sizes.empty()) {
    return BED_BLOCK_COUNT;
  }
  if (block_starts.empty()) {
    return BED_BLOCK_SIZES;
  }
  assert(!extra_tokens.empty());  // NOLINT
  return BED12 + static_cast<size_t>(std::count(extra_tokens.begin(), extra_tokens.end(), '\t'));
}

bool BED::empty() const { return chrom.empty(); }

std::string BED::to_string() const noexcept { return fmt::to_string(*this); }

uint64_t BED::hash(XXH_INLINE_XXH3_state_t* state, uint64_t seed) const {
  auto handle_errors = [&](const auto& status) {
    if (status == XXH_ERROR || !state) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Failed to hash the following BED record: {}"), this->to_string()));
    }
  };

  handle_errors(XXH3_64bits_reset_withSeed(state, seed));
  handle_errors(XXH3_64bits_update(state, this->chrom.data(), this->chrom.size() * sizeof(char)));
  handle_errors(XXH3_64bits_update(state, &this->chrom_start, sizeof(decltype(this->chrom_start))));
  handle_errors(XXH3_64bits_update(state, &this->chrom_end, sizeof(decltype(this->chrom_end))));

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  return static_cast<uint64_t>(XXH3_64bits_digest(state));
  DISABLE_WARNING_POP
}

Parser::Parser(boost::filesystem::path path_to_bed, BED::Dialect bed_standard,
               bool enforce_std_compliance)
    // For now we always skip the header
    : _reader(std::move(path_to_bed)),
      _dialect(bed_standard),
      _enforce_std_compliance(enforce_std_compliance),
      _num_lines_read(this->skip_header()) {
  if (this->_dialect == BED::autodetect) {
    this->_dialect = BED::detect_standard(this->_buff);
  }
}

BED Parser::parse_next() {
  if (this->_reader.eof() && this->_buff.empty()) {
    return BED{};
  }

  BED record{this->_buff, this->_num_records_parsed++, this->_dialect,
             this->_enforce_std_compliance};
  if (this->_reader.eof()) {
    this->_buff.clear();
    return record;
  }

  // Look for the next non-empty line
  while (this->_reader.getline(this->_buff)) {
    ++this->_num_lines_read;
    if (!this->_buff.empty()) {
      break;
    }
  }

  return record;
}

std::vector<BED> Parser::parse_n(size_t num_records) {
  if (this->_reader.path().empty()) {
    return std::vector<BED>{};
  }
  assert(this->_reader.is_open());  // NOLINT

  using record_idx_t = size_t;
  using line_num_t = size_t;
  absl::flat_hash_map<BED, std::pair<record_idx_t, line_num_t>> records;

  if (num_records != std::numeric_limits<decltype(num_records)>::max()) {
    records.reserve(num_records);
  }

  for (auto record = this->parse_next(); !record.empty() && this->_num_records_parsed < num_records;
       record = this->parse_next()) {
    if (this->_dialect != BED::none && record.num_fields() < this->_dialect) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Expected BED record with at least {} fields, got {} at line {} of file {}"),
          this->_dialect, record.size(), this->_num_lines_read, this->_reader.path()));
    }

    if (auto [node, new_insertion] = records.try_emplace(
            record, std::make_pair(this->_num_records_parsed - 1, this->_num_lines_read));
        !new_insertion) {
      const auto& other_record = node->first;
      const auto& other_record_line = node->second.second;
      throw std::runtime_error(fmt::format(
          FMT_STRING("Detected duplicate record at line {} of file {}. First occurrence was at "
                     "line {}.\n - First occurrence:  '{}'\n - Second occurrence: '{}'"),
          this->_num_lines_read, this->_reader.path(), other_record_line, record.to_string(),
          other_record.to_string()));
    }
  }

  std::vector<BED> _records(records.size());
  for (const auto& [record, idx] : records) {
    assert(idx.first < _records.size());  // NOLINT
    _records[idx.first] = record;
  }

  return _records;
}

BED_tree<> Parser::parse_n_in_interval_tree(size_t num_records) {
  if (this->_reader.path().empty()) {
    return BED_tree<>{};
  }
  assert(this->_reader.is_open());  // NOLINT

  using line_num_t = size_t;
  absl::flat_hash_map<BED, line_num_t> records;
  BED_tree<> intervals;

  for (auto record = this->parse_next(); !record.empty() && this->_num_records_parsed < num_records;
       record = this->parse_next()) {
    if (this->_dialect != BED::none && record.num_fields() < this->_dialect) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Expected BED record with at least {} fields, got {} at line {} of file {}"),
          this->_dialect, record.size(), this->_num_lines_read, this->_reader.path()));
    }

    if (auto [node, new_insertion] = records.try_emplace(record, this->_num_lines_read);
        !new_insertion) {
      const auto& other_record_line = node->second;
      const auto& other_record = node->first;
      throw std::runtime_error(fmt::format(
          FMT_STRING("Detected duplicate record at line {} of file {}. First occurrence was at "
                     "line {}.\n - First occurrence:  '{}'\n - Second occurrence: '{}'"),
          this->_num_lines_read, this->_reader.path(), other_record_line, record.to_string(),
          other_record.to_string()));
    }

    intervals.emplace(std::move(record));
  }
  intervals.index();

  return intervals;
}

std::string Parser::validate(size_t nrecords) {
  try {
    (void)parse_n(nrecords);
  } catch (const std::runtime_error& e) {
    return e.what();
  }
  return "";
}

std::vector<BED> Parser::parse_all() { return parse_n(std::numeric_limits<size_t>::max()); }

BED_tree<> Parser::parse_all_in_interval_tree() {
  return parse_n_in_interval_tree(std::numeric_limits<size_t>::max());
}

void Parser::reset() {
  if (boost::filesystem::status(this->_reader.path()).type() == boost::filesystem::fifo_file) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("BED::Parser::reset() was called on a file that is a FIFO: file path '{}'"),
        this->_reader.path()));
  }
  if (!this->_reader.is_open()) {
    throw std::runtime_error("BED::Parser::reset() was called on a closed file!");
  }

  try {
    this->_reader.reset();
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "An error occourred while seeking to the begin of file {} using BED::Parser::reset()"),
        this->_reader.path()));
  }
  this->_buff.clear();
  this->_num_lines_read = 0;
  this->_num_records_parsed = 0;

  this->skip_header();
}

size_t Parser::skip_header() {
  if (!this->_reader.is_open()) {
    return 0;
  }
  assert(this->_num_records_parsed == 0);  // NOLINT
  auto num_header_lines = 0ULL;
  assert(this->_reader.is_open());  // NOLINT
  while (this->_reader.getline(this->_buff)) {
    if (this->_buff.empty()) {  // Skip empty lines
      ++num_header_lines;
      continue;
    }
    if (this->_buff.front() == '#' || absl::StrContains(this->_buff, "track") ||
        absl::StrContains(this->_buff, "browser")) {  // Skip header line(s)
      ++num_header_lines;
      continue;
    }
    return num_header_lines;
  }

  return num_header_lines;
}

}  // namespace modle::bed
