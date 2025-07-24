// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// IWYU pragma: no_include <ext/alloc_traits.h>

#include "modle/bed/bed.hpp"

#include <absl/strings/str_split.h>
#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/std.h>

#include <algorithm>
#include <cassert>
#include <exception>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/common/fmt_helpers.hpp"
#include "modle/common/numeric_utils.hpp"
#include "modle/common/string_utils.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"
#include "modle/common/utils.hpp"
#include "modle/compressed_io/compressed_io.hpp"

namespace modle::bed {

std::string RGB::to_string() const { return fmt::to_string(*this); }

bool RGB::operator==(const modle::bed::RGB& other) const noexcept {
  return r == other.r && g == other.g && b == other.b;
}

bool RGB::operator!=(const modle::bed::RGB& other) const noexcept { return !(*this == other); }

void BED::parse_strand_or_throw(const std::vector<std::string_view>& toks, u8 idx, char& field) {
  const auto tok = utils::strip_quote_pairs(toks[idx]);
  const auto match = bed_strand_encoding.find(tok);
  if (match == bed_strand_encoding.end()) {
    throw std::runtime_error(fmt::format("unrecognized strand \"{}\"", tok));
  }
  field = *match.second;
}

void BED::parse_rgb_or_throw(const std::vector<std::string_view>& toks, u8 idx, RGB& field) {
  const auto tok = utils::strip_quote_pairs(toks[idx]);
  if (tok == "0") {
    field = RGB{0, 0, 0};
    return;
  }
  const std::vector<std::string_view> channels = absl::StrSplit(tok, ',');
  if (channels.size() != 3) {
    throw std::runtime_error(
        fmt::format("RGB: expected 3 fields, got {}: \"{}\"", channels.size(), tok));
  }
  utils::parse_numeric_or_throw(channels, 0, field.r);
  utils::parse_numeric_or_throw(channels, 1, field.g);
  utils::parse_numeric_or_throw(channels, 2, field.b);
}

RGB BED::parse_rgb_or_throw(const std::vector<std::string_view>& toks, u8 idx) {
  RGB buff{};
  BED::parse_rgb_or_throw(toks, idx, buff);
  return buff;
}

BED::Dialect BED::detect_standard(std::string_view line) {
  return BED::detect_standard(absl::StrSplit(line, absl::ByAnyChar("\t ")));
}

BED::Dialect BED::detect_standard(const std::vector<std::string_view>& toks) {
  if (toks.size() < BED3) {
    throw std::runtime_error(
        fmt::format("expected at least 3 fields, got {}.\nRecord that caused the error: \"{}\"",
                    toks.size(), fmt::join(toks, "\t")));
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
  assert(toks.size() >= BED3);
  if (standard == none) {
    return;
  }

  const auto& ntoks = toks.size();
  if (standard == autodetect && detect_standard(toks) == none) {
    throw std::runtime_error(fmt::format(
        "Invalid BED record detected: expected 3, 4, 5, 6 or 12 fields, got {}.\nRefer to "
        "https://bedtools.readthedocs.io/_end/latest/content/general-usage.html#bed-format for "
        "the BED format specification",
        ntoks));
  }

  if (detect_standard(toks) < standard) {
    throw std::runtime_error(fmt::format(
        "Invalid BED record detected: Expected BED record with at least {} fields, got {}",
        static_cast<std::underlying_type_t<Dialect>>(standard), ntoks));
  }
}

void BED::parse_chrom(const std::vector<std::string_view>& toks) {
  chrom = utils::strip_quote_pairs(toks[BED_CHROM_IDX]);
}

void BED::parse_chrom_start(const std::vector<std::string_view>& toks) {
  utils::parse_numeric_or_throw(toks, BED_CHROM_START_IDX, chrom_start);
}

bool BED::parse_chrom_end(const std::vector<std::string_view>& toks) {
  utils::parse_numeric_or_throw(toks, BED_CHROM_END_IDX, chrom_end);
  if (chrom_start > chrom_end) {
    throw std::runtime_error(
        fmt::format("Invalid BED record detected: chrom_start > chrom_end: chrom=\"{}\"; "
                    "start={}; end={}",
                    chrom, chrom_start, chrom_end));
  }
  return _standard == BED3;
}

bool BED::parse_name(const std::vector<std::string_view>& toks) {
  assert(_standard >= BED4);
  name = utils::strip_quote_pairs(toks[BED_NAME_IDX]);
  return _standard == BED4;
}

bool BED::parse_score(const std::vector<std::string_view>& toks, bool validate) {
  assert(_standard >= BED5);
  utils::parse_numeric_or_throw(toks, BED_SCORE_IDX, score);
  if (_standard != none && validate && (score < 0 || score > 1000)) {
    throw std::runtime_error(fmt::format(
        "Invalid BED record detected: score field should be between 0.0 and 1000.0, is {}.",
        score));
  }
  return _standard == BED5;
}

bool BED::parse_strand(const std::vector<std::string_view>& toks) {
  assert(_standard >= BED6);
  parse_strand_or_throw(toks, BED_STRAND_IDX, strand);
  return _standard == BED6;
}

bool BED::parse_thick_start(const std::vector<std::string_view>& toks, bool validate) {
  assert(_standard >= 7);
  utils::parse_numeric_or_throw(toks, BED_THICK_START_IDX, thick_start);
  if (validate && thick_start < chrom_start) {
    throw std::runtime_error(
        fmt::format("Invalid BED record detected: thick_start < chrom_start: chrom=\"{}\"; "
                    "overlap_start={}; thick_start={}",
                    chrom, chrom_start, thick_start));
  }
  return _standard == none && toks.size() == BED_THICK_START;
}

bool BED::parse_thick_end(const std::vector<std::string_view>& toks, bool validate) {
  assert(_standard >= 8);
  utils::parse_numeric_or_throw(toks, BED_THICK_END_IDX, thick_end);
  if (validate && thick_end > chrom_end) {
    throw std::runtime_error(
        fmt::format("Invalid BED record detected: thick_end > chrom_end: chrom=\"{}\"; "
                    "overlap_start={}; thick_start={}",
                    chrom, chrom_end, thick_end));
  }
  if (validate && thick_start > thick_end) {
    throw std::runtime_error(
        fmt::format("Invalid BED record detected: thick_start > thick_end: chrom=\"{}\"; "
                    "thick_start={}; thick_end={}",
                    chrom, chrom_start, chrom_end));
  }
  return _standard == none && toks.size() == BED_THICK_END;
}

bool BED::parse_item_rgb(const std::vector<std::string_view>& toks) {
  assert(_standard >= BED9);
  rgb = std::make_unique<RGB>(parse_rgb_or_throw(toks, BED_ITEM_RGB_IDX));
  return _standard == BED9;
}

bool BED::parse_block_count(const std::vector<std::string_view>& toks) {
  assert(_standard >= 10);
  utils::parse_numeric_or_throw(toks, BED_BLOCK_COUNT_IDX, block_count);
  return _standard == none && toks.size() == BED_BLOCK_COUNT;
}

bool BED::parse_block_sizes(const std::vector<std::string_view>& toks) {
  assert(_standard >= 11);
  utils::parse_vect_of_numbers_or_throw(toks, BED_BLOCK_SIZES_IDX, block_sizes, block_count);
  return _standard == none && toks.size() == BED_BLOCK_SIZES;
}
bool BED::parse_block_starts(const std::vector<std::string_view>& toks) {
  assert(_standard >= 12);
  utils::parse_vect_of_numbers_or_throw(toks, BED_BLOCK_STARTS_IDX, block_starts, block_count);
  return _standard == BED12;
}
void BED::parse_extra_tokens(const std::vector<std::string_view>& toks) {
  assert(_standard == none);
  assert(toks.size() >= BED12);
  // Copy non-whitespace tokens
  extra_tokens = fmt::format(FMT_COMPILE("{}"), fmt::join(toks.begin() + BED12, toks.end(), "\t"));
}

BED::Dialect BED::str_to_dialect(std::string_view s) {
  assert(bed::str_to_bed_dialect_mappings.contains(s));
  return bed::str_to_bed_dialect_mappings.at(s);
}

std::string BED::dialect_to_str(Dialect d) {
  assert(bed::bed_dialect_to_str_mappings.contains(d));
  return std::string{bed::bed_dialect_to_str_mappings.at(d)};
}

BED::BED(std::string_view chrom_, bp_t chrom_start_, bp_t chrom_end_)
    : chrom(chrom_), chrom_start(chrom_start_), chrom_end(chrom_end_), _standard(BED3) {}

BED::BED(BED::Dialect d) : _standard(d) {}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
BED::BED(std::string_view record, usize id_, BED::Dialect bed_standard, bool validate) : _id(id_) {
  std::vector<std::string_view> toks;
  for (std::string_view tok :
       absl::StrSplit(strip_trailing_whitespace(record), absl::ByAnyChar("\t "))) {
    if (!tok.empty()) {
      toks.push_back(tok);
    }
  }
  if (bed_standard == autodetect) {
    if ((bed_standard = detect_standard(toks)) == none && validate) {
      throw std::runtime_error(fmt::format(
          "Expected 3, 4, 5, 6 or 12 fields, got {}.\n"
          "Refer to "
          "https://bedtools.readthedocs.io/_end/latest/content/general-usage.html#bed-format "
          "for the BED format specification.\n"
          "Record that caused the error: \"{}\"",
          toks.size(), fmt::join(toks, "\t")));
    }
  }

  validate_record(toks, bed_standard);
  _standard = bed_standard;

  try {
    parse_chrom(toks);
    parse_chrom_start(toks);
    if (parse_chrom_end(toks)) {
      return;
    }
    if (parse_name(toks)) {
      return;
    }
    if (parse_score(toks, validate)) {
      return;
    }
    if (parse_strand(toks)) {
      return;
    }
    if (parse_thick_start(toks, validate)) {
      assert(bed_standard == none);
      return;
    }
    if (parse_thick_end(toks, validate)) {
      assert(bed_standard == none);
      return;
    }
    if (parse_item_rgb(toks)) {
      return;
    }
    if (parse_block_count(toks)) {
      assert(bed_standard == none);
      return;
    }
    if (parse_block_sizes(toks)) {
      assert(bed_standard == none);
      return;
    }
    if (parse_block_starts(toks)) {
      return;
    }
    parse_extra_tokens(toks);

  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format("An error occurred while parsing the following BED record \"{}\":\n  {}",
                    record, e.what()));
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
  return chrom == other.chrom && chrom_start == other.chrom_start && chrom_end == other.chrom_end;
}

bool BED::operator<(const BED& other) const noexcept {
  if (chrom != other.chrom) {
    return chrom < other.chrom;
  }
  if (chrom_start != other.chrom_start) {
    return chrom_start < other.chrom_start;
  }
  return chrom_end < other.chrom_end;
}

usize BED::num_fields() const noexcept {
  assert(_standard != autodetect);
  if (_standard != none) {
    return static_cast<usize>(_standard);
  }
  if (thick_end == -1ULL) {
    return BED_THICK_START;
  }
  if (!rgb) {
    return BED_THICK_END;
  }
  assert(block_count != -1ULL);
  if (block_sizes.empty()) {
    return BED_BLOCK_COUNT;
  }
  if (block_starts.empty()) {
    return BED_BLOCK_SIZES;
  }
  assert(!extra_tokens.empty());
  return BED12 + static_cast<usize>(std::count(extra_tokens.begin(), extra_tokens.end(), '\t'));
}

bool BED::empty() const { return chrom.empty(); }

u64 BED::hash(XXH3_state_t* state, u64 seed) const {
  auto handle_errors = [&](const auto& status) {
    if (status == XXH_ERROR || !state) {
      throw std::runtime_error(fmt::format("Failed to hash the following BED record: {}", *this));
    }
  };

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USED_BUT_MARKED_UNUSED
  handle_errors(XXH3_64bits_reset_withSeed(state, seed));
  handle_errors(XXH3_64bits_update(state, chrom.data(), chrom.size() * sizeof(char)));
  handle_errors(XXH3_64bits_update(state, &chrom_start, sizeof(decltype(chrom_start))));
  handle_errors(XXH3_64bits_update(state, &chrom_end, sizeof(decltype(chrom_end))));

  return utils::conditional_static_cast<u64>(XXH3_64bits_digest(state));
  DISABLE_WARNING_POP
}

Parser::Parser(const std::filesystem::path& path_to_bed, BED::Dialect bed_standard,
               bool enforce_std_compliance)
    // For now we always skip the header
    : _reader(path_to_bed),
      _dialect(bed_standard),
      _enforce_std_compliance(enforce_std_compliance),
      _num_lines_read(skip_header()) {
  if (_dialect == BED::autodetect) {
    try {
      _dialect = BED::detect_standard(_buff);
    } catch (const std::runtime_error& e) {
      throw std::runtime_error(
          fmt::format("An error occurred while parsing file {}: {}", path_to_bed, e.what()));
    }
  }
}

BED Parser::parse_next() {
  if (_reader.eof() && _buff.empty()) {
    return BED{};
  }

  BED record{_buff, _num_records_parsed++, _dialect, _enforce_std_compliance};
  if (_reader.eof()) {
    _buff.clear();
    return record;
  }

  // Look for the next non-empty line
  while (_reader.getline(_buff)) {
    ++_num_lines_read;
    if (!_buff.empty()) {
      break;
    }
  }

  return record;
}

std::vector<BED> Parser::parse_n(usize num_records) {
  if (_reader.path().empty()) {
    return std::vector<BED>{};
  }
  assert(_reader.is_open());

  struct RecordMetadata {
    usize record_idx;
    usize line_num;
  };

  phmap::btree_map<BED, RecordMetadata> records;
  for (auto record = parse_next(); !record.empty() && _num_records_parsed < num_records;
       record = parse_next()) {
    if (_dialect != BED::none && record.num_fields() < _dialect) {
      throw std::runtime_error(
          fmt::format("Expected BED record with at least {} fields, got {} at line {} of file {}",
                      static_cast<std::underlying_type_t<BED::Dialect>>(_dialect), record.size(),
                      _num_lines_read, _reader.path()));
    }

    if (auto [node, new_insertion] = records.try_emplace(
            std::move(record), RecordMetadata{_num_records_parsed - 1, _num_lines_read});
        !new_insertion) {
      const auto& other_record = node->first;
      const auto& other_record_line = node->second.line_num;
      throw std::runtime_error(
          fmt::format("Detected duplicate record at line {} of file {}. First occurrence was at "
                      "line {}.\n - First occurrence:  \"{}\"\n - Second occurrence: \"{}\"",
                      _num_lines_read, _reader.path(), other_record_line, record, other_record));
    }
  }

  std::vector<BED> _records(records.size());
  for (const auto& [record, idx] : records) {
    assert(idx.record_idx < _records.size());
    _records[idx.record_idx] = record;
  }

  return _records;
}

BED_tree<> Parser::parse_n_in_interval_tree(usize num_records) {
  if (_reader.path().empty()) {
    return BED_tree<>{};
  }
  assert(_reader.is_open());

  using line_num_t = usize;
  phmap::btree_map<BED, line_num_t> records;
  BED_tree<> intervals;

  for (auto record = parse_next(); !record.empty() && _num_records_parsed < num_records;
       record = parse_next()) {
    if (_dialect != BED::none && record.num_fields() < _dialect) {
      throw std::runtime_error(
          fmt::format("Expected BED record with at least {} fields, got {} at line {} of file {}",
                      static_cast<std::underlying_type_t<BED::Dialect>>(_dialect), record.size(),
                      _num_lines_read, _reader.path()));
    }

    if (auto [node, new_insertion] = records.try_emplace(record, _num_lines_read); !new_insertion) {
      const auto& other_record_line = node->second;
      const auto& other_record = node->first;
      throw std::runtime_error(
          fmt::format("Detected duplicate record at line {} of file {}. First occurrence was at "
                      "line {}.\n - First occurrence:  \"{}\"\n - Second occurrence: \"{}\"",
                      _num_lines_read, _reader.path(), other_record_line, record, other_record));
    }

    intervals.emplace(std::move(record));
  }
  intervals.index();

  return intervals;
}

std::string Parser::validate(usize nrecords) {
  try {
    std::ignore = parse_n(nrecords);
  } catch (const std::runtime_error& e) {
    return e.what();
  }
  return "";
}

std::vector<BED> Parser::parse_all() { return parse_n((std::numeric_limits<usize>::max)()); }

BED_tree<> Parser::parse_all_in_interval_tree() {
  return parse_n_in_interval_tree((std::numeric_limits<usize>::max)());
}

void Parser::reset() {
  if (std::filesystem::status(_reader.path()).type() == std::filesystem::file_type::fifo) {
    throw std::runtime_error(
        fmt::format("BED::Parser::reset() was called on a file that is a FIFO: file path \"{}\"",
                    _reader.path()));
  }
  if (!_reader.is_open()) {
    throw std::runtime_error("BED::Parser::reset() was called on a closed file!");
  }

  try {
    _reader.reset();
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format("An error occurred while seeking to the begin "
                    "of file {} using BED::Parser::reset(): {}",
                    _reader.path(), e.what()));
  }
  _buff.clear();
  _num_lines_read = 0;
  _num_records_parsed = 0;

  skip_header();
}

usize Parser::skip_header() {
  if (!_reader.is_open()) {
    return 0;
  }
  assert(_num_records_parsed == 0);
  usize num_header_lines = 0L;
  assert(_reader.is_open());
  while (_reader.getline(_buff)) {
    if (_buff.empty()) {  // Skip empty lines
      ++num_header_lines;
      continue;
    }
    if (_buff.front() == '#' || str_contains(_buff, "track") ||
        str_contains(_buff, "browser")) {  // Skip header line(s)
      ++num_header_lines;
      continue;
    }
    return num_header_lines;
  }

  return num_header_lines;
}

}  // namespace modle::bed
