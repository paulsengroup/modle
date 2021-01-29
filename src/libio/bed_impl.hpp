#pragma once

#include <absl/container/flat_hash_map.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_join.h>
#include <absl/strings/str_split.h>
#include <fmt/format.h>  // for format, system_error

#include <cassert>
#include <cerrno>
#include <filesystem>  // for is_fifo
#include <iosfwd>      // for size_t
#include <new>
#include <stdexcept>  // for runtime_error
#include <string>     // string, getline
#include <string_view>
#include <vector>

#include "modle/utils.hpp"

namespace modle::bed {

std::string RGB::to_string() const { return absl::StrCat(r, ",", g, ",", b); }
bool RGB::empty() const { return (static_cast<uint16_t>(r) + g + b) == 0; }

void BED::parse_strand_or_throw(const std::vector<std::string_view>& toks, uint8_t idx,
                                char& field) {
  if (!bed_strand_encoding.contains(toks[idx])) {
    throw std::runtime_error(fmt::format("Unrecognized strand '{}'", toks[idx]));
  }
  field = bed_strand_encoding.at(toks[idx]);
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

BED::BED(std::string_view record, BED::Standard bed_standard) {
  std::vector<std::string_view> toks;
  for (std::string_view tok : absl::StrSplit(record, absl::ByAnyChar("\t "))) {
    if (!tok.empty()) {
      toks.push_back(tok);
    }
  }
  auto ntoks = toks.size();
  if ((bed_standard != BED::Standard::none && ntoks < bed_standard) ||
      (ntoks < BED3 || (ntoks > BED6 && ntoks != BED9 && ntoks < BED12) || ntoks > BED12)) {
    if (bed_standard != BED::Standard::none) {
      throw std::runtime_error(
          fmt::format("Expected {} fields, got {}.\nRecord that caused the error: '{}'",
                      bed_standard, ntoks, record));
    }
    throw std::runtime_error(
        fmt::format("Expected 3, 4, 5, 6 or 12 fields, got {}.\nRefer to "
                    "https://bedtools.readthedocs.io/en/latest/content/"
                    "general-usage.html#bed-format for the BED format specification.\nRecord "
                    "that caused the error: '{}'",
                    ntoks, record));
  }

  this->chrom = toks[BED_CHROM_IDX];
  try {
    utils::parse_numeric_or_throw(toks, BED_CHROM_START_IDX, this->chrom_start);
    utils::parse_numeric_or_throw(toks, BED_CHROM_END_IDX, this->chrom_end);
    if (this->chrom_start > this->chrom_end) {
      throw std::runtime_error(fmt::format(
          "Invalid BED record detected: chrom_start > chrom_end: chrom='{}'; start={}; end={}",
          this->chrom, this->chrom_start, this->chrom_end));
    }
    if (ntoks == BED_CHROM_END || bed_standard == BED3) {
      this->_size = BED3;
      return;
    }
    this->name = toks[BED_NAME_IDX];
    if (ntoks == BED_NAME || bed_standard == BED4) {
      this->_size = BED4;
      return;
    }
    utils::parse_real_or_throw(toks, BED_SCORE_IDX, this->score);
    // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
    if (this->score < 0 || this->score > 1000) {
      throw std::runtime_error(fmt::format(
          "Invalid BED record detected: score field should be between 0 and 1000, is %f.",
          this->score));
    }
    if (ntoks == BED_SCORE || bed_standard == BED5) {
      this->_size = BED5;
      return;
    }
    parse_strand_or_throw(toks, BED_STRAND_IDX, this->strand);
    if (ntoks == BED_STRAND || bed_standard == BED6) {
      this->_size = BED6;
      return;
    }

    utils::parse_numeric_or_throw(toks, BED_THICK_START_IDX, this->thick_start);
    if (this->thick_start < this->chrom_start) {
      throw std::runtime_error(
          fmt::format("Invalid BED record detected: thick_start < chrom_start: chrom='{}'; "
                      "start={}; thick_start={}",
                      this->chrom, this->chrom_start, this->thick_start));
    }
    utils::parse_numeric_or_throw(toks, BED_THICK_END_IDX, this->thick_end);
    if (this->thick_end > this->chrom_end) {
      throw std::runtime_error(
          fmt::format("Invalid BED record detected: thick_end > chrom_end: chrom='{}'; "
                      "start={}; thick_start={}",
                      this->chrom, this->chrom_end, this->thick_end));
    }
    if (this->thick_start > this->thick_end) {
      throw std::runtime_error(
          fmt::format("Invalid BED record detected: thick_start > thick_end: chrom='{}'; "
                      "thick_start={}; thick_end={}",
                      this->chrom, this->chrom_start, this->chrom_end));
    }
    parse_rgb_or_throw(toks, BED_ITEM_RGB_IDX, this->rgb);
    if (ntoks == BED_ITEM_RGB || bed_standard == BED9) {
      this->_size = BED9;
      return;
    }
    utils::parse_numeric_or_throw(toks, BED_BLOCK_COUNT_IDX, this->block_count);
    utils::parse_vect_of_numbers_or_throw(toks, BED_BLOCK_SIZES_IDX, this->block_sizes,
                                          this->block_count);
    utils::parse_vect_of_numbers_or_throw(toks, BED_BLOCK_STARTS_IDX, this->block_starts,
                                          this->block_count);
    this->_size = BED12;
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        "An error occurred while parsing the following BED record '{}':\n  {}", record, e.what()));
  }
}

bool BED::operator==(const BED& other) const {
  return this->chrom == other.chrom && this->chrom_start == other.chrom_start &&
         this->chrom_end == other.chrom_end;
}

bool BED::operator<(const BED& other) const {
  if (this->chrom != other.chrom) {
    return this->chrom < other.chrom;
  }
  if (this->chrom_start != other.chrom_start) {
    return this->chrom_start < other.chrom_start;
  }
  return this->chrom_end < other.chrom_end;
}

uint8_t BED::size() const { return this->_size; }

bool BED::empty() const { return chrom.empty(); }

std::string BED::to_string() const {
  switch (this->size()) {
    case BED3:
      return absl::StrCat(chrom, "\t", chrom_start, "\t", chrom_end);
    case BED4:
      return absl::StrCat(chrom, "\t", chrom_start, "\t", chrom_end, "\t", name);
    case BED5:
      return absl::StrCat(chrom, "\t", chrom_start, "\t", chrom_end, "\t", score);
    case BED6:
      return absl::StrCat(chrom, "\t", chrom_start, "\t", chrom_end, "\t", name, "\t", score, "\t",
                          std::string(1, strand));
    case BED9:
      return absl::StrCat(chrom, "\t", chrom_start, "\t", chrom_end, "\t", name, "\t", score, "\t",
                          std::string(1, strand), "\t", thick_start, "\t", thick_end, "\t",
                          rgb.to_string());
    case BED12:
      return absl::StrCat(chrom, "\t", chrom_start, "\t", chrom_end, "\t", name, "\t", score, "\t",
                          std::string(1, strand), "\t", thick_start, "\t", thick_end, "\t",
                          rgb.to_string(), "\t", block_count, "\t", absl::StrJoin(block_sizes, ","),
                          "\t", absl::StrJoin(block_starts, ","));
    default:
      throw std::runtime_error(
          fmt::format("If you see this error, please report it to the developers on GitHub: "
                      "BED::to_string() reached the "
                      "default case. This should not be possible! BED::size() == {}",
                      this->size()));
  }
}

Parser::Parser(std::string path_to_bed, BED::Standard bed_standard)
    // For now we always skip the header
    : _path_to_bed(std::move(path_to_bed)),
      _skip_header(true),
      _standard(bed_standard) /*,
      _ncols(_standard) */
{
  this->_fp.open(this->_path_to_bed);
  if (!this->_fp) {
    throw fmt::system_error(errno, "Unable to open file '{}' for reading", this->_path_to_bed);
  }
}

Parser::Parser(std::string_view path_to_bed, BED::Standard bed_standard)
    : Parser(std::string{path_to_bed}, bed_standard) {}

std::vector<BED> Parser::parse_n(std::size_t nrecords, bool throw_on_duplicates) {
  absl::flat_hash_map<BED, uint64_t> unique_records;
  if (nrecords != std::numeric_limits<decltype(nrecords)>::max()) {
    unique_records.reserve(nrecords);
  }
  uint8_t ncols = 0;
  assert(this->_fp.is_open() && this->_fp.good());
  for (auto i = 1UL; std::getline(this->_fp, this->_buff) && i < nrecords; ++i) {
    if (this->_buff.empty()) {
      continue;  // Skip empty lines
    }
    if (this->_skip_header && unique_records.empty() &&  // Skip header line(s)
        (this->_buff.front() == '#' || this->_buff.find("track") != std::string::npos ||
         this->_buff.find("browser") != std::string::npos)) {
      continue;
    }
    BED record(this->_buff, this->_standard);
    if (unique_records.empty()) {
      ncols = record.size();
    }
    if (record.size() != ncols) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected {} fields, got {} at line {} of file '{}'"), ncols,
                      record.size(), i, this->_path_to_bed));
    }
    if (throw_on_duplicates && unique_records.contains(record)) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Detected duplicate entry at line {}. First occurrence was at "
                                 "line {} of file '{}': record: '{}'"),
                      i, unique_records.at(record), this->_path_to_bed, this->_buff));
    }
    unique_records.emplace(std::move(record), i);
  }

  if (!this->_fp && !this->_fp.eof()) {
    throw fmt::system_error(errno, "An error occurred while reading file '{}'", this->_path_to_bed);
  }

  std::vector<BED> records(unique_records.size());
  std::transform(unique_records.begin(), unique_records.end(), records.begin(),
                 [](const auto& node) { return node.first; });
  return records;
}

std::string Parser::validate(std::size_t nrecords, bool throw_on_duplicates) {
  try {
    (void)parse_n(nrecords, throw_on_duplicates);
  } catch (const std::runtime_error& e) {
    return e.what();
  }
  return "";
}

std::vector<BED> Parser::parse_all(bool throw_on_duplicates) {
  return parse_n(std::numeric_limits<std::size_t>::max(), throw_on_duplicates);
}

void Parser::reset() {
  if (std::filesystem::is_fifo(this->_path_to_bed)) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("BEDParser::reset() was called on a file that is a FIFO: file path '{}'"),
        this->_path_to_bed));
  }
  if (!this->_fp.is_open()) {
    throw std::runtime_error("BedParser::reset() was called on a closed file!");
  }
  assert(this->_fp.good());
  this->_fp.clear();
  this->_fp.seekg(0);
  this->_buff.clear();
}
}  // namespace modle::bed
