#include "modle/parsers.hpp"

#include <charconv>
#include <cinttypes>
#include <filesystem>
#include <limits>
#include <utility>

#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "fmt/printf.h"
#include "modle/utils.hpp"

namespace modle {

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
  std::vector<std::string_view> channels = absl::StrSplit(toks[idx], ',');
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
                      "default case. This should not be possible! BED::size() = {}",
                      this->size()));
  }
}

BEDParser::BEDParser(std::string path_to_bed, BED::Standard bed_standard)
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

BEDParser::BEDParser(std::string_view path_to_bed, BED::Standard bed_standard)
    // For now we always skip the header
    : _path_to_bed(std::string(path_to_bed)),
      _skip_header(true),
      _standard(bed_standard) /*,
      _ncols(_standard) */
{
  this->_fp.open(this->_path_to_bed);
  if (!this->_fp) {
    throw fmt::system_error(errno, "Unable to open file '{}' for reading", this->_path_to_bed);
  }
}

std::vector<BED> BEDParser::parse_all(bool throw_on_duplicates) {
  absl::flat_hash_map<BED, uint64_t> unique_records;
  uint8_t ncols = 0;
  assert(this->_fp.is_open() && this->_fp.good());
  for (auto i = 1UL; std::getline(this->_fp, this->_buff); ++i) {
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
      throw std::runtime_error(fmt::format("Expected %lu fields, got %lu at line %lu of file '%s'.",
                                           ncols, record.size(), i, this->_path_to_bed));
    }
    if (throw_on_duplicates && unique_records.contains(record)) {
      throw std::runtime_error(fmt::format(
          "Detected duplicate entry at line %lu. First occurrence was at line %lu: record: '%s'", i,
          unique_records.at(record), this->_buff));
    }
    unique_records.emplace(std::move(record), i);
  }

  if (!this->_fp && !this->_fp.eof()) {
    throw fmt::system_error(errno, "An error occurred while reading file '{}'", this->_path_to_bed);
  }

  std::vector<BED> records;
  records.reserve(unique_records.size());
  for (const auto& [record, _] : unique_records) {
    records.push_back(record);
  }
  return records;
}

void BEDParser::reset() {
  if (std::filesystem::is_fifo(this->_path_to_bed)) {
    throw std::runtime_error(
        fmt::format("BEDParser::reset() was called on a file that is a FIFO: file path '%s'",
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

ChrSizeParser::ChrSizeParser(std::string path_to_chr_sizes) : _path(std::move(path_to_chr_sizes)) {}
ChrSizeParser::ChrSizeParser(std::string_view path_to_chr_sizes)
    : _path(std::string(path_to_chr_sizes)) {}

std::vector<ChrSize> ChrSizeParser::parse(char sep) {
  std::string buff;
  std::vector<std::string> tokens;
  this->_f = std::ifstream(this->_path);
  for (auto i = 1UL; this->_f.good(); ++i) {
    if (std::getline(this->_f, buff); buff.empty()) {
      continue;
    }

    if (!this->_f) {
      throw fmt::system_error(errno, "IO error while reading file '{}'", this->_path);
    }
    tokens = absl::StrSplit(buff, sep);
    assert(!tokens.empty());  // This should only happen at EOF, which is handled elsewhere
    try {
      if (tokens.size() != 2) {
        throw std::runtime_error(
            fmt::format("Expected 2 tokens, got {}: '{}'", tokens.size(), buff));
      }
      if (ChrSize record(tokens); this->_chrs.contains(record)) {
        throw std::runtime_error(fmt::format("Found multiple records for chr '{}'", record.name));
      } else {
        this->_chrs.emplace(std::move(record));
      }

    } catch (const std::runtime_error& e) {
      this->_errors.push_back(
          fmt::format("Encountered a malformed record at line {} of file '{}': {}.\n "
                      "Line that triggered the error:\n'{}'",
                      i, this->_path, e.what(), buff.data()));
    }
  }
  if (!this->_f && !this->_f.eof()) {
    throw fmt::system_error(errno, "IO error while reading file '{}'", this->_path);
  }
  if (!this->_errors.empty()) {
    throw std::runtime_error(
        fmt::format("The following error(s) occurred while parsing file '%s':\n - %s", this->_path,
                    absl::StrJoin(this->_errors, "\n - ")));
  }
  return {this->_chrs.begin(), this->_chrs.end()};
}

ChrSize::ChrSize(std::vector<std::string>& toks) {
  assert(toks.size() == 2);
  this->name = std::move(toks[0]);
  if (int64_t n = std::stoi(toks[1]); n > 0) {
    this->size = static_cast<uint64_t>(n);
  } else {
    throw std::runtime_error(fmt::format(
        "Sequence '{}' has an invalid length of {} (length cannot be <= 0)", this->name, n));
  }
}

ChrSize::ChrSize(std::string_view chr_name, uint64_t length) : name(chr_name), size(length) {}

bool ChrSize::operator==(const ChrSize& other) const {
  return this->name == other.name && this->size == other.size;
}

bool ChrSize::operator<(const ChrSize& other) const { return this->size < other.size; }

template <typename I, typename R>
Contact::Contact(I b1, I b2, R contacts) {
  static_assert(std::is_integral<I>());
  static_assert(std::is_arithmetic<R>());
  if (b1 < 0) {
    throw std::runtime_error(fmt::format(
        "Exception thrown in Contact constructor: b1 should be a positive integral number, is '{}'",
        b1));
  }
  if (b2 < 0) {
    throw std::runtime_error(fmt::format(
        "Exception thrown in Contact constructor: b2 should be a positive integral number, is '{}'",
        b2));
  }
  if (contacts < 0) {
    throw std::runtime_error(fmt::format(
        "Exception thrown in Contact constructor: b2 should be a positive number, is '{}'",
        contacts));
  }
  this->bin1 = static_cast<uint64_t>(b1);
  this->bin2 = static_cast<uint64_t>(b2);
  this->contacts = static_cast<double>(contacts);
}

ContactsParser::ContactsParser(std::string contact_file) : _path(std::move(contact_file)) {
  if (!std::filesystem::exists(this->_path)) {
    throw std::runtime_error(fmt::format("File '{}' does not exist", this->_path));
  }
  if (std::filesystem::is_directory(this->_path)) {
    throw std::runtime_error(fmt::format("'{}' is a directory. Expected a file", this->_path));
  }
}

bool ContactsParser::parse_next(std::string& buff, Contact& c, std::string_view sep) {
  assert(this->_f);
  if (std::getline(this->_f, buff); buff.empty()) {
    if (this->_f.eof()) {
      return false;
    }
    throw std::logic_error(
        "If you see this error, please report it to the developers on GitHub: "
        "ContactsParser::parse_next() "
        "reached a branch that should not be reachable!");
  }
  if (!this->_f && !this->_f.eof()) {
    throw fmt::system_error(errno, "IO error while reading file {}", this->_path);
  }
  std::vector<std::string_view> tmp = absl::StrSplit(buff, sep);
  if (tmp.size() != 3) {
    throw std::runtime_error(
        fmt::format("Malformed record in file '{}': expected 3 fields, got {} in line '{}'",
                    this->_path, tmp.size(), buff));
  }
  utils::parse_numeric_or_throw(tmp[0], c.bin1);
  utils::parse_numeric_or_throw(tmp[1], c.bin2);
  utils::parse_real_or_throw(tmp[2], c.contacts);
  return false;
}

ContactsParser::MatrixProperties ContactsParser::get_matrix_properties(std::string_view sep) {
  Contact c{};
  ContactsParser::MatrixProperties prop{};
  std::string buff;
  while (this->parse_next(buff, c, sep)) {
    prop.min_bin1 = std::min(prop.min_bin1, c.bin1);
    prop.min_bin2 = std::min(prop.min_bin2, c.bin2);
    prop.max_bin1 = std::max(prop.max_bin1, c.bin1);
    prop.max_bin2 = std::max(prop.max_bin2, c.bin2);
    if (prop.bin_size == 0) {
      prop.bin_size = c.bin1 > c.bin2 ? c.bin1 - c.bin2 : c.bin2 - c.bin1;
      continue;
    }
    if (uint64_t new_bin_size = c.bin1 > c.bin2 ? c.bin1 - c.bin2 : c.bin2 - c.bin1;
        new_bin_size != prop.bin_size) {
      throw std::runtime_error(
          fmt::format("Expected bin of size {}, got {}", prop.bin_size, new_bin_size));
    }
  }
  if (!this->_f && !this->_f.eof()) {
    throw fmt::system_error(errno, "IO error while reading file {}", this->_path);
  }
  this->_f.clear();
  this->_f.seekg(0);
  return prop;
}

ContactMatrix<uint32_t> ContactsParser::parse_into_contact_matrix(uint64_t width,
                                                                  std::string_view sep) {
  std::string buff;
  Contact c{};
  const auto matrix_prop = this->get_matrix_properties(sep);
  ContactMatrix<uint32_t> m((matrix_prop.max_bin1 - matrix_prop.min_bin1) / matrix_prop.bin_size,
                            width / matrix_prop.bin_size);
  this->_f = std::ifstream(this->_path);
  if (!this->_f) {
    throw std::runtime_error(fmt::format("Unable to open file '%s' for reading.", this->_path));
  }
  try {
    while (this->parse_next(buff, c, sep)) {
      uint64_t row = (c.bin1 - matrix_prop.min_bin1) / matrix_prop.bin_size;
      uint64_t col = (c.bin2 - matrix_prop.min_bin2) / matrix_prop.bin_size;
      assert(c.contacts >= 0);
      m.set(row, col, static_cast<uint32_t>(std::lround(c.contacts)));
    }
    if (!this->_f.eof()) {
      throw std::runtime_error("An IO error occurred");
    }
  } catch (const std::exception& err) {
    throw std::runtime_error(
        fmt::format("An error occurred while reading file '%s': %s.", this->_path, err.what()));
  }
  return m;
}

}  // namespace modle