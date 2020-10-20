#include "modle/parsers.hpp"

#include <filesystem>

#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "absl/time/clock.h"

namespace modle {

bool SimpleBED::empty() const { return this->chr.empty(); }
std::string SimpleBED::to_string() const {
  return absl::StrFormat("%s\t%lu\t%lu\t%f\t%c\n", this->chr, this->chr_start, this->chr_end,
                         this->score, this->strnd);
}

SimpleBEDParser::SimpleBEDParser(std::string_view path_to_bed) : _path_to_bed(path_to_bed) {
  if (std::filesystem::is_regular_file(this->_path_to_bed) ||
      std::filesystem::is_fifo(this->_path_to_bed)) {
    this->_bed = std::ifstream(this->_path_to_bed);
  } else {
    throw std::runtime_error(absl::StrFormat("Invalid file '%s'.", this->_path_to_bed));
  }
}

SimpleBED SimpleBEDParser::parse_next() {
  do {
    this->_bed.getline(this->_buff.data(), this->_buff.size());
    if (this->_bed.eof()) return {};
    if (this->_bed.bad())
      throw std::runtime_error(
          absl::StrFormat("An IO error while parsing file '%s'.", this->_path_to_bed));
  } while (this->_buff[0] == '\0' || this->_buff[0] == '#');

  std::vector<std::string> toks = absl::StrSplit(this->_buff.begin(), "\t");
  if (toks.size() < 5) {
    absl::FPrintF(stderr, "toks.size()=%lu\n", toks.size());
    throw std::runtime_error(
        absl::StrFormat("Malformed BED record encountered at line %lu in file '%s': expected "
                        "'chrom\tchromStart\tchromEnd\tscore\tstrand', got '%s'.",
                        this->_line, this->_path_to_bed, absl::StrJoin(toks, "\t")));
  }
  ++this->_line;
  SimpleBED record;
  record.chr = std::move(toks[0]);
  record.chr_start = std::stoul(toks[1]);
  record.chr_end = std::stoul(toks[2]);
  record.score = std::stod(toks[3]);
  if (toks[4] != "+" && toks[4] != "-" && toks[4] != ".") {
    throw std::runtime_error(absl::StrFormat(
        "Malformed BED record encountered at line %lu of file '%s': expected one of '+', '-' or "
        "'.' for field 'strand', got '%s'.",
        this->_line, this->_path_to_bed, toks[4]));
  }
  record.strnd = SimpleBED::strand(toks[4].front());

  return record;
}

std::vector<SimpleBED> SimpleBEDParser::parse_all() {
  std::vector<SimpleBED> records;
  while (true) {
    if (auto record = this->parse_next(); !record.empty()) {
      records.emplace_back(std::move(record));
    } else {
      break;
    }
  }
  return records;
}

ChrSizeParser::ChrSizeParser(std::string path_to_chr_sizes) : _path(std::move(path_to_chr_sizes)) {}

std::vector<ChrSize> ChrSizeParser::parse(char sep) {
  std::string buff;
  std::vector<std::string> tokens;
  this->_f = std::ifstream(this->_path);
  for (auto i = 1UL; this->_f.good(); ++i) {
    if (std::getline(this->_f, buff); buff.empty()) {
      if (this->_f.eof()) break;
      continue;
    }

    if (this->_f.bad()) {
      throw std::runtime_error(
          absl::StrFormat("An IO error while reading file '%s'.", this->_path));
    }
    tokens = absl::StrSplit(buff, sep);
    assert(!tokens.empty());  // This should only happen at EOF, which is handled elsewhere
    try {
      if (tokens.size() != 2) {
        throw std::runtime_error(absl::StrFormat("Expected 2 tokens, found %lu.", tokens.size()));
      }
      if (ChrSize record(tokens); this->_chrs.contains(record)) {
        throw std::runtime_error(
            absl::StrFormat("Found multiple entries with id '%s'.", record.name));
      } else {
        this->_chrs.emplace(std::move(record));
      }

    } catch (const std::runtime_error& e) {
      this->_errors.push_back(
          absl::StrFormat("Encountered a malformed record at line '%lu' of file '%s': %s.\n "
                          "Content of the line that triggered the error:\n'%s'",
                          i, this->_path, e.what(), buff.data()));
    }
  }
  this->_f.close();
  if (!this->_errors.empty()) {
    throw std::runtime_error(
        absl::StrFormat("The following error(s) occurred while parsing file '%s':\n - %s",
                        this->_path, absl::StrJoin(this->_errors, "\n - ")));
  }
  return {this->_chrs.begin(), this->_chrs.end()};
}

ChrSize::ChrSize(std::vector<std::string>& toks) {
  assert(toks.size() == 2);
  this->name = std::move(toks[0]);
  if (int64_t n = std::stoi(toks[1]); n > 0) {
    this->size = static_cast<uint64_t>(n);
  } else {
    throw std::runtime_error(
        absl::StrFormat("Sequence '%s' has a length of %d that is <= 0", this->name, n));
  }
}

ChrSize::ChrSize(std::string_view chr_name, uint64_t length) : name(chr_name), size(length) {}

bool ChrSize::operator==(const ChrSize& other) const {
  return this->name == other.name && this->size == other.size;
}

bool ChrSize::operator<(const ChrSize& other) const { return this->size < other.size; }

}  // namespace modle