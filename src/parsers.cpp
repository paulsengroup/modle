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
  const auto t0 = absl::Now();
  while (true) {
    if (auto record = this->parse_next(); !record.empty()) {
      records.emplace_back(std::move(record));
    } else {
      break;
    }
  }
  absl::FPrintF(stderr, "Parsed %lu records in %s\n", this->_line,
                absl::FormatDuration(absl::Now() - t0));
  return records;
}

}  // namespace modle