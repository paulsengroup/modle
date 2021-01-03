#pragma once

#include <filesystem>

#include "absl/strings/str_split.h"
#include "catch2/catch.hpp"
#include "modle/bed.hpp"

namespace modle::test::bed {
using namespace modle::bed;

inline void compare_bed_records_with_file(std::vector<BED> records, const std::string &bed_file) {
  auto fp = std::ifstream(bed_file);
  std::vector<std::string> lines;
  lines.reserve(records.size());
  std::string buff;
  while (std::getline(fp, buff)) {
    if (buff.front() == '#' || buff.find("track") != std::string::npos ||
        buff.find("browser") != std::string::npos)
      continue;
    lines.push_back(buff);
  }

  REQUIRE(records.size() == lines.size());
  std::sort(records.begin(), records.end());
  std::sort(lines.begin(), lines.end(), [&](const std::string_view &a, const std::string_view &b) {
    const std::vector<std::string_view> toksa = absl::StrSplit(a, '\t');
    const std::vector<std::string_view> toksb = absl::StrSplit(b, '\t');
    const auto &chra = toksa[0];
    const auto &chrb = toksb[0];
    const auto &starta = toksa[1];
    const auto &startb = toksb[1];
    const auto &enda = toksa[2];
    const auto &endb = toksb[2];
    if (chra != chrb) return chra < chrb;
    if (starta != startb)
      return std::stoull(starta.data(), nullptr) < std::stoull(startb.data(), nullptr);
    assert(enda != endb);
    return std::stoull(enda.data(), nullptr) < std::stoull(endb.data(), nullptr);
  });

  for (auto i = 0UL; i < records.size(); ++i) {
    CHECK(records[i].to_string() == lines[i]);
  }
}

TEST_CASE("BED Parser simple", "[parsers][BED][short]") {
  const std::string bed_file = "../../test/data/unit_tests/sample.bed6";
  REQUIRE(std::filesystem::exists(bed_file));
  auto p = bed::Parser(bed_file);
  auto records = p.parse_all();
  CHECK(records.size() == 9);
  std::sort(records.begin(), records.end());
  CHECK(records[0].chrom == "chr7");
  CHECK(records[0].chrom_start == 127471196);
  CHECK(records[0].chrom_end == 127472363);
  CHECK(records[0].score == 0);
  CHECK(records[0].strand == '+');
}

TEST_CASE("BED Parser simple: BED6 -> BED3", "[parsers][BED][short]") {
  const std::string bed_file = "../../test/data/unit_tests/sample.bed6";
  REQUIRE(std::filesystem::exists(bed_file));
  auto p = bed::Parser(bed_file, BED::BED3);
  auto records = p.parse_all();
  std::sort(records.begin(), records.end());
  CHECK(records[0].chrom == "chr7");
  CHECK(records[0].chrom_start == 127471196);
  CHECK(records[0].chrom_end == 127472363);
  CHECK(records[0].score == 0);
  CHECK(records[0].strand == '.');
}

TEST_CASE("BED Parser: BED6", "[parsers][BED][medium]") {
  const std::string bed_file = "../../test/data/unit_tests/sample.bed6";
  auto p = bed::Parser(bed_file);
  auto records = p.parse_all();
  compare_bed_records_with_file(records, bed_file);
}

TEST_CASE("BED Parser: BED9", "[parsers][BED][medium]") {
  std::string bed_file = "../../test/data/unit_tests/sample.bed9";
  auto p = bed::Parser(bed_file);
  auto records = p.parse_all();
  compare_bed_records_with_file(records, bed_file);
}

}  // namespace modle::test::bed