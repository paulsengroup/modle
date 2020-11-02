#include "gtest/gtest.h"
#include "modle/parsers.hpp"
#include "absl/strings/str_split.h"

namespace modle {

void compare_bed_records_with_file(std::vector<BED> records, const std::string &bed_file) {
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

  ASSERT_EQ(records.size(), lines.size());
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
    EXPECT_EQ(records[i].to_string(), lines[i]);
  }
}

TEST(modle_test_suite, bed_parser_simple_test) {
  std::string bed_file = "data/sample.bed6";
  auto p = BEDParser(bed_file);
  auto records = p.parse_all();
  EXPECT_EQ(records.size(), 9);
  std::sort(records.begin(), records.end());
  EXPECT_STREQ(records[0].chrom.c_str(), "chr7");
  EXPECT_EQ(records[0].chrom_start, 127471196);
  EXPECT_EQ(records[0].chrom_end, 127472363);
  EXPECT_DOUBLE_EQ(records[0].score, 0);
  EXPECT_EQ(records[0].strand, '+');
}

TEST(modle_test_suite, bed_parser_simple_test_bed3) {
  std::string bed_file = "data/sample.bed6";
  auto p = BEDParser(bed_file, BED::BED3);
  auto records = p.parse_all();
  std::sort(records.begin(), records.end());
  EXPECT_STREQ(records[0].chrom.c_str(), "chr7");
  EXPECT_EQ(records[0].chrom_start, 127471196);
  EXPECT_EQ(records[0].chrom_end, 127472363);
  EXPECT_DOUBLE_EQ(records[0].score, 0);
  EXPECT_EQ(records[0].strand, '.');
}

TEST(modle_test_suite, bed_parser_short_test) {
  std::string bed_file = "data/sample.bed6";
  auto p = BEDParser(bed_file);
  auto records = p.parse_all();
  compare_bed_records_with_file(records, bed_file);
}

TEST(modle_test_suite, bed_parser_medium_test) {
  std::string bed_file = "data/GM12878_CTCF_orientation.bed";
  auto p = BEDParser(bed_file);
  auto records = p.parse_all();
  compare_bed_records_with_file(records, bed_file);
}

}  // namespace modle