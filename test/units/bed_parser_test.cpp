#include "gtest/gtest.h"
#include "modle/parsers.hpp"

void compare_bed_records_with_file(std::vector<modle::BED> records, const std::string& bed_file) {
  auto fp = std::ifstream("data/sample.bed6");
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
  std::sort(lines.begin(), lines.end());

  for (auto i = 0UL; i < records.size(); ++i) {
    EXPECT_EQ(records[i].to_string(), lines[i]);
  }
}

TEST(modle_test_suite, bed_parser_simple_test) {
  auto p = modle::BEDParser("data/sample.bed6");
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
  auto p = modle::BEDParser("data/sample.bed6", modle::BED::BED3);
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
  auto p = modle::BEDParser(bed_file);
  auto records = p.parse_all();
  compare_bed_records_with_file(records, bed_file);
}