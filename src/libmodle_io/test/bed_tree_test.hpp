#pragma once

#include <absl/strings/str_split.h>

#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr, SourceLineInfo
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include "modle/bed.hpp"
#include "modle/compressed_io.hpp"

namespace modle::test::bed {
using namespace modle::bed;

TEST_CASE("BED Tree simple", "[BED][io][long]") {
  const std::string all_intervals = "test/data/unit_tests/H1_ctcf_all_chroms.bed.gz";

  Parser p(all_intervals, bed::BED::BED3);

  const auto records = p.parse_all();
  p.reset();
  const auto intervals = p.parse_all_in_interval_tree();

  for (const auto& record : records) {
    CHECK(intervals.contains(record.chrom));
    CHECK(intervals.count_overlaps(record) == 1);
    const auto overlaps = intervals.find_overlaps(record);
    CHECK(overlaps.size() == 1);
    if (!overlaps.empty()) {
      CHECK(overlaps.front() == record);
    }
  }
}

TEST_CASE("BED Tree multiple overlaps", "[BED][io][long]") {
  const std::string all_intervals = "test/data/unit_tests/H1_ctcf_all_chroms.bed.gz";
  const std::string counts_per_interval =
      "test/data/unit_tests/H1_ctcf_all_chroms_per_interval.tsv.gz";

  Parser p(all_intervals, bed::BED::BED3);
  const auto intervals = p.parse_all_in_interval_tree();

  compressed_io::Reader r(counts_per_interval);
  std::string buff;

  std::vector<std::string_view> toks;
  while (r.getline(buff)) {
    toks = absl::StrSplit(buff, '\t');
    size_t num_expected_overlaps, start, end;  // NOLINT
    const auto chrom = std::string{toks.front()};
    utils::parse_numeric_or_throw(toks[1], start);
    utils::parse_numeric_or_throw(toks[2], end);
    utils::parse_numeric_or_throw(toks[3], num_expected_overlaps);
    CHECK(intervals.contains(chrom));
    CHECK(intervals.count_overlaps(chrom, start, end) == num_expected_overlaps);

    for (const auto& record : intervals.find_overlaps(chrom, start, end)) {
      CHECK(record.chrom == chrom);
      CHECK(record.chrom_start < end);
      CHECK(record.chrom_end >= start);
    }
  }
}

}  // namespace modle::test::bed
