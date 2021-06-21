#pragma once

#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr, SourceLineInfo
#include <string>
#include <vector>

#include "modle/bed.hpp"

namespace modle::test::bed {
using namespace modle::bed;

TEST_CASE("BED Tree simple", "[BED][short]") {
  const std::string all_intervals = "test/data/unit_tests/H1_ctcf_all_chroms.bed";

  Parser p(all_intervals, bed::BED::BED3);

  const auto records = p.parse_all();
  p.reset();
  const auto intervals = p.parse_all_in_interval_tree();

  std::vector<const BED*> buff;
  for (const auto& record : records) {
    CHECK(intervals.contains(record.chrom));
    CHECK(intervals.count_overlaps(record) == 1);
    intervals.find_overlaps(record, buff);
    CHECK(buff.size() == 1);
  }
}

}  // namespace modle::test::bed
