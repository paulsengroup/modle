// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/strings/str_split.h>  // for StrSplit, Splitter
#include <absl/types/span.h>         // for Span
#include <fmt/format.h>              // for format

#include <algorithm>  // for max, min
#include <catch2/catch_test_macros.hpp>
#include <filesystem>   // for path
#include <string>       // for basic_string, operator==, string
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "modle/bed/bed.hpp"        // for BED, BED_tree, BED_tree::contains, BED_tree::count_...
#include "modle/common/common.hpp"  // for usize, u8
#include "modle/common/numeric_utils.hpp"         // for parse_numeric_or_throw
#include "modle/compressed_io/compressed_io.hpp"  // for Reader

namespace modle::test::bed {
using namespace modle::bed;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BED Tree simple", "[BED][io][long]") {
  const std::filesystem::path all_intervals = "test/data/unit_tests/H1_hg38_CTCFs_filtered.bed.gz";

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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BED Tree multiple overlaps", "[BED][io][long]") {
  const std::filesystem::path all_intervals = "test/data/unit_tests/H1_ctcf_all_chroms.bed.gz";
  const std::filesystem::path counts_per_interval =
      "test/data/unit_tests/H1_ctcf_all_chroms_per_interval.tsv.gz";

  const BED_tree<> intervals{all_intervals, bed::BED::BED3};

  compressed_io::Reader r(counts_per_interval);
  std::string buff;

  std::vector<std::string_view> toks;
  while (r.getline(buff)) {
    toks = absl::StrSplit(buff, '\t');
    usize num_expected_overlaps, start, end;
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
