// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "modle/bed/bed.hpp"
#include "modle/common/common.hpp"
#include "modle/common/numeric_utils.hpp"
#include "modle/common/string_utils.hpp"
#include "modle/compressed_io/compressed_io.hpp"

namespace modle::bed::test {

[[maybe_unused]] static const std::filesystem::path& data_dir() {
  static const std::filesystem::path data_dir{"test/data/unit_tests"};
  return data_dir;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BED Tree simple", "[BED][io][long]") {
  const auto all_intervals =
      data_dir() / "genomic_intervals" / "H1_ctcf_binding_sites_filtered.bed.xz";

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
  const auto parent = data_dir() / "genomic_intervals";
  const auto all_intervals = parent / "H1_ctcf_binding_sites.bed.xz";
  const auto counts_per_interval = parent / "H1_ctcf_binding_sites_counts_per_interval.bed.xz";

  const BED_tree<> intervals{all_intervals, bed::BED::BED3};

  compressed_io::Reader r(counts_per_interval);
  std::string buff;

  std::vector<std::string_view> toks;
  while (r.getline(buff)) {
    toks = str_split(buff, '\t');
    usize num_expected_overlaps{};
    usize start{};
    usize end{};
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

}  // namespace modle::bed::test
