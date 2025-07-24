// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "modle/bed/bed.hpp"
#include "modle/common/common.hpp"
#include "modle/common/string_utils.hpp"
#include "modle/compressed_io/compressed_io.hpp"
#include "modle/test/tmpdir.hpp"

namespace modle::test {
inline const TmpDir testdir{true};  // NOLINT(cert-err58-cpp)
}  // namespace modle::test

namespace modle::bed::test {

constexpr auto &testdir = modle::test::testdir;

[[maybe_unused]] static const std::filesystem::path &data_dir() {
  static const std::filesystem::path data_dir{"test/data/unit_tests"};
  return data_dir;
}

static void compare_bed_records_with_file(std::vector<BED> records, const std::string &bed_file) {
  compressed_io::Reader r(bed_file);
  std::vector<std::string> lines;
  lines.reserve(records.size());
  std::string buff;
  while (r.getline(buff)) {
    if (buff.front() == '#' || str_contains(buff, "track") || str_contains(buff, "browser")) {
      continue;
    }
    lines.push_back(buff);
  }

  REQUIRE(records.size() == lines.size());
  std::sort(records.begin(), records.end());
  std::sort(lines.begin(), lines.end(), [&](const std::string_view &a, const std::string_view &b) {
    const auto toksa = str_split(a, '\t');
    const auto toksb = str_split(b, '\t');
    const auto &chra = toksa[0];
    const auto &chrb = toksb[0];
    const auto &starta = toksa[1];
    const auto &startb = toksb[1];
    const auto &enda = toksa[2];
    const auto &endb = toksb[2];
    if (chra != chrb) {
      return chra < chrb;
    }
    if (starta != startb) {
      return std::stoull(starta.data(), nullptr) < std::stoull(startb.data(), nullptr);
    }
    return std::stoull(enda.data(), nullptr) < std::stoull(endb.data(), nullptr);
  });

  for (usize i = 0; i < records.size(); ++i) {
    CHECK(fmt::to_string(records[i]) == lines[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BED: strip quotes", "[parsers][BED][io][short]") {
  SECTION("valid") {
    constexpr std::string_view line{
        "chr1\t"
        "0\t"
        "10\t"
        "\"name\"\t"
        "0.0\t"
        "\"+\"\t"
        "0\t"
        "1\t"
        "\"0,0,0\""};

    const bed::BED record(line);
    CHECK(record.chrom == "chr1");
    CHECK(record.chrom_start == 0);
    CHECK(record.chrom_end == 10);
    CHECK(record.name == "name");
    CHECK(record.score == 0.0);
    CHECK(record.strand == '+');
    CHECK(record.thick_start == 0);
    CHECK(record.thick_end == 1);
    CHECK(*record.rgb == RGB{});

    CHECK(bed::BED("\"chr1\t0\t1").chrom == "\"chr1");
  }

  SECTION("invalid") {
    CHECK_THROWS(bed::BED("chr1\t\"0\"\t1"));
    CHECK_THROWS(bed::BED("chr1\t0\t1\t.\t\"0.0\""));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BED Parser CRLF", "[parsers][BED][io][short]") {
  const auto bed_file = testdir() / "crlf.bed";

  const usize num_records = 3;

  {
    compressed_io::Writer w(bed_file, compressed_io::Writer::NONE);
    for (usize i = 0; i < num_records; ++i) {
      w.write(fmt::format("chr{}\t0\t1\r\n", i));
    }
  }

  const auto records = bed::Parser(bed_file).parse_all();
  REQUIRE(records.size() == num_records);
  for (usize i = 0; i < num_records; ++i) {
    CHECK(bed::BED(fmt::format("chr{}\t0\t1", i)) == records[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BED Parser simple", "[parsers][BED][io][short]") {
  const auto bed_file = data_dir() / "genomic_intervals" / "intervals.bed6.xz";
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BED Parser simple: BED6 -> BED3", "[parsers][BED][io][short]") {
  const auto bed_file = data_dir() / "genomic_intervals" / "intervals.bed6.xz";
  auto p = bed::Parser(bed_file, BED::BED3);
  auto records = p.parse_all();
  std::sort(records.begin(), records.end());
  CHECK(records[0].chrom == "chr7");
  CHECK(records[0].chrom_start == 127471196);
  CHECK(records[0].chrom_end == 127472363);
  CHECK(records[0].score == 0);
  CHECK(records[0].strand == '.');
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BED Parser: BED6", "[parsers][BED][io][medium]") {
  const auto bed_file = data_dir() / "genomic_intervals" / "intervals.bed6.xz";
  auto p = bed::Parser(bed_file);
  auto records = p.parse_all();
  REQUIRE(records.size() == 9);
  compare_bed_records_with_file(records, bed_file);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BED Parser: BED9", "[parsers][BED][io][long]") {
  const auto bed_file = data_dir() / "genomic_intervals" / "intervals.bed9";
  auto p = bed::Parser(bed_file);
  auto records = p.parse_all();
  REQUIRE(records.size() == 62580);
  compare_bed_records_with_file(records, bed_file);
}

}  // namespace modle::bed::test
