// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/bigwig/bigwig.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>  // for path
#include <numeric>
#include <string>  // for string, basic_string, operator==, char_traits, stoull
#include <vector>  // for vector

#include "modle/common/common.hpp"              // for usize
#include "modle/test/self_deleting_folder.hpp"  // for SelfDeletingFolder

namespace modle::test {
inline const SelfDeletingFolder testdir{true};  // NOLINT(cert-err58-cpp)
}  // namespace modle::test

/*
# Test file was generated with the following Python code:

import numpy as np
import pyBigWig
import os

chroms = {"chr1": 1_000_000, "chr2": 500_000}
bin_size = 100_000

os.remove("/tmp/test.bw")
with pyBigWig.open("/tmp/test.bw", "w") as f:
    f.addHeader([(chrom, size) for chrom, size in chroms.items()], maxZooms=1)
    for chrom, size in chroms.items():
        starts = np.arange(0, size, bin_size)
        values = np.arange(0, len(starts), dtype=float)

        f.addEntries(chrom, starts.tolist(), values=values.tolist(), span=bin_size, step=bin_size)
*/

namespace modle::bigwig::test {

constexpr auto& testdir = modle::test::testdir;

[[maybe_unused]] static const std::filesystem::path& data_dir() {
  static const std::filesystem::path data_dir{"test/data/unit_tests"};
  return data_dir;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("bigwig::Reader", "[bigwig][io][short]") {
  const auto path = data_dir() / "bigwig" / "test.bw";

  const bp_t bin_size = 100'000;
  io::bigwig::Reader r{path};
  REQUIRE(!!r);

  SECTION("chromosomes") {
    const auto& chroms = r.chromosomes();
    CHECK(chroms.size() == 2);
    CHECK(chroms.contains("chr1"));
    CHECK(chroms.contains("chr2"));
  }

  SECTION("read values") {
    auto values = r.read_values("chr1", 0, 100);
    CHECK(values.size() == 100);
    for (const auto& v : values) {
      CHECK(v == 0);
    }
    values = r.read_values("chr2", bin_size, bin_size + 10);
    CHECK(values.size() == 10);
    for (const auto& v : values) {
      CHECK(v == 1);
    }
  }

  SECTION("read intervals") {
    const auto intervals = r.get_intervals("chr1", 0, bin_size * 10);
    REQUIRE(!!intervals);
    CHECK(intervals->l == 10);

    for (u32 i = 0; i < intervals->l; ++i) {
      CHECK(intervals->start[i] == bin_size * i);           // NOLINT
      CHECK(intervals->end[i] == bin_size * (i + 1));       // NOLINT
      CHECK(intervals->value[i] == static_cast<float>(i));  // NOLINT
    }
  }

  SECTION("read stats") {
    const auto values = r.stats<mean>("chr1", 0, 4 * bin_size, 2 * bin_size);
    REQUIRE(values.size() == 2);

    CHECK_THAT(values[0], Catch::Matchers::WithinRel(0.5));
    CHECK_THAT(values[1], Catch::Matchers::WithinRel(2.5));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("bigwig::Writer", "[bigwig][io][short]") {
  const auto path = testdir() / "test.bw";

  const std::vector<std::string> chrom_names{"chr1", "chr2"};
  const std::vector<u32> chrom_sizes{1000, 500};
  const u32 bin_size = 10;
  {
    io::bigwig::Writer w{path};
    REQUIRE(!!w);
    w.write_chromosomes(chrom_names, chrom_sizes);

    for (usize i = 0; i < chrom_names.size(); ++i) {
      std::vector<float> values(chrom_sizes[i] / bin_size);
      std::iota(values.begin(), values.end(), 0.0F);

      w.write_range(chrom_names[i], absl::MakeConstSpan(values), bin_size, bin_size);
    }
  }

  io::bigwig::Reader r{path};
  const auto& chroms = r.chromosomes();
  CHECK(chroms.size() == 2);
  CHECK(chroms.contains("chr1"));
  CHECK(chroms.contains("chr2"));

  for (const auto& chrom : chroms) {
    const auto intervals = r.get_intervals(chrom.first, 0, chrom.second);
    REQUIRE(intervals->l == ((chrom.second + bin_size - 1) / bin_size));

    for (u32 i = 0; i < intervals->l; ++i) {
      CHECK(intervals->start[i] == bin_size * i);           // NOLINT
      CHECK(intervals->end[i] == bin_size * (i + 1));       // NOLINT
      CHECK(intervals->value[i] == static_cast<float>(i));  // NOLINT
    }
  }
}

}  // namespace modle::bigwig::test
