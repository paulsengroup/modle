#pragma once

#include <absl/types/span.h>  // for Span
#include <fmt/format.h>       // for print

#include <catch2/catch.hpp>  // for operator""_catch_sr, AssertionHandler, SourceLineInfo, Str...
#include <cstddef>           // IWYU pragma: keep for size_t
#include <cstdint>           // for uint64_t, int32_t, uint8_t
#include <filesystem>        // for operator/, path, create_directories, is_empty, remove, rem...
#include <memory>            // for allocator
#include <stdexcept>         // for runtime_error
#include <string_view>       // for string_view

#include "modle/common/smartdir.hpp"  // IWYU pragma: keep
#include "modle/contacts.hpp"         // for ContactMatrix
#include "modle/cooler.hpp"  // for Cooler, Cooler::READ_ONLY, Cooler::COOL, Cooler::WRITE_ONLY

namespace modle::test {
const extern SmartDir testdir;  // NOLINT
}
namespace modle::test::cooler {
using namespace modle::cooler;

inline const std::filesystem::path data_dir{"test/data/unit_tests"};  // NOLINT

TEST_CASE("cooler ctor", "[io][cooler][short]") {
  const auto test_file = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.1000kb.cool";

  {
    auto c = Cooler(test_file, Cooler::READ_ONLY);  // NOLINT
    CHECK(c.is_read_only());
  }

  H5::H5File f(test_file.string(), H5F_ACC_RDONLY);
  CHECK(Cooler::validate_file_format(f, Cooler::COOL));
  CHECK_THROWS_WITH(!Cooler::validate_file_format(f, Cooler::MCOOL, Cooler::READ_ONLY, 1'000'000),
                    Catch::Matchers::Contains("Expected format flavor MCOOL, found COOL"));
}

TEST_CASE("CMatrix to cooler", "[io][cooler][short]") {
  const auto test_file = testdir() / "cmatrix_to_cooler.cool";
  std::filesystem::create_directories(testdir());

  constexpr std::string_view chrom = "chr0";
  constexpr uint64_t start = 0;
  constexpr uint64_t end = 5'000'000;
  constexpr uint64_t bin_size = 1000;
  constexpr uint64_t nrows = 25;
  constexpr uint64_t ncols = end / bin_size;

  auto c = Cooler(test_file, Cooler::WRITE_ONLY, bin_size, chrom.size());

  ContactMatrix<int32_t> cmatrix(nrows, ncols, true);
  c.write_or_append_cmatrix_to_file(cmatrix, chrom, start, end, end);

  std::filesystem::remove(test_file);
}

TEST_CASE("Cooler to CMatrix", "[io][cooler][short]") {
  const auto test_file = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.1000kb.cool";

  constexpr std::string_view chrom = "chr7";
  constexpr uint64_t end = 159'138'663;
  constexpr uint64_t bin_size = 1'000'000;
  constexpr uint64_t nrows = 25;
  constexpr uint64_t ncols = (end / bin_size) + static_cast<uint64_t>(end % bin_size != 0);

  auto c = Cooler(test_file, Cooler::READ_ONLY);

  auto cmatrix = c.cooler_to_cmatrix(chrom, nrows, {0, -1}, true, false);
  CHECK(cmatrix.nrows() == nrows);
  CHECK(cmatrix.ncols() == ncols);
}

TEST_CASE("Cooler to CMatrix and CMatrix to Cooler", "[io][cooler][short]") {
  const auto test_file_in = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.1000kb.cool";
  const auto test_file_out = testdir() / "cmatrix_to_cooler.cool";

  auto c1 = Cooler(test_file_in, Cooler::READ_ONLY);

  constexpr std::string_view chrom = "chr1";
  constexpr uint64_t start = 0;
  constexpr uint64_t end = 249'250'621;
  constexpr uint64_t nrows = 25;
  const auto bin_size = c1.get_bin_size();
  const uint64_t ncols = (end / bin_size) + static_cast<uint64_t>(end % bin_size != 0);

  const auto cmatrix1 = c1.cooler_to_cmatrix(chrom, nrows, {start, end}, true, false);
  REQUIRE(cmatrix1.nrows() == nrows);
  REQUIRE(cmatrix1.ncols() == ncols);

  {
    auto c2 = Cooler(test_file_out, Cooler::WRITE_ONLY, bin_size, chrom.size());
    c2.write_or_append_cmatrix_to_file(cmatrix1, chrom, start, end, end);
  }
  auto c3 = Cooler(test_file_out, Cooler::READ_ONLY);

  const auto cmatrix2 = c3.cooler_to_cmatrix(chrom, nrows);
  const auto& v1 = cmatrix1.get_raw_count_vector();
  const auto& v2 = cmatrix2.get_raw_count_vector();
  REQUIRE(v1.size() == v2.size());

  size_t mismatches = 0;
  for (auto i = 0UL; i < v1.size(); ++i) {
    CHECK(v1[i] == v2[i]);
    mismatches += static_cast<uint8_t>(v1[i] != v2[i]);
  }

  CHECK(mismatches == 0);

  std::filesystem::remove(test_file_out);
}
/*
TEST_CASE("Cooler testing balanced matrix", "[io][cooler][short]") {
  // const auto test_file = data_dir / "4DNFI5RMAGF9_test.mcool";
  const auto test_file = data_dir / "4DNFIXP4QG5B_10kb.cool";

  // constexpr uint64_t start = 0;
  constexpr uint64_t end = 248'956'422;
  constexpr uint64_t bin_size = 10'000;
  constexpr uint64_t ncols = (end / bin_size) + (end % bin_size != 0);
  constexpr uint64_t nrows = ncols;

  cooler::Cooler c1(test_file, cooler::Cooler::READ_ONLY, bin_size);

  const auto balanced_cmatrix = c1.cooler_to_cmatrix("chr1", nrows, false, true);
  const auto raw_cmatrix = c1.cooler_to_cmatrix("chr1", nrows, false, false);

  REQUIRE(raw_cmatrix.nrows() == nrows);
  REQUIRE(raw_cmatrix.ncols() == ncols);
  REQUIRE(raw_cmatrix.nrows() == balanced_cmatrix.nrows());
  REQUIRE(raw_cmatrix.ncols() == balanced_cmatrix.ncols());

  std::vector<uint64_t> balanced_row_sum(ncols, 0);
  std::vector<uint64_t> raw_row_sum(ncols, 0);
  std::vector<uint64_t> balanced_col_sum(ncols, 0);
  std::vector<uint64_t> raw_col_sum(ncols, 0);

  for (auto i = 0UL; i < ncols; ++i) {
    for (auto j = i; j < (i + nrows) && j < ncols; ++j) {
      balanced_row_sum[j] += balanced_cmatrix.get(j, i);
      raw_row_sum[j] += raw_cmatrix.get(j, i);
      balanced_col_sum[i] += balanced_cmatrix.get(j, i);
      raw_col_sum[i] += raw_cmatrix.get(j, i);
    }
  }

  double avg_col_diff = 0;
  double avg_row_diff = 0;
  size_t nnz_cols = 0;
  size_t nnz_rows = 0;

  for (auto i = 0UL; i < nrows; ++i) {
    if (balanced_row_sum[i] != 0) {
      // fmt::print(stderr, "i={}; c={}\n", i, balanced_row_sum[i]);
    }
  }

  for (auto i = 0UL; i < raw_col_sum.size(); ++i) {
    if (raw_col_sum[i] != 0 || balanced_col_sum[i] != 0) {
      avg_col_diff +=
          static_cast<int64_t>(raw_col_sum[i]) - static_cast<int64_t>(balanced_col_sum[i]);
      nnz_rows++;
    }
  }
  avg_col_diff /= static_cast<double>(nnz_cols);
  avg_row_diff /= static_cast<double>(nnz_rows);

  fmt::print(
      stderr,
      FMT_STRING(
          "nrows={}\nncols={}\nbalanced_tot_contacts={}\nraw_tot_contacts={}\navg_col_diff={}\navg_"
          "row_diff={}\n"),
      nrows, ncols, balanced_cmatrix.get_tot_contacts(), raw_cmatrix.get_tot_contacts(),
      avg_col_diff, avg_row_diff);
}
 */

}  // namespace modle::test::cooler
