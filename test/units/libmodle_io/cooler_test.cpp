// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/cooler/cooler.hpp"

#include <absl/strings/str_split.h>  // for StrSplit, Splitter
#include <absl/types/span.h>         // for Span
#include <fmt/format.h>              // for format, make_format_args, vformat_to
#include <spdlog/logger.h>           // for logger
#include <spdlog/sinks/null_sink.h>  // for null_sink_mt
#include <spdlog/spdlog-inl.h>       // for set_default_logger

#include <algorithm>                        // for max, max_element, transform
#include <boost/filesystem/operations.hpp>  // for remove, create_directories
#include <boost/filesystem/path.hpp>        // for operator/, path
#include <catch2/catch.hpp>                 // for operator""_catch_sr, AssertionHandler, Source...
#include <exception>                        // for exception
#include <memory>                           // for make_shared, allocator_traits<>::value_type
#include <stdexcept>                        // for runtime_error
#include <string_view>                      // for string_view

#include "modle/common/common.hpp"                // for u64, u32, usize, i64, u8
#include "modle/common/smartdir.hpp"              // for SmartDir
#include "modle/common/utils.hpp"                 // for parse_numeric_or_throw
#include "modle/compressed_io/compressed_io.hpp"  // for Reader
#include "modle/contacts.hpp"                     // for ContactMatrix

namespace modle::test {
[[maybe_unused]] static const boost::filesystem::path& testdir(bool cleanup_on_exit = true) {
  static const SmartDir dir{cleanup_on_exit};
  return dir();
}
}  // namespace modle::test

namespace modle::test::cooler {
using namespace modle::cooler;
using Cooler = Cooler<>;
// using default_sink_t = spdlog::sinks::stderr_color_sink_mt;
using default_sink_t = spdlog::sinks::null_sink_mt;

inline const boost::filesystem::path data_dir{"test/data/unit_tests"};

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("cooler ctor", "[io][cooler][short]") {
  const auto test_file = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.1000kb.cool";
  spdlog::set_default_logger(
      std::make_shared<spdlog::logger>("main_logger", std::make_shared<default_sink_t>()));
  {
    auto c = Cooler(test_file, Cooler::IO_MODE::READ_ONLY);
    CHECK(c.is_read_only());
  }

  H5::H5File f(test_file.string(), H5F_ACC_RDONLY);
  CHECK(Cooler::validate_file_format(f, Cooler::FLAVOR::COOL));
  CHECK_THROWS_WITH(!Cooler::validate_file_format(f, Cooler::FLAVOR::MCOOL,
                                                  Cooler::IO_MODE::READ_ONLY, 1'000'000),
                    Catch::Matchers::Contains("Expected format flavor MCOOL, found COOL"));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix to cooler", "[io][cooler][short]") {
  const auto input_file = data_dir / "cmatrix_001.tsv.gz";
  const auto output_file = testdir() / "cmatrix_to_cooler.cool";
  boost::filesystem::create_directories(testdir());
  spdlog::set_default_logger(
      std::make_shared<spdlog::logger>("main_logger", std::make_shared<default_sink_t>()));

  constexpr std::string_view chrom = "chr0";
  const u64 start = 0;
  const u64 end = 300'000;
  const u64 bin_size = 100'000;
  const u64 nrows = 3;
  const u64 ncols = (end + bin_size - 1) / bin_size;

  ContactMatrix<> cmatrix1{};
  cmatrix1.unsafe_import_from_txt(input_file);

  Cooler(output_file, Cooler::IO_MODE::WRITE_ONLY, bin_size, chrom.size())
      .write_or_append_cmatrix_to_file(cmatrix1, chrom, start, end, end);

  const auto cmatrix2 =
      Cooler(output_file, Cooler::IO_MODE::READ_ONLY).cooler_to_cmatrix(chrom, nrows);

  for (usize i = 0; i < ncols; ++i) {
    for (usize j = 0; j < ncols; ++j) {
      CHECK(cmatrix1.get(i, j) == cmatrix2.get(i, j));
    }
  }

  boost::filesystem::remove(output_file);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix to cooler - multiple chromosomes", "[io][cooler][short]") {
  const auto input_file = data_dir / "cmatrix_001.tsv.gz";
  const auto output_file = testdir() / "cmatrix_to_cooler.cool";
  boost::filesystem::create_directories(testdir());
  spdlog::set_default_logger(
      std::make_shared<spdlog::logger>("main_logger", std::make_shared<default_sink_t>()));

  constexpr std::string_view chrom1 = "chr1";
  constexpr std::string_view chrom2 = "chr2";
  const u64 start = 0;
  const u64 end = 300'000;
  const u64 bin_size = 100'000;
  const u64 nrows = 3;
  const u64 ncols = (end + bin_size - 1) / bin_size;

  ContactMatrix<> cmatrix1{};
  cmatrix1.unsafe_import_from_txt(input_file);

  {
    auto c = Cooler(output_file, Cooler::IO_MODE::WRITE_ONLY, bin_size, chrom1.size());
    c.write_or_append_cmatrix_to_file(cmatrix1, chrom1, start, end, end);
    c.write_or_append_cmatrix_to_file(cmatrix1, chrom2, start, end, end);
  }

  const auto cmatrix2 =
      Cooler(output_file, Cooler::IO_MODE::READ_ONLY).cooler_to_cmatrix(chrom1, nrows);
  const auto cmatrix3 =
      Cooler(output_file, Cooler::IO_MODE::READ_ONLY).cooler_to_cmatrix(chrom2, nrows);

  for (usize i = 0; i < ncols; ++i) {
    for (usize j = 0; j < ncols; ++j) {
      CHECK(cmatrix1.get(i, j) == cmatrix2.get(i, j));
      CHECK(cmatrix1.get(i, j) == cmatrix3.get(i, j));
    }
  }

  boost::filesystem::remove(output_file);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler to CMatrix", "[io][cooler][short]") {
  const auto test_file = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.1000kb.cool";
  const auto reference_file =
      data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.1000kb_chr7.tsv.gz";
  spdlog::set_default_logger(
      std::make_shared<spdlog::logger>("main_logger", std::make_shared<default_sink_t>()));

  constexpr std::string_view chrom = "chr7";
  const u64 end = 159'138'663;
  const u64 bin_size = 1'000'000;
  const u64 ncols = (end + bin_size - 1) / bin_size;
  const auto nrows = ncols;

  auto c = Cooler(test_file, Cooler::IO_MODE::READ_ONLY);

  auto cmatrix = c.cooler_to_cmatrix(chrom, nrows, {0, -1}, true, false);
  CHECK(cmatrix.nrows() == nrows);
  CHECK(cmatrix.ncols() == ncols);

  modle::compressed_io::Reader reader(reference_file);
  std::vector<std::string_view> toks;
  std::vector<u32> reference_row(ncols, 0);

  for (usize i = 0; i < cmatrix.ncols(); ++i) {
    const auto line = reader.getline();
    toks = absl::StrSplit(line, '\t');
    REQUIRE(toks.size() == ncols);
    std::transform(toks.begin(), toks.end(), reference_row.begin(),
                   [](const auto tok) { return utils::parse_numeric_or_throw<u32>(tok); });
    for (usize j = 0; j < cmatrix.ncols(); ++j) {
      CHECK(cmatrix.get(i, j) == reference_row[j]);
    }
  }

  CHECK(cmatrix.get_tot_contacts() == 5753412);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler to CMatrix and CMatrix to Cooler", "[io][cooler][short]") {
  const auto test_file_in = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.1000kb.cool";
  const auto test_file_out = testdir() / "cmatrix_to_cooler.cool";
  spdlog::set_default_logger(
      std::make_shared<spdlog::logger>("main_logger", std::make_shared<default_sink_t>()));

  auto c1 = Cooler(test_file_in, Cooler::IO_MODE::READ_ONLY);

  constexpr std::string_view chrom = "chr1";
  const u64 start = 0;
  const u64 end = 249'250'621;
  const u64 nrows = 25;
  const auto bin_size = c1.get_bin_size();
  const u64 ncols = (end + bin_size - 1) / bin_size;

  const auto cmatrix1 = c1.cooler_to_cmatrix(chrom, nrows, {start, end}, true, false);
  REQUIRE(cmatrix1.nrows() == nrows);
  REQUIRE(cmatrix1.ncols() == ncols);

  {
    auto c2 = Cooler(test_file_out, Cooler::IO_MODE::WRITE_ONLY, bin_size, chrom.size());
    c2.write_or_append_cmatrix_to_file(cmatrix1, chrom, start, end, end);
  }
  auto c3 = Cooler(test_file_out, Cooler::IO_MODE::READ_ONLY);

  const auto cmatrix2 = c3.cooler_to_cmatrix(chrom, nrows);
  const auto v1 = cmatrix1.get_raw_count_vector();
  const auto v2 = cmatrix2.get_raw_count_vector();
  REQUIRE(v1.size() == v2.size());

  usize mismatches = 0;
  for (usize i = 0; i < v1.size(); ++i) {
    CHECK(v1[i] == v2[i]);
    mismatches += static_cast<u8>(v1[i] != v2[i]);
  }

  CHECK(mismatches == 0);

  boost::filesystem::remove(test_file_out);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler to CMatrix and CMatrix to Cooler - all chromosomes", "[io][cooler][long]") {
  const auto test_file_in = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.1000kb.cool";
  const auto test_file_out = testdir() / "cmatrix_to_cooler_all_chromosomes.cool";
  spdlog::set_default_logger(
      std::make_shared<spdlog::logger>("main_logger", std::make_shared<default_sink_t>()));

  auto c1 = Cooler(test_file_in, Cooler::IO_MODE::READ_ONLY);

  const auto bin_size = c1.get_bin_size();
  const usize nrows = 25;

  const auto chrom_names = c1.get_chrom_names();
  const auto chrom_sizes = c1.get_chrom_sizes();

  std::vector<ContactMatrix<>> matrices;
  const auto max_chrom_name_size =
      std::max_element(chrom_names.begin(), chrom_names.end(), [](const auto& a, const auto& b) {
        return a.size() < b.size();
      })->size();

  {
    REQUIRE(chrom_names.size() == chrom_sizes.size());
    auto c2 = Cooler(test_file_out, Cooler::IO_MODE::WRITE_ONLY, bin_size, max_chrom_name_size + 1);
    for (usize i = 0; i < chrom_names.size(); ++i) {
      matrices.emplace_back(c1.cooler_to_cmatrix(chrom_names[i], nrows, {0, -1}, true, false));
      c2.write_or_append_cmatrix_to_file(matrices.back(), chrom_names[i], i64(0), chrom_sizes[i],
                                         chrom_sizes[i]);
    }
  }

  auto c3 = Cooler(test_file_out, Cooler::IO_MODE::READ_ONLY);
  for (usize i = 0; i < chrom_names.size(); ++i) {
    const auto matrix = c3.cooler_to_cmatrix(chrom_names[i], nrows, {0, -1}, true, false);
    const auto& reference = matrices[i];
    REQUIRE(matrix.ncols() == reference.ncols());
    REQUIRE(matrix.nrows() == reference.nrows());
    for (usize j = 0; j < reference.ncols(); ++j) {
      for (usize k = 0; k < reference.ncols(); ++k) {
        CHECK(matrix.get(j, k) == reference.get(j, k));
      }
    }
  }
  boost::filesystem::remove(test_file_out);
}

/*
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler testing balanced matrix", "[io][cooler][short]") {
  // const auto test_file = data_dir / "4DNFI5RMAGF9_test.mcool";
  const auto test_file = data_dir / "4DNFIXP4QG5B_10kb.cool";

  // const u64 start = 0;
  const u64 end = 248'956'422;
  const u64 bin_size = 10'000;
  const u64 ncols = (end / bin_size) + (end % bin_size != 0);
  const u64 nrows = ncols;

  Cooler::FLAVOR::COOLer c1(test_file, Cooler::FLAVOR::COOLer<>::IO_MODE::READ_ONLY, bin_size);

  const auto balanced_cmatrix = c1.cooler_to_cmatrix("chr1", nrows, false, true);
  const auto raw_cmatrix = c1.cooler_to_cmatrix("chr1", nrows, false, false);

  REQUIRE(raw_cmatrix.nrows() == nrows);
  REQUIRE(raw_cmatrix.ncols() == ncols);
  REQUIRE(raw_cmatrix.nrows() == balanced_cmatrix.nrows());
  REQUIRE(raw_cmatrix.ncols() == balanced_cmatrix.ncols());

  std::vector<u64> balanced_row_sum(ncols, 0);
  std::vector<u64> raw_row_sum(ncols, 0);
  std::vector<u64> balanced_col_sum(ncols, 0);
  std::vector<u64> raw_col_sum(ncols, 0);

  for (usize  i = 0; i < ncols; ++i) {
    for (auto j = i; j < (i + nrows) && j < ncols; ++j) {
      balanced_row_sum[j] += balanced_cmatrix.get(j, i);
      raw_row_sum[j] += raw_cmatrix.get(j, i);
      balanced_col_sum[i] += balanced_cmatrix.get(j, i);
      raw_col_sum[i] += raw_cmatrix.get(j, i);
    }
  }

  double avg_col_diff = 0;
  double avg_row_diff = 0;
  usize nnz_cols = 0;
  usize nnz_rows = 0;

  for (usize  i = 0; i < nrows; ++i) {
    if (balanced_row_sum[i] != 0) {
      // fmt::print(stderr, "i={}; c={}\n", i, balanced_row_sum[i]);
    }
  }

  for (usize  i = 0; i < raw_col_sum.size(); ++i) {
    if (raw_col_sum[i] != 0 || balanced_col_sum[i] != 0) {
      avg_col_diff +=
          static_cast<i64>(raw_col_sum[i]) - static_cast<i64>(balanced_col_sum[i]);
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler to CMatrix subrange", "[io][cooler][short]") {
  const auto test_file_in = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.1000kb.cool";
  spdlog::set_default_logger(
      std::make_shared<spdlog::logger>("main_logger", std::make_shared<default_sink_t>()));

  auto c = Cooler(test_file_in, Cooler::IO_MODE::READ_ONLY);

  constexpr std::string_view chrom = "chr1";
  const u64 start1 = 0;
  const u64 end1 = 249'250'621;
  const u64 start2 = 50'000'000;
  const u64 end2 = 75'000'000;
  const u64 nrows = 25;
  const auto bin_size = c.get_bin_size();

  const auto cmatrix2 = c.cooler_to_cmatrix(chrom, nrows, {start2, end2}, true, false);

  const auto cmatrix1 = c.cooler_to_cmatrix(chrom, nrows, {start1, end1}, true, false);
  REQUIRE(cmatrix2.nrows() == 25);
  REQUIRE(cmatrix2.ncols() == 25);

  for (auto p1 = start2; p1 < end2; p1 += bin_size) {
    for (auto p2 = start2; p2 < end2; p2 += bin_size) {
      const auto b11 = p1 / bin_size;
      const auto b12 = p2 / bin_size;

      const auto b21 = (p1 - start2) / bin_size;
      const auto b22 = (p2 - start2) / bin_size;
      CHECK(cmatrix1.get(b11, b12) == cmatrix2.get(b21, b22));
    }
  }
}

}  // namespace modle::test::cooler
