#pragma once

#include <catch2/catch.hpp>
#include <filesystem>
#include <string_view>

#include "modle/contacts.hpp"
#include "modle/cooler.hpp"

namespace modle::test::cooler {
using namespace modle::cooler;

inline const std::filesystem::path test_dir{"/tmp/modle/unit_tests"};  // NOLINT
inline const std::filesystem::path data_dir{"test/data/cooler"};       // NOLINT

TEST_CASE("cooler ctor", "[io][cooler][short]") {
  const auto test_file = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.100kb.cool";

  {
    auto c = Cooler(test_file.string(), Cooler::READ_ONLY, 100'000);  // NOLINT
    CHECK(c.is_read_only());
  }

  H5::H5File f(test_file.string(), H5F_ACC_RDONLY);
  CHECK(Cooler::detect_file_flavor(f) == Cooler::COOL);
  CHECK(Cooler::validate_file_format(f, Cooler::COOL));
  CHECK_THROWS_WITH(!Cooler::validate_file_format(f, Cooler::MCOOL, Cooler::READ_ONLY, 1000),
                    Catch::Matchers::Contains("Expected format flavor MCOOL, found COOL"));
}

TEST_CASE("CMatrix to cooler", "[io][cooler][short]") {
  const auto test_file = test_dir / "cmatrix_to_cooler.cool";
  std::filesystem::create_directories(test_dir);

  constexpr uint64_t start = 0;
  constexpr uint64_t end = 5'000'000;
  constexpr uint64_t bin_size = 1000;
  constexpr uint64_t nrows = 25;
  constexpr uint64_t ncols = end / bin_size;

  auto c = Cooler(test_file.string(), Cooler::WRITE_ONLY, bin_size);

  ContactMatrix<int32_t> cmatrix(nrows, ncols, true);
  c.write_cmatrix_to_file(cmatrix, "chr0", start, end, end);

  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}

TEST_CASE("Cooler to CMatrix", "[io][cooler][short]") {
  const auto test_file = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.100kb.cool";

  // constexpr uint64_t bin_size = 100'000;
  constexpr uint64_t nrows = 25;
  constexpr uint64_t ncols = 1'591 + 1;

  auto c = Cooler(test_file.string(), Cooler::READ_ONLY);

  auto cmatrix = c.cooler_to_cmatrix("chr7", nrows);
  CHECK(cmatrix.n_rows() == nrows);
  CHECK(cmatrix.n_cols() == ncols);
}

TEST_CASE("Cooler to CMatrix and CMatrix to Cooler", "[io][cooler][short]") {
  const auto test_file_in = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.100kb.cool";
  const auto test_file_out = test_dir / "cmatrix_to_cooler.cool";
  std::filesystem::create_directories(test_dir);

  constexpr uint64_t start = 0;
  constexpr uint64_t end = 249'250'621;
  constexpr uint64_t bin_size = 100'000;
  constexpr uint64_t nrows = 25;
  constexpr uint64_t ncols = (end / bin_size) + (end % bin_size != 0);

  auto c1 = Cooler(test_file_in.string(), Cooler::READ_ONLY);

  const auto cmatrix1 = c1.cooler_to_cmatrix("chr1", nrows);
  REQUIRE(cmatrix1.n_rows() == nrows);
  REQUIRE(cmatrix1.n_cols() == ncols);

  auto c2 = Cooler(test_file_out.string(), Cooler::WRITE_ONLY, bin_size);
  c2.write_cmatrix_to_file(cmatrix1, "chr1", start, end, end);
  const auto cmatrix2 = c2.cooler_to_cmatrix("chr1", nrows);
  const auto& v1 = cmatrix1.get_raw_count_vector();
  const auto& v2 = cmatrix2.get_raw_count_vector();
  REQUIRE(v1.size() == v2.size());

  std::size_t mismatches = 0;
  for (auto i = 0UL; i < v1.size(); ++i) {
    CHECK(v1[i] == v2[i]);
    mismatches += static_cast<uint8_t>(v1[i] != v2[i]);
  }

  CHECK(mismatches == 0);

  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}

}  // namespace modle::test::cooler