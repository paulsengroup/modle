#pragma once

#include <absl/strings/str_split.h>
#include <fmt/format.h>

#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <catch2/catch.hpp>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

#include "modle/contacts.hpp"
#include "modle/utils.hpp"

namespace modle::test::cmatrix {

TEST_CASE("CMatrix simple", "[cmatrix][short]") {
  ContactMatrix<uint32_t> c(10, 100);  // NOLINT
  CHECK(c.get(0, 0) == 0);
  c.increment(0, 0);
  CHECK(c.get(0, 0) == 1);
  c.increment(0, 0);
  CHECK(c.get(0, 0) == 2);
  c.subtract(0, 0, 2);
  CHECK(c.get(0, 0) == 0);
}

[[nodiscard]] inline std::vector<std::vector<uint32_t>> load_matrix_from_file(
    const std::string& path_to_file, const std::string& sep = "\t") {
  std::vector<std::vector<uint32_t>> m;
  std::ifstream f(path_to_file);
  std::string line;
  std::string buff;
  uint32_t n;  // NOLINT
  while (std::getline(f, line)) {
    std::vector<uint32_t> v;

    for (const auto& tok : absl::StrSplit(line, sep)) {
      modle::utils::parse_numeric_or_throw(tok, n);
      v.push_back(n);
    }
    m.emplace_back(std::move(v));
    if (!f && !f.eof()) {
      throw fmt::system_error(errno, "Unable to open file '{}'", path_to_file);
    }
  }
  return m;
}

TEST_CASE("CMatrix 10x200", "[cmatrix][medium]") {
  const std::string file_path = "../../test/data/symm_matrix_200_10.tsv";
  REQUIRE(std::filesystem::exists(file_path));
  const auto m1 = load_matrix_from_file(file_path);
  ContactMatrix<uint32_t> m2(10, 200);  // NOLINT
  for (auto i = 0UL; i < m1.size(); ++i) {
    for (auto j = 0UL; j < m1[i].size(); ++j) {
      if (m1[i][j] != 0 && j >= i) {
        m2.set(i, j, m1[i][j]);
      }
    }
  }

  const auto m3 = m2.generate_symmetric_matrix();
  for (auto i = 0UL; i < m1.size(); ++i) {
    for (auto j = 0UL; j < m1[0].size(); ++j) {
      CHECK(m1[i][j] == m3[i][j]);
    }
  }
}

TEST_CASE("CMatrix Mask", "[cmatrix][short]") {
  ContactMatrix<uint32_t> m(10, 20);  // NOLINT
  auto mask = m.generate_mask_for_bins_without_contacts();
  REQUIRE(mask.size() == m.ncols());
  CHECK(mask.none());  // Matrix is full of zeros: bitmask should also be all zeros

  for (auto i = 0UL; i < m.ncols(); ++i) {
    for (auto j = i; j < m.ncols(); ++j) {
      // Even row/columns are all zeros: Bitmask should have zeros at pos corresponding to even idx
      m.set(i, j, !(i % 2 == 0 || j % 2 == 0));
    }
  }
  // m.print(true);

  mask = m.generate_mask_for_bins_without_contacts();
  REQUIRE(mask.size() == m.ncols());
  for (auto i = 0UL; i < mask.size(); ++i) {
    CHECK((i % 2 != 0) == mask[i]);
  }

  for (auto i = 0UL; i < m.ncols(); ++i) {
    for (auto j = i; j < m.ncols(); ++j) {
      // Each row/column will have half of the fields filled with zeros, and the remaining half with
      // ones: All bits in the bitmask should be set to 1
      m.set(i, j, !(i % 2 == 0 && j % 2 == 0));
    }
  }

  mask = m.generate_mask_for_bins_without_contacts();
  REQUIRE(mask.size() == m.ncols());

  for (auto i = 0UL; i < mask.size(); ++i) {
    CHECK(mask[i]);
  }
  // m.print(true);
}

/*
TEST_CASE("CMatrix 50x5'000", "[cmatrix][long]") {
  const std::string file_path = "../../test/data/symm_matrix_5000_50.tsv.bz2";
  REQUIRE(std::filesystem::exists(file_path));
  test_with_large_matrix(file_path);  // NOLINT
}

TEST_CASE("CMatrix 50x50'000", "[cmatrix][long]") {
  const std::string file_path = "../../test/data/symm_matrix_50000_50.tsv.bz2";
  REQUIRE(std::filesystem::exists(file_path));
  test_with_large_matrix(file_path);  // NOLINT
}
 */
}  // namespace modle::test::cmatrix
