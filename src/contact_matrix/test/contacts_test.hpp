#pragma once

#include <absl/strings/str_split.h>  // for SplitIterator, Splitter, StrSplit
#include <fmt/format.h>              // for string_view, system_error

#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset, dynamic_bitset<>::ref...
#include <boost/filesystem/path.hpp>                // for exists
#include <catch2/catch.hpp>                         // for AssertionHandler, operator""_catch_sr
#include <cerrno>                                   // for errno
#include <cstdint>                                  // for uint32_t
#include <ext/alloc_traits.h>                       // for __alloc_traits<>::value_type
#include <memory>                                   // for allocator_traits<>::value_type
#include <stdexcept>                                // for runtime_error
#include <string>                                   // for string, getline
#include <utility>                                  // for move
#include <vector>                                   // for vector, allocator

#include "modle/common/utils.hpp"   // for ndebug_defined, parse_numeric_or_throw
#include "modle/compressed_io.hpp"  // for Reader
#include "modle/contacts.hpp"       // for ContactMatrix

namespace modle::test::cmatrix {

TEST_CASE("CMatrix simple", "[cmatrix][short]") {
  ContactMatrix<> c(10, 100);  // NOLINT
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
  compressed_io::Reader r(path_to_file);
  std::string line;
  std::string buff;
  uint32_t n;  // NOLINT
  while (r.getline(line)) {
    std::vector<uint32_t> v;

    for (const auto& tok : absl::StrSplit(line, sep)) {
      modle::utils::parse_numeric_or_throw(tok, n);
      v.push_back(n);
    }
    m.emplace_back(std::move(v));
  }
  REQUIRE(r.eof());
  return m;
}

TEST_CASE("CMatrix 10x200", "[cmatrix][medium]") {
  const std::string file_path = "test/data/unit_tests/symm_matrix_200_10.tsv.gz";
  REQUIRE(boost::filesystem::exists(file_path));
  const auto m1 = load_matrix_from_file(file_path);
  ContactMatrix<> m2(10, 200);  // NOLINT
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
  ContactMatrix<> m(10, 20);  // NOLINT
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

TEST_CASE("CMatrix in/decrement", "[cmatrix][short]") {
  ContactMatrix<> m(10, 20);  // NOLINT
  REQUIRE_NOTHROW(m.increment(0, 0));
  REQUIRE_NOTHROW(m.increment(15, 15));  // NOLINT

  CHECK(m.get_tot_contacts() == 2);  // NOLINT
  CHECK(m.get(0, 0) == 1);

  REQUIRE_NOTHROW(m.decrement(0, 0));
  CHECK(m.get_tot_contacts() == 1);
  CHECK(m.get(0, 0) == 0);

  REQUIRE(m.get_n_of_missed_updates() == 0);
  REQUIRE_NOTHROW(m.increment(11, 0));  // NOLINT
  CHECK(m.get(0, 0) == 0);
  CHECK(m.get_n_of_missed_updates() == 1);
  CHECK(m.get_tot_contacts() == 1);

#ifndef NDEBUG
  CHECK_THROWS_WITH(
      m.increment(25, 25),  // NOLINT
      Catch::Contains("caught an attempt to access element past the end of the contact matrix"));
  CHECK(m.get_n_of_missed_updates() == 1);
  CHECK(m.get_tot_contacts() == 1);

  CHECK_THROWS_WITH(
      m.decrement(25, 25),  // NOLINT
      Catch::Contains("caught an attempt to access element past the end of the contact matrix"));
  CHECK(m.get_n_of_missed_updates() == 1);
  CHECK(m.get_tot_contacts() == 1);
#endif
}

TEST_CASE("CMatrix in/decrement vector", "[cmatrix][short]") {
  ContactMatrix<> m(10, 20);  // NOLINT
  // clang-format off
  const std::vector<std::pair<size_t, size_t>> pixels{{ 0,  0},   // NOLINT
                                                          { 5, 10},         // NOLINT
                                                          { 2,  3},         // NOLINT
                                                          {15,  0},         // NOLINT
                                                          { 7,  1},         // NOLINT
                                                          {25, 25}};        // NOLINT
  // clang-format on
  auto pixels_ = pixels;

  m.increment(absl::MakeSpan(pixels_.data(), pixels_.size() - 1));

  CHECK(m.get_n_of_missed_updates() == 1);
  CHECK(m.get_tot_contacts() == pixels.size() - 2);  // NOLINT
  for (auto i = 0UL; i < pixels.size() - 1; ++i) {
    const auto& [row, col] = pixels[i];
    if (row == 15UL && col == 0UL) {  // NOLINT
      CHECK(m.get(row, col) == 0);
      continue;
    }
    CHECK(m.get(row, col) == 1);
  }

#ifndef NDEBUG
  pixels_ = pixels;
  // Setting thresh to 0 forces this function to switch to a branch that is supposed to be more
  CHECK_THROWS_WITH(  // efficient for large  buffers
      m.increment(absl::MakeSpan(pixels_), 0),
      Catch::Contains("caught an attempt to access element past the end of the contact matrix"));
  CHECK(m.get_n_of_missed_updates() == 2);                 // NOLINT
  CHECK(m.get_tot_contacts() == 2 * (pixels.size() - 2));  // NOLINT
  for (const auto& [row, col] : pixels) {
    if (row == 15UL && col == 0UL) {  // NOLINT
      CHECK(m.get(row, col) == 0);
      continue;
    }
    if (row == 25UL && col == 25UL) {  // NOLINT
      continue;
    }
    CHECK(m.get(row, col) == 2);  // NOLINT
  }

#endif
}

}  // namespace modle::test::cmatrix
