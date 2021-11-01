// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/strings/str_split.h>  // for SplitIterator, Splitter, StrSplit
#include <fmt/format.h>              // for format

#include <algorithm>                                // for max
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset, dynamic_bitset<>::ref...
#include <boost/filesystem/operations.hpp>          // for exists
#include <boost/filesystem/path.hpp>                // for path
#include <boost/process.hpp>
#include <catch2/catch.hpp>  // for operator""_catch_sr, AssertionHandler
#include <stdexcept>         // for runtime_error
#include <string>            // for string
#include <utility>           // for pair, move
#include <vector>            // for vector, allocator

#include "modle/common/common.hpp"  // for u32
#include "modle/common/utils.hpp"   // for parse_numeric_or_throw
#include "modle/compressed_io.hpp"  // for Reader
#include "modle/contacts.hpp"       // for ContactMatrix

namespace modle::test::cmatrix {

inline const boost::filesystem::path data_dir{"test/data/unit_tests"};  // NOLINT

// clang-format off
constexpr auto SCIPY_GAUSSIAN_BLUR_CMD{
    "#!/usr/bin/env python3\n"
    "from scipy.ndimage import gaussian_filter\n"
    "import numpy as np\n"
    "from sys import stdin\n"
    "if __name__ == \"__main__\":\n"
    "    shape = [{:d}, {:d}]\n"
    "    sigma = {:.16e}\n"
    "    buff = stdin.read().replace(\"\\n\", \",\")[:-1]\n"
    "    m = np.fromstring(buff, sep=\",\", dtype=float)\n"
    "    m = gaussian_filter(m.reshape(shape), sigma, truncate={})\n"
    "    print(\",\".join([str(n) for n in m.flatten()]))\n"};
// clang-format on

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix simple", "[cmatrix][short]") {
  ContactMatrix<> c(10, 100);  // NOLINT
  CHECK(c.unsafe_get(0, 0) == 0);
  c.increment(0, 0);
  CHECK(c.unsafe_get(0, 0) == 1);
  c.increment(0, 0);
  CHECK(c.unsafe_get(0, 0) == 2);
  c.subtract(0, 0, 2);
  CHECK(c.unsafe_get(0, 0) == 0);
}

[[nodiscard]] inline std::vector<std::vector<u32>> load_matrix_from_file(
    const std::string& path_to_file, const std::string& sep = "\t") {
  std::vector<std::vector<u32>> m;
  compressed_io::Reader r(path_to_file);
  std::string line;
  std::string buff;
  u32 n;  // NOLINT
  while (r.getline(line)) {
    std::vector<u32> v;

    for (const auto& tok : absl::StrSplit(line, sep)) {
      modle::utils::parse_numeric_or_throw(tok, n);
      v.push_back(n);
    }
    m.emplace_back(std::move(v));
  }
  REQUIRE(r.eof());
  return m;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix 10x200", "[cmatrix][medium]") {
  const auto input_file = data_dir / "symm_matrix_200_10.tsv.gz";
  REQUIRE(boost::filesystem::exists(input_file));
  const auto m1 = load_matrix_from_file(input_file.string());
  ContactMatrix<> m2(10, 200);  // NOLINT
  for (usize i = 0; i < m1.size(); ++i) {
    for (usize j = 0; j < m1[i].size(); ++j) {
      if (m1[i][j] != 0 && j >= i) {
        m2.set(i, j, m1[i][j]);
      }
    }
  }

  const auto m3 = m2.unsafe_generate_symmetric_matrix();
  for (usize i = 0; i < m1.size(); ++i) {
    for (usize j = 0; j < m1[0].size(); ++j) {
      CHECK(m1[i][j] == m3[i][j]);
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix Mask", "[cmatrix][short]") {
  ContactMatrix<> m(10, 20);  // NOLINT
  auto mask = m.unsafe_generate_mask_for_bins_without_contacts();
  REQUIRE(mask.size() == m.ncols());
  CHECK(mask.none());  // Matrix is full of zeros: bitmask should also be all zeros

  for (usize i = 0; i < m.ncols(); ++i) {
    for (auto j = i; j < m.ncols(); ++j) {
      // Even row/columns are all zeros: Bitmask should have zeros at pos corresponding to even idx
      m.set(i, j, !(i % 2 == 0 || j % 2 == 0));  // NOLINT(readability-implicit-bool-conversion)
    }
  }
  // m.print(true);

  mask = m.unsafe_generate_mask_for_bins_without_contacts();
  REQUIRE(mask.size() == m.ncols());
  for (usize i = 0; i < mask.size(); ++i) {
    CHECK((i % 2 != 0) == mask[i]);
  }

  for (usize i = 0; i < m.ncols(); ++i) {
    for (auto j = i; j < m.ncols(); ++j) {
      // Each row/column will have half of the fields filled with zeros, and the remaining half with
      // ones: All bits in the bitmask should be set to 1
      m.set(i, j, !(i % 2 == 0 && j % 2 == 0));  // NOLINT(readability-implicit-bool-conversion)
    }
  }

  mask = m.unsafe_generate_mask_for_bins_without_contacts();
  REQUIRE(mask.size() == m.ncols());

  for (usize i = 0; i < mask.size(); ++i) {
    CHECK(mask[i]);
  }
  // m.print(true);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix in/decrement", "[cmatrix][short]") {
  ContactMatrix<> m(10, 20);  // NOLINT
  m.increment(0, 0);
  m.increment(15, 15);  // NOLINT

  CHECK(m.get_tot_contacts() == 2);  // NOLINT
  CHECK(m.unsafe_get(0, 0) == 1);

  m.decrement(0, 0);
  CHECK(m.get_tot_contacts() == 1);
  CHECK(m.unsafe_get(0, 0) == 0);

  REQUIRE(m.get_n_of_missed_updates() == 0);
  m.increment(11, 0);  // NOLINT
  CHECK(m.unsafe_get(0, 0) == 0);
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix unsafe_get w/ block", "[cmatrix][short]") {
  ContactMatrix<u32> m1(100, 100);  // NOLINT
  // Fill the upper left corner
  for (u32 i = 0; i < 3; ++i) {
    for (u32 j = i; j < 3; ++j) {
      m1.set(i, j, i + j);
    }
  }

  // Fill a region away from the diagonal
  for (u32 i = 20; i < 25; ++i) {    // NOLINT
    for (u32 j = 25; j < 30; ++j) {  // NOLINT
      m1.set(i, j, 1);
    }
  }

  // Fill the lower right corner
  for (u32 i = 97; i < 100; ++i) {        // NOLINT
    for (u32 j = 97; j < 100; ++j) {      // NOLINT
      m1.set(i, j, (i - 97) + (j - 97));  // NOLINT
    }
  }

  // m1.print(true);

  CHECK(m1.unsafe_get(0, 0, 5) == 30);
  CHECK(m1.unsafe_get(22, 27, 5) == 25);
  CHECK(m1.unsafe_get(99, 99, 5) == 70);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix unsafe_get w/ block small", "[cmatrix][short]") {
  const auto reference_file = data_dir / "contacts_chr1_bs9_small.tsv";
  const auto input_file = data_dir / "contacts_chr1_raw_small.tsv";

  const usize block_size = 9;

  const auto reference_matrix = [&]() {
    ContactMatrix<> m;
    m.unsafe_import_from_txt(reference_file);
    return m;
  }();

  const auto input_matrix = [&]() {
    ContactMatrix<> m;
    m.unsafe_import_from_txt(input_file);
    return m;
  }();

  REQUIRE(input_matrix.nrows() == reference_matrix.nrows());
  REQUIRE(input_matrix.ncols() == reference_matrix.ncols());

  for (usize i = 0; i < input_matrix.nrows(); ++i) {
    for (usize j = 0; j < input_matrix.ncols(); ++j) {
      CHECK(input_matrix.unsafe_get(i, j, block_size) == reference_matrix.unsafe_get(i, j));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix get column", "[cmatrix][short]") {
  ContactMatrix<> c(10, 100);  // NOLINT

  const usize col = 25;  // Set a column of pixels to 1
  for (usize i = 0; i < c.nrows(); ++i) {
    c.set(col - i, col, static_cast<u32>(i));
  }

  REQUIRE(c.get_tot_contacts() == 45);  // NOLINT
  REQUIRE(c.get_n_of_missed_updates() == 0);

  std::vector<u32> buff;
  // Test getting a column of pixels
  c.unsafe_get_column(col, buff);
  REQUIRE(buff.size() == c.nrows());
  for (usize i = 0; i < c.nrows(); ++i) {
    CHECK(buff[i] == i);
  }

  // Test getting a column of pixel skipping the first 5 values (starting from the diagonal)
  const usize offset = 5;
  c.unsafe_get_column(col, buff, offset);

  REQUIRE(buff.size() == c.nrows() - offset);
  for (usize i = 0; i < c.nrows() - offset; ++i) {
    CHECK(buff[i] == i + offset);
  }

  // Test getting the first column of pixels
  c.set(0, 0, 1);
  c.unsafe_get_column(0, buff);
  REQUIRE(buff.size() == c.nrows());
  CHECK(buff.front() == 1);
  for (usize i = 1; i < c.nrows(); ++i) {
    CHECK(buff[i] == 0);
  }

  // Test getting the last column of pixels (which is truncated)
  c.set(c.ncols(), c.ncols(), 1);
  c.unsafe_get_column(c.ncols(), buff);
  REQUIRE(buff.size() == 1);
  CHECK(buff.front() == 1);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix get row", "[cmatrix][short]") {
  ContactMatrix<> c(10, 100);  // NOLINT

  const usize row = 25;  // Set a row of pixels to 1
  for (usize i = 0; i < c.nrows(); ++i) {
    c.set(row, row + i, static_cast<u32>(i));
  }

  REQUIRE(c.get_tot_contacts() == 45);  // NOLINT
  REQUIRE(c.get_n_of_missed_updates() == 0);

  std::vector<u32> buff;
  c.unsafe_get_row(row, buff);

  // Test getting a row of pixels
  REQUIRE(buff.size() == c.nrows());
  for (usize i = 0; i < c.nrows(); ++i) {
    CHECK(buff[i] == i);
  }

  // Test getting a row of pixel skipping the first 5 values (starting from the diagonal)
  const usize offset = 5;
  c.unsafe_get_row(row, buff, offset);

  REQUIRE(buff.size() == c.nrows() - offset);
  for (usize i = 0; i < c.nrows() - offset; ++i) {
    CHECK(buff[i] == offset + i);
  }

  // Test getting the first row of pixels
  for (usize i = 0; i < c.nrows(); ++i) {
    c.set(0, i, static_cast<u32>(i));
  }

  c.unsafe_get_row(0, buff);
  REQUIRE(buff.size() == 10);
  for (usize i = 0; i < c.nrows(); ++i) {
    CHECK(buff[i] == i);
  }

  // Test getting the last row of pixels (which is truncated)
  c.set(c.ncols() - 1, c.ncols() - 1, 1);
  c.unsafe_get_row(c.ncols() - 1, buff);
  REQUIRE(buff.size() == 1);
  CHECK(buff.front() == 1);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix blur", "[cmatrix][short]") {
  const auto reference_file = data_dir / "contacts_chr1_bs9_small.tsv";
  const auto input_matrix = [&]() {
    ContactMatrix<> m;
    m.unsafe_import_from_txt(reference_file);
    return m;
  }();

  auto compute_reference_matrix = [&input_matrix](auto shape, auto sigma, auto trunc) {
    boost::process::ipstream stdout_stream;
    boost::process::opstream stdin_stream;
    auto py = boost::process::child(
        boost::process::search_path("python3").string(), "-c",
        fmt::format(FMT_STRING(SCIPY_GAUSSIAN_BLUR_CMD), shape, shape, sigma, trunc),
        boost::process::std_in<stdin_stream, boost::process::std_out> stdout_stream);
    assert(py.running());  // NOLINT
    {
      std::vector<contacts_t> buff(input_matrix.ncols());
      for (usize i = 0; i < input_matrix.ncols(); ++i) {
        buff.clear();
        for (usize j = 0; j < input_matrix.ncols(); ++j) {
          buff.push_back(input_matrix.unsafe_get(i, j));
        }
        const auto sbuff = fmt::format(FMT_STRING("{}\n"), fmt::join(buff, ","));
        stdin_stream.write(sbuff.data(), static_cast<std::streamsize>(sbuff.size()));
      }
      stdin_stream.flush();
    }
    stdin_stream.pipe().close();

    const auto reference_matrix = [&]() {
      std::string sbuff;
      std::getline(stdout_stream, sbuff);

      std::vector<double> buff;
      for (const auto& tok : absl::StrSplit(sbuff, ',')) {
        buff.push_back(utils::parse_numeric_or_throw<double>(tok));
      }
      assert(buff.size() == input_matrix.nrows() * input_matrix.ncols());  // NOLINT
      ContactMatrix<double> m(input_matrix.nrows(), input_matrix.ncols());
      for (usize i = 0; i < input_matrix.nrows(); ++i) {
        for (auto j = i; j < input_matrix.ncols(); ++j) {
          m.set(i, j, buff[(i * m.nrows()) + j]);
        }
      }
      return m;
    }();
    py.wait();
    return reference_matrix;
  };

  constexpr std::array<double, 3> sigmas{0.5, 1.0, 1.5};
  constexpr std::array<double, 3> cutoffs{3.0, 3.0, 3.0};

  for (usize i = 0; i < sigmas.size(); ++i) {
    const auto reference_matrix =
        compute_reference_matrix(input_matrix.ncols(), sigmas[i], cutoffs[i]);
    const auto blurred_matrix = input_matrix.blur(sigmas[i]);

    for (usize j = 4; j < input_matrix.nrows(); ++j) {       // NOLINT
      for (auto k = j; k < input_matrix.ncols() - 4; ++k) {  // NOLINT
        CHECK(Approx(reference_matrix.unsafe_get(j, k)) == blurred_matrix.unsafe_get(j, k));
      }
    }
  }
}

}  // namespace modle::test::cmatrix
