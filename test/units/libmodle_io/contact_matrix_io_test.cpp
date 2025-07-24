// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <hictk/cooler.hpp>
#include <random>

#include "modle/common/random.hpp"
#include "modle/io/contact_matrix_dense.hpp"
#include "modle/test/tmpdir.hpp"

namespace modle::test {
inline const TmpDir testdir{true};  // NOLINT(cert-err58-cpp)
}  // namespace modle::test

namespace modle::io::test {

constexpr auto& testdir = modle::test::testdir;

template <class N>
[[nodiscard]] static ContactMatrixDense<N> init_dense_matrix(random::PRNG_t& rand_eng,
                                                             std::size_t nrows, std::size_t ncols,
                                                             double nnz_fraction = 0.5) {
  ContactMatrixDense<N> m{nrows, ncols};
  for (std::size_t i = 0; i < ncols; ++i) {
    for (std::size_t j = i; j < ncols && j - i < nrows; ++j) {
      if (random::bernoulli_trial{nnz_fraction}(rand_eng)) {
        m.set(i, j, random::uniform_int_distribution<std::uint32_t>{1, 100}(rand_eng));
      }
    }
  }
  REQUIRE(m.get_n_of_missed_updates() == 0);
  return m;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense to Cooler Roundtrip (square)", "[io][matrix][short]") {
  const auto test_file = testdir() / "contact_matrix_dense_to_cooler_roundtrip_square.cool";

  constexpr std::size_t nrows = 50;
  constexpr std::size_t ncols = 50;
  constexpr bp_t bin_size = 10;
  constexpr bp_t chrom_size = bin_size * static_cast<bp_t>(ncols);

  std::random_device rd{};
  auto rand_eng = random::PRNG(rd());

  const auto m1 = init_dense_matrix<std::uint32_t>(rand_eng, nrows, ncols);

  const hictk::Reference chroms{{0, "chr1", chrom_size}};
  {
    auto f = init_cooler_file<std::int32_t>(test_file, false, chroms,
                                            hictk::cooler::Attributes::init(bin_size));
    append_contact_matrix_to_cooler(f, "chr1", m1);
    REQUIRE(m1.get_nnz() == f.attributes().nnz);
  }

  const auto m2 = read_contact_matrix_from_cooler<std::uint32_t>(test_file, "chr1");
  CHECK(m1.get_nnz() == m2.get_nnz());

  for (std::size_t i = 0; i < ncols; ++i) {
    for (std::size_t j = i; j < ncols; ++j) {
      CHECK(m1.get(i, j) == m2.get(i, j));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense to Cooler Roundtrip", "[io][matrix][long]") {
  const auto test_file = testdir() / "contact_matrix_dense_to_cooler_roundtrip.cool";

  constexpr std::size_t nrows = 100;
  constexpr std::size_t ncols = 1000;
  constexpr bp_t bin_size = 10;
  constexpr bp_t chrom_size = bin_size * static_cast<bp_t>(ncols);

  std::random_device rd{};
  auto rand_eng = random::PRNG(rd());

  const auto m1 = init_dense_matrix<std::uint32_t>(rand_eng, nrows, ncols);
  const hictk::Reference chroms{{0, "chr1", chrom_size}};
  {
    auto f = init_cooler_file<std::int32_t>(test_file, false, chroms,
                                            hictk::cooler::Attributes::init(bin_size));
    append_contact_matrix_to_cooler(f, "chr1", m1);
    REQUIRE(m1.get_nnz() == f.attributes().nnz);
  }

  const auto m2 = read_contact_matrix_from_cooler<std::uint32_t>(test_file, "chr1");
  CHECK(m1.get_nnz() == m2.get_nnz());

  for (std::size_t i = 0; i < ncols; ++i) {
    for (std::size_t j = i; j < ncols; ++j) {
      CHECK(m1.get(i, j) == m2.get(i, j));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense to Cooler Roundtrip (FP)", "[io][matrix][long]") {
  const auto test_file = testdir() / "contact_matrix_dense_fp_to_cooler_roundtrip.cool";

  constexpr std::size_t nrows = 100;
  constexpr std::size_t ncols = 1000;
  constexpr bp_t bin_size = 10;
  constexpr bp_t chrom_size = bin_size * static_cast<bp_t>(ncols);

  std::random_device rd{};
  auto rand_eng = random::PRNG(rd());

  const auto m1 = init_dense_matrix<double>(rand_eng, nrows, ncols);
  const hictk::Reference chroms{{0, "chr1", chrom_size}};
  {
    auto f = init_cooler_file<double>(test_file, false, chroms,
                                      hictk::cooler::Attributes::init<double>(bin_size));
    append_contact_matrix_to_cooler(f, "chr1", m1);
    REQUIRE(m1.get_nnz() == f.attributes().nnz);
  }

  const auto m2 = read_contact_matrix_from_cooler<double>(test_file, "chr1");
  CHECK(m1.get_nnz() == m2.get_nnz());

  for (std::size_t i = 0; i < ncols; ++i) {
    for (std::size_t j = i; j < ncols; ++j) {
      CHECK_THAT(m1.get(i, j), Catch::Matchers::WithinRel(m2.get(i, j)));
    }
  }
}

}  // namespace modle::io::test
