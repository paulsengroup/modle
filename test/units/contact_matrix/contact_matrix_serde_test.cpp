// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/contact_matrix_serde.hpp"

#include <absl/time/clock.h>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <limits>
#include <stdexcept>  // for runtime_error
#include <string>     // for string
#include <type_traits>
#include <utility>  // for pair, move
#include <vector>   // for vector, allocator

#include "modle/common/common.hpp"  // for u32
#include "modle/contact_matrix_dense.hpp"
#include "modle/contact_matrix_sparse.hpp"
#include "modle/test/self_deleting_folder.hpp"

namespace modle::test {
inline const SelfDeletingFolder testdir{true};  // NOLINT(cert-err58-cpp)
}  // namespace modle::test

namespace modle::cmatrix::test {

constexpr auto& testdir = modle::test::testdir;

[[maybe_unused]] static const std::filesystem::path& data_dir() {
  static const std::filesystem::path data_dir{"test/data/unit_tests"};
  return data_dir;
}

template <class ContactMatrixT>
static void create_random_matrix(ContactMatrixT& m, usize nnz, u64 seed = 8336046165695760686ULL) {
  using N = typename ContactMatrixT::value_type;
  assert(nnz <= m.npixels());

  auto rand_eng = random::PRNG(seed);

  auto contact_gen = [&rand_eng]() {
    if constexpr (std::is_floating_point_v<N>) {
      return random::uniform_real_distribution<N>{1, 65553}(rand_eng);
    } else {
      u64 max_ = std::min(u64(65553), static_cast<u64>((std::numeric_limits<N>::max)()));
      return random::uniform_int_distribution<N>{1, static_cast<N>(max_)}(rand_eng);
    }
  };

  auto row_gen = [&]() {
    return random::uniform_int_distribution<usize>{0, m.ncols() - 1}(rand_eng);
  };

  auto col_gen = [&]() {
    return random::uniform_int_distribution<usize>{0, m.nrows() - 1}(rand_eng);
  };

  do {
    for (usize i = m.unsafe_get_nnz(); i < nnz; ++i) {
      const auto row = row_gen();
      const auto col = std::min(row + col_gen(), m.ncols() - 1);
      if (row < m.ncols()) {
        m.set(row, col, contact_gen());
        assert(m.get_n_of_missed_updates() == 0);
      }
    }
  } while (m.unsafe_get_nnz() < nnz);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense serialization", "[cmatrix][serialization][medium]") {
  const auto serialized_matrix_path = testdir() / "cmatrix_dense.serialized.zst";

  constexpr bp_t chrom_size = 250'000'000;
  constexpr bp_t diag_width = 3'000'000;
  constexpr bp_t bin_size = 50'000;
  constexpr usize nnz = 100'000;

  constexpr usize ncols = (chrom_size + bin_size - 1) / bin_size;
  constexpr usize nrows = (diag_width + bin_size - 1) / bin_size;

  using T = contacts_t;

  ContactMatrixDense<T> m1(nrows, ncols);
  ContactMatrixDense<T> m2{};

  create_random_matrix(m1, nnz);
  ContactMatrixSerde<T> serde{};

  const auto t0 = absl::Now();
  serde.write(serialized_matrix_path, m1, "test");
  const auto t1 = absl::Now();

  serde.read(serialized_matrix_path, m2);
  const auto t2 = absl::Now();

  spdlog::debug(
      FMT_STRING(
          "serialization took {}; deserialization took: {} ({} pixels, {} nnz,  {:.2f}MB)\n"),
      absl::FormatDuration(t1 - t0), absl::FormatDuration(t2 - t1), m1.npixels(), m1.get_nnz(),
      double(std::filesystem::file_size(serialized_matrix_path)) / 1.0e6);

  REQUIRE(m1.nrows() == m2.nrows());
  REQUIRE(m1.ncols() == m2.ncols());

  CHECK(m1.get_nnz() == m2.get_nnz());
  CHECK(m1.get_tot_contacts() == m2.get_tot_contacts());
  CHECK(m1.get_n_of_missed_updates() == m2.get_n_of_missed_updates());

  const auto& buff1 = m1.get_raw_count_vector();
  const auto& buff2 = m2.get_raw_count_vector();
  REQUIRE(buff1.size() == buff2.size());

  for (usize i = 0; i < buff1.size(); ++i) {
    CHECK(buff1[i] == buff2[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixSparse serialization", "[cmatrix][serialization][long]") {
  const auto serialized_matrix_path = testdir() / "cmatrix_sparse.serialized.zst";

  constexpr bp_t chrom_size = 250'000'000;
  constexpr bp_t diag_width = 30'000'000;
  constexpr bp_t bin_size = 50'000;
  constexpr usize nnz = 100'000;

  constexpr usize ncols = (chrom_size + bin_size - 1) / bin_size;
  constexpr usize nrows = (diag_width + bin_size - 1) / bin_size;

  using T = contacts_t;

  ContactMatrixSparse<T> m1(nrows, ncols);
  ContactMatrixSparse<T> m2{};
  create_random_matrix(m1, nnz);

  ContactMatrixSerde<T> serde{};

  const auto t0 = absl::Now();
  serde.write(serialized_matrix_path, m1, "test");
  const auto t1 = absl::Now();

  serde.read(serialized_matrix_path, m2);
  const auto t2 = absl::Now();

  spdlog::debug(
      FMT_STRING(
          "serialization took {}; deserialization took: {} ({} pixels, {} nnz,  {:.2f}MB)\n"),
      absl::FormatDuration(t1 - t0), absl::FormatDuration(t2 - t1), m1.npixels(), m1.get_nnz(),
      double(std::filesystem::file_size(serialized_matrix_path)) / 1.0e6);

  REQUIRE(m1.nrows() == m2.nrows());
  REQUIRE(m1.ncols() == m2.ncols());

  CHECK(m1.get_nnz() == m2.get_nnz());
  CHECK(m1.get_tot_contacts() == m2.get_tot_contacts());
  CHECK(m1.get_n_of_missed_updates() == m2.get_n_of_missed_updates());

  for (usize i = 0; i < nrows; ++i) {
    for (usize j = 0; j < ncols; ++j) {
      CHECK(m1.get(i, j) == m2.get(i, j));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense serialization (destructive)", "[cmatrix][serialization][medium]") {
  const auto serialized_matrix_path = testdir() / "cmatrix_dense_destructive.serialized.zst";

  constexpr bp_t chrom_size = 250'000'000;
  constexpr bp_t diag_width = 3'000'000;
  constexpr bp_t bin_size = 50'000;
  constexpr usize nnz = 100'000;

  constexpr usize ncols = (chrom_size + bin_size - 1) / bin_size;
  constexpr usize nrows = (diag_width + bin_size - 1) / bin_size;

  using T = contacts_t;

  ContactMatrixDense<T> m1(nrows, ncols);
  create_random_matrix(m1, nnz);
  const auto m2 = m1;

  ContactMatrixSerde<T> serde{};

  const auto t0 = absl::Now();
  serde.write_destructive(serialized_matrix_path, m1, "test");
  const auto t1 = absl::Now();
  CHECK(m1.nrows() == 0);
  CHECK(m1.ncols() == 0);
  CHECK(m1.get_nnz() == 0);

  serde.read(serialized_matrix_path, m1);
  const auto t2 = absl::Now();

  spdlog::debug(
      FMT_STRING(
          "serialization took {}; deserialization took: {} ({} pixels, {} nnz,  {:.2f}MB)\n"),
      absl::FormatDuration(t1 - t0), absl::FormatDuration(t2 - t1), m1.npixels(), m1.get_nnz(),
      double(std::filesystem::file_size(serialized_matrix_path)) / 1.0e6);

  REQUIRE(m1.nrows() == m2.nrows());
  REQUIRE(m1.ncols() == m2.ncols());

  CHECK(m1.get_nnz() == m2.get_nnz());
  CHECK(m1.get_tot_contacts() == m2.get_tot_contacts());
  CHECK(m1.get_n_of_missed_updates() == m2.get_n_of_missed_updates());

  for (usize i = 0; i < nrows; ++i) {
    for (usize j = 0; j < ncols; ++j) {
      CHECK(m1.get(i, j) == m2.get(i, j));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixSparse serialization (destructive)", "[cmatrix][serialization][long]") {
  const auto serialized_matrix_path = testdir() / "cmatrix_sparse_destructive.serialized.zst";

  constexpr bp_t chrom_size = 250'000'000;
  constexpr bp_t diag_width = 30'000'000;
  constexpr bp_t bin_size = 50'000;
  constexpr usize nnz = 100'000;

  constexpr usize ncols = (chrom_size + bin_size - 1) / bin_size;
  constexpr usize nrows = (diag_width + bin_size - 1) / bin_size;

  using T = contacts_t;

  ContactMatrixSparse<T> m1(nrows, ncols);
  create_random_matrix(m1, nnz);
  const auto m2 = m1;

  ContactMatrixSerde<T> serde{};

  const auto t0 = absl::Now();
  serde.write_destructive(serialized_matrix_path, m1, "test");
  const auto t1 = absl::Now();
  CHECK(m1.nrows() == 0);
  CHECK(m1.ncols() == 0);
  CHECK(m1.get_nnz() == 0);

  serde.read(serialized_matrix_path, m1);
  const auto t2 = absl::Now();

  spdlog::debug(
      FMT_STRING(
          "serialization took {}; deserialization took: {} ({} pixels, {} nnz,  {:.2f}MB)\n"),
      absl::FormatDuration(t1 - t0), absl::FormatDuration(t2 - t1), m1.npixels(), m1.get_nnz(),
      double(std::filesystem::file_size(serialized_matrix_path)) / 1.0e6);

  REQUIRE(m1.nrows() == m2.nrows());
  REQUIRE(m1.ncols() == m2.ncols());

  CHECK(m1.get_nnz() == m2.get_nnz());
  CHECK(m1.get_tot_contacts() == m2.get_tot_contacts());
  CHECK(m1.get_n_of_missed_updates() == m2.get_n_of_missed_updates());

  for (usize i = 0; i < nrows; ++i) {
    for (usize j = 0; j < ncols; ++j) {
      CHECK(m1.get(i, j) == m2.get(i, j));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixSparse to ContactMatrixDense", "[cmatrix][serialization][medium]") {
  const auto serialized_matrix_path = testdir() / "cmatrix_sparse_to_dense.serialized.zst";

  constexpr bp_t chrom_size = 250'000'000;
  constexpr bp_t diag_width = 3'000'000;
  constexpr bp_t bin_size = 50'000;
  constexpr usize nnz = 100'000;

  constexpr usize ncols = (chrom_size + bin_size - 1) / bin_size;
  constexpr usize nrows = (diag_width + bin_size - 1) / bin_size;

  using T = contacts_t;

  ContactMatrixSparse<T> m1(nrows, ncols);
  ContactMatrixDense<T> m2{};
  create_random_matrix(m1, nnz);

  ContactMatrixSerde<T> serde{};

  const auto t0 = absl::Now();
  serde.write(serialized_matrix_path, m1, "test");
  const auto t1 = absl::Now();

  serde.read(serialized_matrix_path, m2);
  const auto t2 = absl::Now();

  spdlog::debug(
      FMT_STRING(
          "serialization took {}; deserialization took: {} ({} pixels, {} nnz,  {:.2f}MB)\n"),
      absl::FormatDuration(t1 - t0), absl::FormatDuration(t2 - t1), m1.npixels(), m1.get_nnz(),
      double(std::filesystem::file_size(serialized_matrix_path)) / 1.0e6);

  REQUIRE(m1.nrows() == m2.nrows());
  REQUIRE(m1.ncols() == m2.ncols());

  CHECK(m1.get_nnz() == m2.get_nnz());
  CHECK(m1.get_tot_contacts() == m2.get_tot_contacts());
  CHECK(m1.get_n_of_missed_updates() == m2.get_n_of_missed_updates());

  for (usize i = 0; i < nrows; ++i) {
    for (usize j = 0; j < ncols; ++j) {
      CHECK(m1.get(i, j) == m2.get(i, j));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense to ContactMatrixSparse", "[cmatrix][serialization][medium]") {
  const auto serialized_matrix_path = testdir() / "cmatrix_dense_to_sparse.serialized.zst";

  constexpr bp_t chrom_size = 250'000'000;
  constexpr bp_t diag_width = 3'000'000;
  constexpr bp_t bin_size = 50'000;
  constexpr usize nnz = 100'000;

  constexpr usize ncols = (chrom_size + bin_size - 1) / bin_size;
  constexpr usize nrows = (diag_width + bin_size - 1) / bin_size;

  using T = contacts_t;

  ContactMatrixDense<> m1(nrows, ncols);
  create_random_matrix(m1, nnz);

  ContactMatrixSerde<T> serde{};

  const auto t0 = absl::Now();
  serde.write(serialized_matrix_path, m1, "test");
  const auto t1 = absl::Now();

  ContactMatrixSparse<> m2{};
  serde.read(serialized_matrix_path, m2);
  const auto t2 = absl::Now();

  spdlog::debug(
      FMT_STRING(
          "serialization took {}; deserialization took: {} ({} pixels, {} nnz,  {:.2f}MB)\n"),
      absl::FormatDuration(t1 - t0), absl::FormatDuration(t2 - t1), m1.npixels(), m1.get_nnz(),
      double(std::filesystem::file_size(serialized_matrix_path)) / 1.0e6);

  REQUIRE(m1.nrows() == m2.nrows());
  REQUIRE(m1.ncols() == m2.ncols());

  CHECK(m1.get_nnz() == m2.get_nnz());
  CHECK(m1.get_tot_contacts() == m2.get_tot_contacts());
  CHECK(m1.get_n_of_missed_updates() == m2.get_n_of_missed_updates());

  for (usize i = 0; i < nrows; ++i) {
    for (usize j = 0; j < ncols; ++j) {
      CHECK(m1.get(i, j) == m2.get(i, j));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense serialization (FP)", "[cmatrix][serialization][medium]") {
  const auto serialized_matrix_path = testdir() / "cmatrix_dense_fp.serialized.zst";

  constexpr bp_t chrom_size = 250'000'000;
  constexpr bp_t diag_width = 3'000'000;
  constexpr bp_t bin_size = 50'000;
  constexpr usize nnz = 100'000;

  constexpr usize ncols = (chrom_size + bin_size - 1) / bin_size;
  constexpr usize nrows = (diag_width + bin_size - 1) / bin_size;

  using T = double;

  ContactMatrixDense<T> m1(nrows, ncols);
  ContactMatrixDense<T> m2{};

  create_random_matrix(m1, nnz);
  ContactMatrixSerde<T> serde{};

  const auto t0 = absl::Now();
  serde.write(serialized_matrix_path, m1, "test");
  const auto t1 = absl::Now();

  serde.read(serialized_matrix_path, m2);
  const auto t2 = absl::Now();

  spdlog::debug(
      FMT_STRING(
          "serialization took {}; deserialization took: {} ({} pixels, {} nnz,  {:.2f}MB)\n"),
      absl::FormatDuration(t1 - t0), absl::FormatDuration(t2 - t1), m1.npixels(), m1.get_nnz(),
      double(std::filesystem::file_size(serialized_matrix_path)) / 1.0e6);

  REQUIRE(m1.nrows() == m2.nrows());
  REQUIRE(m1.ncols() == m2.ncols());

  CHECK(m1.get_nnz() == m2.get_nnz());
  CHECK(Catch::Approx(m1.get_tot_contacts()) == m2.get_tot_contacts());
  CHECK(m1.get_n_of_missed_updates() == m2.get_n_of_missed_updates());

  const auto& buff1 = m1.get_raw_count_vector();
  const auto& buff2 = m2.get_raw_count_vector();
  REQUIRE(buff1.size() == buff2.size());

  for (usize i = 0; i < buff1.size(); ++i) {
    CHECK(buff1[i] == buff2[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixSparse serialization (FP)", "[cmatrix][serialization][long]") {
  const auto serialized_matrix_path = testdir() / "cmatrix_sparse_fp.serialized.zst";

  constexpr bp_t chrom_size = 250'000'000;
  constexpr bp_t diag_width = 30'000'000;
  constexpr bp_t bin_size = 50'000;
  constexpr usize nnz = 100'000;

  constexpr usize ncols = (chrom_size + bin_size - 1) / bin_size;
  constexpr usize nrows = (diag_width + bin_size - 1) / bin_size;

  using T = double;

  ContactMatrixSparse<T> m1(nrows, ncols);
  ContactMatrixSparse<T> m2{};
  create_random_matrix(m1, nnz);

  ContactMatrixSerde<T> serde{};

  const auto t0 = absl::Now();
  serde.write(serialized_matrix_path, m1, "test");
  const auto t1 = absl::Now();

  serde.read(serialized_matrix_path, m2);
  const auto t2 = absl::Now();

  spdlog::debug(
      FMT_STRING(
          "serialization took {}; deserialization took: {} ({} pixels, {} nnz,  {:.2f}MB)\n"),
      absl::FormatDuration(t1 - t0), absl::FormatDuration(t2 - t1), m1.npixels(), m1.get_nnz(),
      double(std::filesystem::file_size(serialized_matrix_path)) / 1.0e6);

  REQUIRE(m1.nrows() == m2.nrows());
  REQUIRE(m1.ncols() == m2.ncols());

  CHECK(m1.get_nnz() == m2.get_nnz());
  CHECK(Catch::Approx(m1.get_tot_contacts()) == m2.get_tot_contacts());
  CHECK(m1.get_n_of_missed_updates() == m2.get_n_of_missed_updates());

  for (usize i = 0; i < nrows; ++i) {
    for (usize j = 0; j < ncols; ++j) {
      CHECK(m1.get(i, j) == m2.get(i, j));
    }
  }
}

}  // namespace modle::cmatrix::test
