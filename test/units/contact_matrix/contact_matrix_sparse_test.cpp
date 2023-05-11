// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/contact_matrix_sparse.hpp"  // for ContactMatrixSparse

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "modle/common/common.hpp"  // for u32

namespace modle::test::cmatrix {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrixSparse simple", "[cmatrix][short]") {
  ContactMatrixSparse<> c(10, 100, ContactMatrixSparse<>::ChunkSize{75});
  CHECK(c.num_blocks() == 15);

  CHECK(c.get(0, 0) == 0);
  c.increment(0, 0);
  CHECK(c.get(0, 0) == 1);
  c.increment(0, 0);
  CHECK(c.get(0, 0) == 2);
  c.subtract(0, 0, 2);
  CHECK(c.get(0, 0) == 0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrixSparse in/decrement", "[cmatrix][short]") {
  ContactMatrixSparse<> m(10, 20, ContactMatrixSparse<>::ChunkSize{50});
  m.increment(0, 0);
  m.increment(15, 15);

  CHECK(m.get_tot_contacts() == 2);
  CHECK(m.get(0, 0) == 1);

  m.decrement(0, 0);
  CHECK(m.get_tot_contacts() == 1);
  CHECK(m.get(0, 0) == 0);

  REQUIRE(m.get_n_of_missed_updates() == 0);
  m.increment(11, 0);
  CHECK(m.get(0, 0) == 0);
  CHECK(m.get_n_of_missed_updates() == 1);
  CHECK(m.get_tot_contacts() == 1);

  if constexpr (utils::ndebug_not_defined()) {
    CHECK_THROWS_WITH(m.increment(25, 25),
                      Catch::Matchers::ContainsSubstring(
                          "detected an out-of-bound read: attempt to access item at"));
    CHECK(m.get_n_of_missed_updates() == 1);
    CHECK(m.get_tot_contacts() == 1);

    CHECK_THROWS_WITH(m.decrement(25, 25),
                      Catch::Matchers::ContainsSubstring(
                          "detected an out-of-bound read: attempt to access item at"));
    CHECK(m.get_n_of_missed_updates() == 1);
    CHECK(m.get_tot_contacts() == 1);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrixSparse test get_nnz", "[cmatrix][short]") {
  ContactMatrixSparse<> m(10, 10);
  m.set(0, 0, 100);
  CHECK(m.get_nnz() == 1);

  m.set(0, 0, 0);
  CHECK(m.get_nnz() == 0);

  m.set(1, 1, 100);
  CHECK(m.get_nnz() == 1);

  m.set(1, 2, 100);
  CHECK(m.get_nnz() == 2);

  m.set(1, 2, 100);
  CHECK(m.get_nnz() == 2);

  m.set(1, 2, 0);
  CHECK(m.get_nnz() == 1);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrixSparse test get_tot_contacts", "[cmatrix][short]") {
  ContactMatrixSparse<> m(10, 10);

  m.set(0, 0, 100);
  CHECK(m.get_tot_contacts() == 100);

  m.set(0, 0, 0);
  CHECK(m.get_tot_contacts() == 0);

  m.set(1, 1, 100);
  CHECK(m.get_tot_contacts() == 100);

  m.set(1, 2, 100);
  CHECK(m.get_tot_contacts() == 200);

  m.set(1, 2, 100);
  CHECK(m.get_tot_contacts() == 200);

  m.set(1, 2, 0);
  CHECK(m.get_tot_contacts() == 100);
}

}  // namespace modle::test::cmatrix
