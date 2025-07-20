// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/contact_matrix_dense.hpp"  // for ContactMatrixDense

#include <fmt/format.h>  // for format

#include <BS_thread_pool.hpp>  // for BS::light_thread_pool
#include <algorithm>           // for generate, max
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>  // for path
#include <stdexcept>   // for runtime_error
#include <string>      // for string
#include <thread>      // for sleep_for
#include <utility>     // for pair, move
#include <vector>      // for vector, allocator

#include "modle/common/common.hpp"  // for u32
#include "modle/compressed_io/compressed_io.hpp"

namespace modle::cmatrix::test {

[[maybe_unused]] static const std::filesystem::path& data_dir() {
  static const std::filesystem::path data_dir{"test/data/unit_tests"};
  return data_dir;
}

// See https://github.com/catchorg/Catch2/blob/v3.2.1/src/catch2/catch_approx.cpp#L27-L32
constexpr double DEFAULT_FP_TOLERANCE =
    static_cast<double>(std::numeric_limits<float>::epsilon() * 100);

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense simple", "[cmatrix][short]") {
  ContactMatrixDense<> c(10, 100);
  CHECK(c.get(0, 0) == 0);
  c.increment(0, 0);
  CHECK(c.get(0, 0) == 1);
  c.increment(0, 0);
  CHECK(c.get(0, 0) == 2);
  c.subtract(0, 0, 2);
  CHECK(c.get(0, 0) == 0);
}

[[nodiscard]] static std::vector<std::vector<u32>> load_matrix_from_file(
    const std::filesystem::path& path_to_matrix, const std::string& sep = "\t") {
  std::vector<std::vector<u32>> m;
  compressed_io::Reader r(path_to_matrix);
  std::string line;
  std::string buff;
  u32 n{};
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
TEST_CASE("ContactMatrixDense 10x200", "[cmatrix][short]") {
  const auto input_file = data_dir() / "contact_matrices" / "symmetric_matrix_200_10.tsv.xz";
  REQUIRE(std::filesystem::exists(input_file));
  const auto m1 = load_matrix_from_file(input_file);
  ContactMatrixDense<> m2(10, 200);
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
TEST_CASE("ContactMatrixDense in/decrement", "[cmatrix][short]") {
  ContactMatrixDense<> m(10, 20);
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
TEST_CASE("ContactMatrixDense get w/ block", "[cmatrix][short]") {
  ContactMatrixDense<u32> m1(100, 100);
  // Fill the upper left corner
  for (u32 i = 0; i < 3; ++i) {
    for (u32 j = i; j < 3; ++j) {
      m1.set(i, j, i + j);
    }
  }

  // Fill a region away from the diagonal
  for (u32 i = 20; i < 25; ++i) {
    for (u32 j = 25; j < 30; ++j) {
      m1.set(i, j, 1);
    }
  }

  // Fill the lower right corner
  for (u32 i = 97; i < 100; ++i) {
    for (u32 j = 97; j < 100; ++j) {
      m1.set(i, j, (i - 97) + (j - 97));
    }
  }

  // m1.print(true);

  CHECK(m1.unsafe_get_block(0, 0, 5) == 30);
  CHECK(m1.unsafe_get_block(22, 27, 5) == 25);
  CHECK(m1.unsafe_get_block(99, 99, 5) == 70);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense get w/ block small", "[cmatrix][short]") {
  const auto parent = data_dir() / "contact_matrices";
  const auto reference_file = parent / "contact_matrix_chr1_bs9_small.tsv.xz";
  const auto input_file = parent / "contact_matrix_chr1_raw_small.tsv.xz";

  const usize block_size = 9;

  const auto reference_matrix = ContactMatrixDense<>::from_txt(reference_file);
  const auto input_matrix = ContactMatrixDense<>::from_txt(input_file);

  REQUIRE(input_matrix.nrows() == reference_matrix.nrows());
  REQUIRE(input_matrix.ncols() == reference_matrix.ncols());

  for (usize i = 0; i < input_matrix.nrows(); ++i) {
    for (usize j = 0; j < input_matrix.ncols(); ++j) {
      CHECK(input_matrix.unsafe_get_block(i, j, block_size) == reference_matrix.get(i, j));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense get column", "[cmatrix][short]") {
  ContactMatrixDense<> c(10, 100);

  const usize col = 25;  // Set a column of pixels to 1
  for (usize i = 0; i < c.nrows(); ++i) {
    c.set(col - i, col, static_cast<u32>(i));
  }

  REQUIRE(c.get_tot_contacts() == 45);
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
  c.set(c.ncols() - 1, c.ncols() - 1, 1);
  c.unsafe_get_column(c.ncols() - 1, buff);
  REQUIRE(buff.size() == 1);
  CHECK(buff.front() == 1);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense get row", "[cmatrix][short]") {
  ContactMatrixDense<> c(10, 100);

  const usize row = 25;  // Set a row of pixels to 1
  for (usize i = 0; i < c.nrows(); ++i) {
    c.set(row, row + i, static_cast<u32>(i));
  }

  REQUIRE(c.get_tot_contacts() == 45);
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

static void contact_matrix_dense_blur_helper(double sigma, double truncate,
                                             BS::light_thread_pool* tpool = nullptr,
                                             double tolerance = DEFAULT_FP_TOLERANCE) {
  const auto path_to_input_matrix =
      data_dir() / "contact_matrices" / "contact_matrix_dense_int_001.tsv.xz";

  const auto path_to_reference_matrix =
      data_dir() / "contact_matrices" / "blurred" /
      fmt::format(FMT_STRING("contact_matrix_dense_int_001_blurred_{:.2f}.tsv.xz"), sigma);

  const auto input_matrix = ContactMatrixDense<double>::from_txt(path_to_input_matrix);
  const auto reference_matrix = ContactMatrixDense<double>::from_txt(path_to_reference_matrix);
  const auto blurred_matrix = input_matrix.blur(sigma, truncate, tpool);

  REQUIRE(reference_matrix.nrows() == blurred_matrix.nrows());
  REQUIRE(reference_matrix.ncols() == blurred_matrix.ncols());

  for (usize i = 0; i < input_matrix.nrows(); ++i) {
    for (auto j = i; j < input_matrix.ncols(); ++j) {
      CHECK_THAT(reference_matrix.get(i, j),
                 Catch::Matchers::WithinRel(blurred_matrix.get(i, j), tolerance));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense blur (sigma=0.01)", "[cmatrix][long]") {
  contact_matrix_dense_blur_helper(0.01, 3.5);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense blur (sigma=0.5)", "[cmatrix][long]") {
  contact_matrix_dense_blur_helper(0.5, 3.5);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense blur (sigma=1.0)", "[cmatrix][long]") {
  contact_matrix_dense_blur_helper(1.0, 3.5);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense blur (sigma=1.5)", "[cmatrix][long]") {
  contact_matrix_dense_blur_helper(1.5, 3.5);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense blur parallel (sigma=5.0)", "[cmatrix][long][exclusive]") {
  BS::light_thread_pool tpool{2};
  contact_matrix_dense_blur_helper(5.0, 3.5, &tpool);
}

static void contact_matrix_dense_dog_helper(double sigma1, double sigma2, double truncate,
                                            BS::light_thread_pool* tpool = nullptr) {
  assert(sigma1 < sigma2);
  const auto path_to_input_matrix =
      data_dir() / "contact_matrices" / "contact_matrix_dense_int_001.tsv.xz";

  const auto path_to_reference_matrix =
      data_dir() / "contact_matrices" / "diff_of_gaussians" /
      fmt::format(FMT_STRING("contact_matrix_dense_int_001_dog_{:.2f}_{:.2f}.tsv.xz"), sigma1,
                  sigma2);

  const auto input_matrix = ContactMatrixDense<double>::from_txt(path_to_input_matrix);

  const auto reference_matrix = ContactMatrixDense<double>::from_txt(path_to_reference_matrix);
  const auto dog_matrix = input_matrix.diff_of_gaussians(sigma1, sigma2, truncate,
                                                         std::numeric_limits<double>::lowest(),
                                                         std::numeric_limits<double>::max(), tpool);

  REQUIRE(reference_matrix.nrows() == dog_matrix.nrows());
  REQUIRE(reference_matrix.ncols() == dog_matrix.ncols());

  for (usize i = 0; i < input_matrix.nrows(); ++i) {
    for (auto j = i; j < input_matrix.ncols(); ++j) {
      CHECK_THAT(reference_matrix.get(i, j),
                 Catch::Matchers::WithinRel(dog_matrix.get(i, j), 1.0e-6));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense difference of Gaussians", "[cmatrix][long]") {
  contact_matrix_dense_dog_helper(1.0, 1.6, 3.0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense difference of Gaussians - parallel", "[cmatrix][long][exclusive]") {
  BS::light_thread_pool tpool{2};
  contact_matrix_dense_dog_helper(1.0, 1.6, 3.0, &tpool);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense test get_nnz", "[cmatrix][short]") {
  ContactMatrixDense<> m(10, 10);
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
TEST_CASE("ContactMatrixDense test get_tot_contacts", "[cmatrix][short]") {
  ContactMatrixDense<> m(10, 10);
  m.set(1, 1, 100);
  CHECK(m.get_tot_contacts() == 100);

  m.set(1, 2, 100);
  CHECK(m.get_tot_contacts() == 200);

  m.set(1, 2, 100);
  CHECK(m.get_tot_contacts() == 200);

  m.set(1, 2, 0);
  CHECK(m.get_tot_contacts() == 100);
}

#ifdef MODLE_ENABLE_SANITIZER_THREAD
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense pixel locking TSAN", "[cmatrix][long]") {
  ContactMatrixDense<i64> m(100, 10);
  std::atomic<bool> stop_sig = false;

  // Increment a random pixel by 1 until stop_sign is flipped
  auto generate_contacts = [&](u64 seed) {
    auto rand_eng = random::PRNG(seed);
    random::uniform_int_distribution<usize> idx_gen(0, m.ncols() - 1);

    usize i = 0;
    while (!stop_sig) {
      m.increment(idx_gen(rand_eng), idx_gen(rand_eng));
      ++i;
    }
    return i;
  };

  const auto nthreads = std::max(u32(2), std::thread::hardware_concurrency());
  BS::light_thread_pool tpool(nthreads);

  // Submit nthreads - 1 tasks generating contacts
  std::vector<std::future<usize>> num_contacts_generated(nthreads - 1);
  std::generate(num_contacts_generated.begin(), num_contacts_generated.end(),
                [&]() { return tpool.submit_task(generate_contacts, u64(random::random_device{}())); });

  // Submit 1 task reading random pixels from the matrix
  volatile i64 dummy_counter = 0;
  auto return_code = tpool.submit_task([&]() {
    auto rand_eng = random::PRNG(u64(random::random_device{}()));
    random::uniform_int_distribution<usize> idx_gen(0, m.ncols() - 1);

    while (!stop_sig) {
      dummy_counter += m.get(idx_gen(rand_eng), idx_gen(rand_eng));
    }
  });

  std::this_thread::sleep_for(std::chrono::seconds(15));
  // std::this_thread::sleep_for(std::chrono::seconds(3600));

  stop_sig = true;
  tpool.wait();

  u64 tot_contacts_expected = 0;
  for (auto& n : num_contacts_generated) {
    tot_contacts_expected += n.get();
  }
  [[maybe_unused]] const auto c = return_code.get();

  CHECK(m.get_tot_contacts() == tot_contacts_expected);
  CHECK(m.get_n_of_missed_updates() == 0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("ContactMatrixDense global locking TSAN", "[cmatrix][long]") {
  ContactMatrixDense<i64> m(100, 10);
  std::atomic<bool> stop_sig = false;

  // Increment a random pixel by 1 until stop_sign is flipped
  auto generate_contacts = [&](u64 seed) {
    auto rand_eng = random::PRNG(seed);
    random::uniform_int_distribution<usize> idx_gen(0, m.ncols() - 1);

    usize i = 0;
    while (!stop_sig) {
      m.increment(idx_gen(rand_eng), idx_gen(rand_eng));
      ++i;
    }
    return i;
  };

  const auto nthreads = std::max(u32(2), std::thread::hardware_concurrency());
  BS::light_thread_pool tpool(nthreads);

  // Submit nthreads - 1 tasks generating contacts
  std::vector<std::future<usize>> num_contacts_generated(nthreads - 1);
  std::generate(num_contacts_generated.begin(), num_contacts_generated.end(),
                [&]() { return tpool.submit_task(generate_contacts, u64(random::random_device{}())); });

  // Submit 1 task computing the calling get_tot_contacts() (which locks the entire matrix)
  u64 tot_contacts = 0;
  auto return_code = tpool.submit_task([&]() {
    auto rand_eng = random::PRNG(u64(random::random_device{}()));
    random::uniform_int_distribution<i64> sleep_gen(0, 100);
    while (!stop_sig) {
      std::this_thread::sleep_for(std::chrono::milliseconds(sleep_gen(rand_eng)));
      CHECK(m.get_n_of_missed_updates() == 0);
      const auto n = m.get_tot_contacts();
      CHECK(n >= tot_contacts);
      tot_contacts = n;
    }
  });

  std::this_thread::sleep_for(std::chrono::seconds(15));
  // std::this_thread::sleep_for(std::chrono::seconds(3600));

  stop_sig = true;
  tpool.wait();

  u64 tot_contacts_expected = 0;
  for (auto& n : num_contacts_generated) {
    tot_contacts_expected += n.get();
  }
  [[maybe_unused]] const auto c = return_code.get();

  CHECK(m.get_tot_contacts() == tot_contacts_expected);
  CHECK(m.get_n_of_missed_updates() == 0);
}
#endif

}  // namespace modle::cmatrix::test
