// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/strings/str_split.h>  // for SplitIterator, Splitter, StrSplit
#include <fmt/format.h>              // for format

#include <BS_thread_pool.hpp>                       // for BS::thread_pool
#include <algorithm>                                // for generate, max
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset, dynamic_bitset<>::ref...
#include <boost/process.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>  // for path
#include <stdexcept>   // for runtime_error
#include <string>      // for string
#include <thread>      // for sleep_for
#include <utility>     // for pair, move
#include <vector>      // for vector, allocator

#include "modle/common/common.hpp"                // for u32
#include "modle/common/utils.hpp"                 // for parse_numeric_or_throw
#include "modle/compressed_io/compressed_io.hpp"  // for Reader
#include "modle/contacts.hpp"                     // for ContactMatrix

namespace modle::test::cmatrix {

[[maybe_unused]] static const std::filesystem::path& data_dir() {
  static const std::filesystem::path data_dir{"test/data/unit_tests"};
  return data_dir;
}

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

constexpr auto SCIPY_GAUSSIAN_DIFFERENCE_CMD{
    "#!/usr/bin/env python3\n"
    "from scipy.ndimage import gaussian_filter\n"
    "import numpy as np\n"
    "from sys import stdin\n"
    "if __name__ == \"__main__\":\n"
    "    shape = [{:d}, {:d}]\n"
    "    sigma1 = {:.16e}\n"
    "    sigma2 = {:.16e}\n"
    "    assert sigma1 < sigma2\n"
    "    trunc = {:.16e}\n"
    "    buff = stdin.read().replace(\"\\n\", \",\")[:-1]\n"
    "    m = np.fromstring(buff, sep=\",\", dtype=float).reshape(shape)\n"
    "    m1 = gaussian_filter(m, sigma1, truncate=trunc)\n"
    "    m2 = gaussian_filter(m, sigma2, truncate=trunc)\n"
    "    print(\",\".join([str(n) for n in (m1 - m2).flatten()]))\n"};
// clang-format on

template <class N>
static void write_cmatrix_to_stream(const ContactMatrix<N>& m, boost::process::opstream& s) {
  std::vector<contacts_t> buff(m.ncols());
  for (usize i = 0; i < m.ncols(); ++i) {
    buff.clear();
    for (usize j = 0; j < m.ncols(); ++j) {
      buff.push_back(m.get(i, j));
    }
    const auto sbuff = fmt::format(FMT_STRING("{}\n"), fmt::join(buff, ","));
    s.write(sbuff.data(), static_cast<std::streamsize>(sbuff.size()));
  }
  s.flush();
  s.pipe().close();
}

template <class N>
[[nodiscard]] static ContactMatrix<N> read_cmatrix_from_stream(const usize ncols, const usize nrows,
                                                               boost::process::ipstream& s) {
  std::string sbuff;
  std::getline(s, sbuff);

  std::vector<double> buff;
  for (const auto& tok : absl::StrSplit(sbuff, ',')) {
    buff.push_back(utils::parse_numeric_or_throw<double>(tok));
  }
  REQUIRE(buff.size() == nrows * ncols);
  ContactMatrix<double> m(nrows, ncols);
  for (usize i = 0; i < nrows; ++i) {
    for (auto j = i; j < ncols; ++j) {
      m.set(i, j, buff[(i * m.nrows()) + j]);
    }
  }
  return m;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix simple", "[cmatrix][short]") {
  ContactMatrix<> c(10, 100);
  CHECK(c.get(0, 0) == 0);
  c.increment(0, 0);
  CHECK(c.get(0, 0) == 1);
  c.increment(0, 0);
  CHECK(c.get(0, 0) == 2);
  c.subtract(0, 0, 2);
  CHECK(c.get(0, 0) == 0);
}

[[nodiscard]] inline std::vector<std::vector<u32>> load_matrix_from_file(
    const std::string& path_to_file, const std::string& sep = "\t") {
  std::vector<std::vector<u32>> m;
  compressed_io::Reader r(path_to_file);
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
TEST_CASE("CMatrix 10x200", "[cmatrix][medium]") {
  const auto input_file = data_dir() / "symm_matrix_200_10.tsv.gz";
  REQUIRE(std::filesystem::exists(input_file));
  const auto m1 = load_matrix_from_file(input_file.string());
  ContactMatrix<> m2(10, 200);
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
TEST_CASE("CMatrix in/decrement", "[cmatrix][short]") {
  ContactMatrix<> m(10, 20);
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
                          "Detected an out-of-bound read: attempt to access item at"));
    CHECK(m.get_n_of_missed_updates() == 1);
    CHECK(m.get_tot_contacts() == 1);

    CHECK_THROWS_WITH(m.decrement(25, 25),
                      Catch::Matchers::ContainsSubstring(
                          "Detected an out-of-bound read: attempt to access item at"));
    CHECK(m.get_n_of_missed_updates() == 1);
    CHECK(m.get_tot_contacts() == 1);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix get w/ block", "[cmatrix][short]") {
  ContactMatrix<u32> m1(100, 100);
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
TEST_CASE("CMatrix get w/ block small", "[cmatrix][short]") {
  const auto reference_file = data_dir() / "contacts_chr1_bs9_small.tsv";
  const auto input_file = data_dir() / "contacts_chr1_raw_small.tsv";

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
      CHECK(input_matrix.unsafe_get_block(i, j, block_size) == reference_matrix.get(i, j));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix get column", "[cmatrix][short]") {
  ContactMatrix<> c(10, 100);

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
TEST_CASE("CMatrix get row", "[cmatrix][short]") {
  ContactMatrix<> c(10, 100);

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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix blur (SciPy)", "[cmatrix][long]") {
  const auto reference_file = data_dir() / "cmatrix_002.tsv.gz";
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
    assert(py.running());

    write_cmatrix_to_stream(input_matrix, stdin_stream);
    return read_cmatrix_from_stream<double>(input_matrix.nrows(), input_matrix.ncols(),
                                            stdout_stream);
  };

  constexpr std::array<double, 3> sigmas{0.5, 1.0, 1.5};
  constexpr std::array<double, 3> cutoffs{3.0, 3.0, 3.0};

  for (usize i = 0; i < sigmas.size(); ++i) {
    const auto reference_matrix =
        compute_reference_matrix(input_matrix.ncols(), sigmas[i], cutoffs[i]);
    const auto blurred_matrix = input_matrix.blur(sigmas[i]);

    for (usize j = 4; j < input_matrix.nrows(); ++j) {
      for (auto k = j; k < input_matrix.ncols() - 4; ++k) {
        CHECK(Catch::Approx(reference_matrix.get(j, k)) == blurred_matrix.get(j, k));
      }
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix blur parallel (SciPy)", "[cmatrix][long]") {
  const auto reference_file = data_dir() / "cmatrix_002.tsv.gz";
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
    assert(py.running());

    write_cmatrix_to_stream(input_matrix, stdin_stream);
    return read_cmatrix_from_stream<double>(input_matrix.nrows(), input_matrix.ncols(),
                                            stdout_stream);
  };

  constexpr std::array<double, 3> sigmas{0.5, 1.0, 1.5};
  constexpr std::array<double, 3> cutoffs{3.0, 3.0, 3.0};

  BS::thread_pool tpool;
  for (usize i = 0; i < sigmas.size(); ++i) {
    const auto reference_matrix =
        compute_reference_matrix(input_matrix.ncols(), sigmas[i], cutoffs[i]);
    const auto blurred_matrix = input_matrix.blur(sigmas[i], 0.005, &tpool);

    for (usize j = 4; j < input_matrix.nrows(); ++j) {
      for (auto k = j; k < input_matrix.ncols() - 4; ++k) {
        CHECK(Catch::Approx(reference_matrix.get(j, k)) == blurred_matrix.get(j, k));
      }
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix difference of gaussians (SciPy)", "[cmatrix][long]") {
  const auto reference_file = data_dir() / "cmatrix_002.tsv.gz";
  const auto input_matrix = [&]() {
    ContactMatrix<> m;
    m.unsafe_import_from_txt(reference_file);
    return m;
  }();

  auto compute_reference_matrix = [&input_matrix](auto shape, auto sigma1, auto sigma2,
                                                  auto trunc) {
    boost::process::ipstream stdout_stream;
    boost::process::opstream stdin_stream;
    auto py = boost::process::child(
        boost::process::search_path("python3").string(), "-c",
        fmt::format(FMT_STRING(SCIPY_GAUSSIAN_DIFFERENCE_CMD), shape, shape, sigma1, sigma2, trunc),
        boost::process::std_in<stdin_stream, boost::process::std_out> stdout_stream);
    assert(py.running());

    write_cmatrix_to_stream(input_matrix, stdin_stream);
    return read_cmatrix_from_stream<double>(input_matrix.nrows(), input_matrix.ncols(),
                                            stdout_stream);
  };

  const double sigma1_ = 1.0;
  const double sigma2_ = 1.6;
  const double trunc_ = 3.0;
  const auto reference_matrix =
      compute_reference_matrix(input_matrix.ncols(), sigma1_, sigma2_, trunc_);
  const auto gauss_diff_matrix = input_matrix.gaussian_diff(sigma1_, sigma2_);

  const auto m1 = input_matrix.blur(sigma1_);
  const auto m2 = input_matrix.blur(sigma2_);

  for (usize j = 4; j < input_matrix.nrows(); ++j) {
    for (auto k = j; k < input_matrix.ncols() - 4; ++k) {
      CHECK(Catch::Approx(reference_matrix.get(j, k)) == gauss_diff_matrix.get(j, k));
      CHECK(Catch::Approx(m1.get(j, k) - m2.get(j, k)) == gauss_diff_matrix.get(j, k));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix difference of gaussians - parallel (SciPy)", "[cmatrix][long]") {
  const auto reference_file = data_dir() / "cmatrix_002.tsv.gz";
  const auto input_matrix = [&]() {
    ContactMatrix<> m;
    m.unsafe_import_from_txt(reference_file);
    return m;
  }();

  auto compute_reference_matrix = [&input_matrix](auto shape, auto sigma1, auto sigma2,
                                                  auto trunc) {
    boost::process::ipstream stdout_stream;
    boost::process::opstream stdin_stream;
    auto py = boost::process::child(
        boost::process::search_path("python3").string(), "-c",
        fmt::format(FMT_STRING(SCIPY_GAUSSIAN_DIFFERENCE_CMD), shape, shape, sigma1, sigma2, trunc),
        boost::process::std_in<stdin_stream, boost::process::std_out> stdout_stream);
    assert(py.running());

    write_cmatrix_to_stream(input_matrix, stdin_stream);
    return read_cmatrix_from_stream<double>(input_matrix.nrows(), input_matrix.ncols(),
                                            stdout_stream);
  };

  const double sigma1_ = 1.0;
  const double sigma2_ = 1.6;
  const double trunc_ = 3.0;

  BS::thread_pool tpool;
  const auto reference_matrix =
      compute_reference_matrix(input_matrix.ncols(), sigma1_, sigma2_, trunc_);
  const auto gauss_diff_matrix =
      input_matrix.gaussian_diff(sigma1_, sigma2_, std::numeric_limits<double>::lowest(),
                                 (std::numeric_limits<double>::max)(), &tpool);

  const auto m1 = input_matrix.blur(sigma1_);
  const auto m2 = input_matrix.blur(sigma2_);

  for (usize j = 4; j < input_matrix.nrows(); ++j) {
    for (auto k = j; k < input_matrix.ncols() - 4; ++k) {
      CHECK(Catch::Approx(reference_matrix.get(j, k)) == gauss_diff_matrix.get(j, k));
      CHECK(Catch::Approx(m1.get(j, k) - m2.get(j, k)) == gauss_diff_matrix.get(j, k));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix test get_nnz", "[cmatrix][short]") {
  ContactMatrix<> m(10, 10);
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
TEST_CASE("CMatrix test get_tot_contacts", "[cmatrix][short]") {
  ContactMatrix<> m(10, 10);
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
TEST_CASE("CMatrix pixel locking TSAN", "[cmatrix][long]") {
  ContactMatrix<i64> m(100, 10);
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
  BS::thread_pool tpool(nthreads);

  // Submit nthreads - 1 tasks generating contacts
  std::vector<std::future<usize>> num_contacts_generated(nthreads - 1);
  std::generate(num_contacts_generated.begin(), num_contacts_generated.end(),
                [&]() { return tpool.submit(generate_contacts, u64(random::random_device{}())); });

  // Submit 1 task reading random pixels from the matrix
  volatile i64 dummy_counter = 0;
  auto return_code = tpool.submit([&]() {
    auto rand_eng = random::PRNG(u64(random::random_device{}()));
    random::uniform_int_distribution<usize> idx_gen(0, m.ncols() - 1);

    while (!stop_sig) {
      dummy_counter += m.get(idx_gen(rand_eng), idx_gen(rand_eng));
    }
  });

  std::this_thread::sleep_for(std::chrono::seconds(15));
  // std::this_thread::sleep_for(std::chrono::seconds(3600));

  stop_sig = true;
  tpool.wait_for_tasks();

  u64 tot_contacts_expected = 0;
  for (auto& n : num_contacts_generated) {
    tot_contacts_expected += n.get();
  }
  [[maybe_unused]] const auto c = return_code.get();

  CHECK(m.get_tot_contacts() == tot_contacts_expected);
  CHECK(m.get_n_of_missed_updates() == 0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("CMatrix global locking TSAN", "[cmatrix][long]") {
  ContactMatrix<i64> m(100, 10);
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
  BS::thread_pool tpool(nthreads);

  // Submit nthreads - 1 tasks generating contacts
  std::vector<std::future<usize>> num_contacts_generated(nthreads - 1);
  std::generate(num_contacts_generated.begin(), num_contacts_generated.end(),
                [&]() { return tpool.submit(generate_contacts, u64(random::random_device{}())); });

  // Submit 1 task computing the calling get_tot_contacts() (which locks the entire matrix)
  u64 tot_contacts = 0;
  auto return_code = tpool.submit([&]() {
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
  tpool.wait_for_tasks();

  u64 tot_contacts_expected = 0;
  for (auto& n : num_contacts_generated) {
    tot_contacts_expected += n.get();
  }
  [[maybe_unused]] const auto c = return_code.get();

  CHECK(m.get_tot_contacts() == tot_contacts_expected);
  CHECK(m.get_n_of_missed_updates() == 0);
}
#endif

}  // namespace modle::test::cmatrix
