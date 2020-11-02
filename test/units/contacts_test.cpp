#include "modle/contacts.hpp"

#include <cstdint>
#include <fstream>
#include <vector>

#include "absl/strings/str_format.h"
#include "absl/strings/str_split.h"
#include "gtest/gtest.h"

namespace modle {

TEST(modle_test_suite, contacts_simple_test) {
  ContactMatrix<uint32_t> c(10, 100);
  EXPECT_EQ(c.get(0, 0), 0);
  c.increment(0, 0);
  EXPECT_EQ(c.get(0, 0), 1);
  c.increment(0, 0);
  EXPECT_EQ(c.get(0, 0), 2);
  c.decrement(0, 0, 2);
  EXPECT_EQ(c.get(0, 0), 0);
}

std::vector<std::vector<uint32_t>> load_matrix_from_file(const std::string& path_to_file,
                                                         const std::string& sep = "\t") {
  std::vector<std::vector<uint32_t>> m;
  std::ifstream f(path_to_file);
  if (f.is_open()) {
    std::string line;
    std::string buff;
    while (std::getline(f, line)) {
      std::vector<uint32_t> v;

      for (const auto& tok : absl::StrSplit(line, sep)) {
        buff = tok;
        v.push_back(std::stoul(buff));
      }
      m.emplace_back(std::move(v));
    }
  } else {
    throw std::runtime_error(absl::StrFormat("Unable to open file '%s'.", path_to_file));
  }
  return m;
}

TEST(modle_test_suite, contacts_small_matrix_test) {
  const auto m1 = load_matrix_from_file("data/symm_matrix_200_10.tsv");
  uint32_t size = 200, diagonal_width = 10;
  ContactMatrix<uint32_t> m2(diagonal_width, size);
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
      EXPECT_EQ(m1[i][j], m3[i][j])
          << "i=" << i << "; j=" << j << "; expected=" << m1[i][j] << "; got=" << m3[i][j] << ";";
    }
  }
}

ContactMatrix<uint32_t> build_matrix_from_file(const std::string& path_to_file, uint32_t size,
                                               uint32_t diagonal_width,
                                               const std::string& sep = "\t") {
  ContactMatrix<uint32_t> m(diagonal_width, size);
  std::ifstream f(path_to_file);
  if (f.is_open()) {
    std::string line;
    std::string buff;
    uint32_t i = 0, j;
    while (std::getline(f, line)) {
      j = 0;
      for (const auto& tok : absl::StrSplit(line, '\t')) {
        //        absl::FPrintF(stderr, "%lu;%lu='%s'\n", i, j, tok);
        if (tok != "0") {
          buff = tok;
          m.set(i, j, std::stoul(buff));
        }
        ++j;
      }
      ++i;
    }
  } else {
    throw std::runtime_error(absl::StrFormat("Unable to open file '%s'.", path_to_file));
  }
  return m;
}

void test_with_huge_matrix(const std::string& path_to_file, uint32_t size, uint32_t diagonal_width,
                           const std::string& sep = "\t") {
  const auto m = build_matrix_from_file(path_to_file, size, diagonal_width, sep);
  std::ifstream f(path_to_file);

  if (f.is_open()) {
    std::string line;
    std::string buff;
    uint32_t i = 0, j;
    while (std::getline(f, line)) {
      j = 0;
      for (const auto& tok : absl::StrSplit(line, '\t')) {
        if (tok != "0") {
          buff = tok;
          EXPECT_EQ(m.get(i, j), std::stoul(buff))
              << "i=" << i << "; j=" << j << "; expected=" << buff << "; got=" << m.get(i, j)
              << ";";
        }
        ++j;
      }
      ++i;
    }
  } else {
    throw std::runtime_error(absl::StrFormat("Unable to open file '%s'.", path_to_file));
  }
}

TEST(modle_test_suite, contacts_large_matrix_test) {
  uint32_t size = 5'000, diagonal_width = 50;
  std::string path_to_file = "data/symm_matrix_5000_50.tsv";
  test_with_huge_matrix(path_to_file, size, diagonal_width);
}

TEST(modle_test_suite, contacts_huge_matrix_test) {
  uint32_t size = 50'000, diagonal_width = 50;
  std::string path_to_file = "data/symm_matrix_50000_50.tsv";
  test_with_huge_matrix(path_to_file, size, diagonal_width);
}

}  // namespace modle