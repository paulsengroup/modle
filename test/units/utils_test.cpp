#include <cstdint>
#include <vector>

#include "gtest/gtest.h"
#include "modle/impl/correlation_utils.hpp"

TEST(computecpp_test_suite, sort_vector_by_idx) {
  const std::vector<uint32_t> v{10, 5, 67, 3, 60, 45, 49, 1000};
  const std::vector<uint32_t> expected{3, 1, 0, 5, 6, 4, 2, 7};
  const auto vi = modle::sort_vector_by_idx(v.begin(), v.end());
  ASSERT_EQ(vi.size(), expected.size());
  for (auto i = 0UL; i < expected.size(); ++i) {
    EXPECT_EQ(vi[i], expected[i]);
  }
}

TEST(computecpp_test_suite, compute_element_ranks) {
  const std::vector<uint32_t> v{10, 5, 67, 3, 60, 45, 49, 1000};
  const std::vector<uint32_t> expected{2, 1, 6, 0, 5, 3, 4, 7};
  const auto vi = modle::compute_element_ranks(v.begin(), v.end());
  ASSERT_EQ(vi.size(), expected.size());
  for (auto i = 0UL; i < expected.size(); ++i) {
    EXPECT_EQ(vi[i], expected[i]);
  }
}
TEST(computecpp_test_suite, compute_element_ranks_tie) {
  const std::vector<uint32_t> v{10, 5, 67, 3, 67, 45, 49, 1000};
  const std::vector<double> expected{2, 1, 5.5, 0, 5.5, 3, 4, 7};
  const auto vi = modle::compute_element_ranks(v.begin(), v.end());
  ASSERT_EQ(vi.size(), expected.size());
  for (auto i = 0UL; i < expected.size(); ++i) {
    EXPECT_EQ(vi[i], expected[i]);
  }
}
