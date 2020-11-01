#include "modle/kendall.hpp"

#include "common.hpp"
#include "gtest/gtest.h"

TEST(modle_test_suite, kendall) {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto& [tau, pv] = modle::compute_kendall(v1, v2);
  EXPECT_FLOAT_EQ(tau, -0.06666666666666667);
  EXPECT_FLOAT_EQ(pv, 0.8618005952380953);
}

TEST(modle_test_suite, kendall_w_tie) {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto& [tau, pv] = modle::compute_kendall(v1, v2);
  EXPECT_FLOAT_EQ(tau, 0.04494665749754947);
  EXPECT_FLOAT_EQ(pv, 0.8574624419592412);
}

TEST(modle_test_suite, kendall_scipy) {
  std::random_device rnd_device;
  std::mt19937 rnd_eng{rnd_device()};
  {
    auto v1 = generate_random_vect(rnd_eng);
    auto v2 = generate_random_vect(rnd_eng);
    const auto& [tau, pv] = modle::compute_kendall(v1, v2);
    const auto& [tau_py, pv_py] = corr_scipy(v1, v2, "kendalltau");
    EXPECT_FLOAT_EQ(tau, tau_py);
    EXPECT_FLOAT_EQ(pv, pv_py);
  }
}

TEST(modle_test_suite, kendall_scipy_long) {
  std::random_device rnd_device;
  std::mt19937 rnd_eng{rnd_device()};
  for (auto i = 0UL; i < 100; ++i) {
    auto v1 = generate_random_vect(rnd_eng);
    auto v2 = generate_random_vect(rnd_eng);
    const auto& [tau, pv] = modle::compute_kendall(v1, v2);
    const auto& [tau_py, pv_py] = corr_scipy(v1, v2, "kendalltau");
    EXPECT_FLOAT_EQ(tau, tau_py);
    EXPECT_FLOAT_EQ(pv, pv_py);
  }
}

TEST(modle_test_suite, kendall_scipy_long_vector) {
  std::random_device rnd_device;
  std::mt19937 rnd_eng{rnd_device()};
  auto v1 = generate_random_vect(rnd_eng, 1'000'000);
  auto v2 = generate_random_vect(rnd_eng, 1'000'000);
  const auto& [tau, pv] = modle::compute_kendall(v1, v2);
  const auto& [tau_py, pv_py] = corr_scipy(v1, v2, "kendalltau");
  EXPECT_FLOAT_EQ(tau, tau_py);
  EXPECT_FLOAT_EQ(pv, pv_py);
}
