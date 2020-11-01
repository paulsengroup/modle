#include "modle/spearman.hpp"

#include "common.hpp"
#include "gtest/gtest.h"

TEST(modle_test_suite, spearman) {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto& [rho, pv] = modle::compute_spearman(v1, v2);
  EXPECT_FLOAT_EQ(rho, -0.16363636363636364);
  EXPECT_FLOAT_EQ(pv, 0.6514773427962428);
}

TEST(modle_test_suite, spearman_w_tie) {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};
  const auto& [rho, pv] = modle::compute_spearman(v1, v2);
  EXPECT_FLOAT_EQ(rho, 0.024316221747202587);
  EXPECT_FLOAT_EQ(pv, 0.9468397049085097);
}

TEST(modle_test_suite, spearman_scipy) {
  std::random_device rnd_device;
  std::mt19937 rnd_eng{rnd_device()};
  {
    auto v1 = generate_random_vect(rnd_eng);
    auto v2 = generate_random_vect(rnd_eng);
    const auto& [rho, pv] = modle::compute_spearman(v1, v2);
    const auto& [rho_py, pv_py] = corr_scipy(v1, v2, "spearmanr");
    EXPECT_FLOAT_EQ(rho, rho_py);
    EXPECT_FLOAT_EQ(pv, pv_py);
  }
}

TEST(modle_test_suite, spearman_scipy_long) {
  std::random_device rnd_device;
  std::mt19937 rnd_eng{rnd_device()};
  for (auto i = 0UL; i < 100; ++i) {
    auto v1 = generate_random_vect(rnd_eng);
    auto v2 = generate_random_vect(rnd_eng);
    const auto& [rho, pv] = modle::compute_spearman(v1, v2);
    const auto& [rho_py, pv_py] = corr_scipy(v1, v2, "spearmanr");
    EXPECT_FLOAT_EQ(rho, rho_py);
    EXPECT_FLOAT_EQ(pv, pv_py);
  }
}

TEST(modle_test_suite, spearman_scipy_long_vector) {
  std::random_device rnd_device;
  std::mt19937 rnd_eng{rnd_device()};
  auto v1 = generate_random_vect(rnd_eng, 1'000'000);
  auto v2 = generate_random_vect(rnd_eng, 1'000'000);
  const auto& [rho, pv] = modle::compute_spearman(v1, v2);
  const auto& [rho_py, pv_py] = corr_scipy(v1, v2, "spearmanr");
  EXPECT_FLOAT_EQ(rho, rho_py);
  EXPECT_FLOAT_EQ(pv, pv_py);
}
