#include "modle/correlation.hpp"

#include <cstdint>
#include <utility>
#include <vector>

#include "catch2/catch.hpp"
#include "common.hpp"

namespace modle::correlation::test {

TEST_CASE("Corr. test: Spearman wo ties", "[correlation][spearman][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  auto ct = CorrelationTest(v1, v2);
  const auto& [rho, pv] = ct.compute_spearman();
  CHECK(Approx(rho).margin(0) == -0.16363636363636364);  // NOLINT
  CHECK(Approx(pv).margin(0) == 0.6514773427962428);     // NOLINT
}

TEST_CASE("Corr. test: Spearman w ties", "[correlation][spearman][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  auto ct = CorrelationTest(v1, v2);
  const auto& [rho, pv] = ct.compute_spearman();
  CHECK(Approx(rho).margin(0) == 0.024316221747202587);  // NOLINT
  CHECK(Approx(pv).margin(0) == 0.9468397049085097);     // NOLINT
}

TEST_CASE("Corr. test: Spearman Scipy", "[correlation][spearman][short]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  auto v1 = generate_random_vect(rnd_eng);
  auto v2 = generate_random_vect(rnd_eng);
  auto ct = CorrelationTest(v1, v2);
  const auto& [rho, pv] = ct.compute_spearman();
  const auto& [rho_py, pv_py] = corr_scipy(v1, v2, "spearmanr");
  CHECK(Approx(rho).margin(0) == rho_py);
  CHECK(Approx(pv).margin(0) == pv_py);
}

TEST_CASE("Corr. test: Spearman Scipy long", "[correlation][spearman][medium]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  for (auto i = 0UL; i < 100; ++i) {  // NOLINT
    auto v1 = generate_random_vect(rnd_eng);
    auto v2 = generate_random_vect(rnd_eng);
    auto ct = CorrelationTest(v1, v2);
    const auto& [rho, pv] = ct.compute_spearman();
    const auto& [rho_py, pv_py] = corr_scipy(v1, v2, "spearmanr");
    CHECK(Approx(rho).margin(0) == rho_py);
    CHECK(Approx(pv).margin(0) == pv_py);
  }
}

TEST_CASE("Corr. test: Spearman Scipy long vect.", "[correlation][spearman][medium]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  auto v1 = generate_random_vect(rnd_eng, 1'000'000);  // NOLINT
  auto v2 = generate_random_vect(rnd_eng, 1'000'000);  // NOLINT
  auto ct = CorrelationTest(v1, v2);
  const auto& [rho, pv] = ct.compute_spearman();
  const auto& [rho_py, pv_py] = corr_scipy(v1, v2, "spearmanr");
  CHECK(Approx(rho).margin(0) == rho_py);
  CHECK(Approx(pv).margin(0) == pv_py);
}

TEST_CASE("Corr. test: Kendall wo ties", "[correlation][kendall][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  auto ct = CorrelationTest(v1, v2);
  const auto& [tau, pv] = ct.compute_kendall();
  CHECK(Approx(tau).margin(0) == -0.06666666666666667);  // NOLINT
  CHECK(Approx(pv).margin(0) == 0.8618005952380953);     // NOLINT
}

TEST_CASE("Corr. test: Kendall w ties", "[correlation][kendall][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  auto ct = CorrelationTest(v1, v2);
  const auto& [tau, pv] = ct.compute_kendall();
  CHECK(Approx(tau).margin(0) == 0.04494665749754947);  // NOLINT
  CHECK(Approx(pv).margin(0) == 0.8574624419592412);    // NOLINT
}

TEST_CASE("Corr. test: Kendall Scipy", "[correlation][kendall][short]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  auto v1 = generate_random_vect(rnd_eng);
  auto v2 = generate_random_vect(rnd_eng);
  auto ct = CorrelationTest(v1, v2);
  const auto& [tau, pv] = ct.compute_kendall();
  const auto& [tau_py, pv_py] = corr_scipy(v1, v2, "kendalltau");
  CHECK(Approx(tau).margin(0) == tau_py);
  CHECK(Approx(pv).margin(0) == pv_py);
}

TEST_CASE("Corr. test: Kendall Scipy long", "[correlation][kendall][medium]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  for (auto i = 0UL; i < 100; ++i) {  // NOLINT
    auto v1 = generate_random_vect(rnd_eng);
    auto v2 = generate_random_vect(rnd_eng);
    auto ct = CorrelationTest(v1, v2);
    const auto& [tau, pv] = ct.compute_kendall();
    const auto& [tau_py, pv_py] = corr_scipy(v1, v2, "kendalltau");
    CHECK(Approx(tau).margin(0) == tau_py);
    CHECK(Approx(pv).margin(0) == pv_py);
  }
}

TEST_CASE("Corr. test: Kendall Scipy long vect.", "[correlation][kendall][medium]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  auto v1 = generate_random_vect(rnd_eng, 1'000'000);  // NOLINT
  auto v2 = generate_random_vect(rnd_eng, 1'000'000);  // NOLINT
  auto ct = CorrelationTest(v1, v2);
  const auto& [tau, pv] = ct.compute_kendall();
  const auto& [tau_py, pv_py] = corr_scipy(v1, v2, "kendalltau");
  CHECK(Approx(tau).margin(0) == tau_py);
  CHECK(Approx(pv).margin(0) == pv_py);
}

}  // namespace modle::correlation::test