#include "modle/correlation.hpp"  // for compute_pearson, compute_pearson_significance, compute_...

#include <boost/filesystem/path.hpp>  // for path
#include <catch2/catch.hpp>  // for Approx, operator==, AssertionHandler, operator""_catch_sr
#include <cstdint>           // for uint32_t
#include <vector>            // for vector, allocator

#include "./common.hpp"  // for generate_random_vect, corr_scipy
#include "modle/common/random.hpp"
#include "modle/common/smartdir.hpp"  // IWYU pragma: keep

namespace modle::test {
constexpr auto cleanup_on_exit{true};     // Useful for debugging
const SmartDir testdir{cleanup_on_exit};  // NOLINT Using auto here upsets GCC8
}  // namespace modle::test

namespace modle::test::correlation {
using namespace modle::correlation;

TEST_CASE("Corr. test: Pearson wo ties", "[correlation][pearson][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  const auto pcc = compute_pearson(v1, v2);
  const auto pv = compute_pearson_significance(pcc, v1.size());
  CHECK(Approx(pcc).margin(0) == -0.033621194725622014);  // NOLINT
  CHECK(Approx(pv).margin(0) == 0.926536715854247);       // NOLINT
}

TEST_CASE("Corr. test: Pearson w ties", "[correlation][pearson][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  const auto pcc = compute_pearson(v1, v2);
  const auto pv = compute_pearson_significance(pcc, v1.size());
  CHECK(Approx(pcc).margin(0) == 0.16426413174421572);  // NOLINT
  CHECK(Approx(pv).margin(0) == 0.6502118872600098);    // NOLINT
}

TEST_CASE("Corr. test: Pearson Scipy", "[correlation][pearson][short]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  auto v1 = generate_random_vect(rnd_eng, 1'000, 0, 15'000);  // NOLINT
  auto v2 = generate_random_vect(rnd_eng, 1'000, 0, 15'000);  // NOLINT
  const auto pcc = compute_pearson(v1, v2);
  const auto pv = compute_pearson_significance(pcc, v1.size());
  const auto& [pcc_py, pv_py] = corr_scipy(v1, v2, "pearsonr", testdir());
  CHECK(Approx(pcc).margin(0) == pcc_py);
  CHECK(Approx(pv).margin(0) == pv_py);
}

TEST_CASE("Corr. test: Pearson Scipy long", "[correlation][pearson][long]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  for (auto i = 0UL; i < 100UL; ++i) {                          // NOLINT
    auto v1 = generate_random_vect(rnd_eng, 1'000, 0, 15'000);  // NOLINT
    auto v2 = generate_random_vect(rnd_eng, 1'000, 0, 15'000);  // NOLINT
    const auto pcc = compute_pearson(v1, v2);
    const auto pv = compute_pearson_significance(pcc, v1.size());
    const auto& [pcc_py, pv_py] = corr_scipy(v1, v2, "pearsonr", testdir());
    CHECK(Approx(pcc).margin(0) == pcc_py);
    CHECK(Approx(pv).margin(0) == pv_py);
  }
}

TEST_CASE("Corr. test: Pearson Scipy long vect.", "[correlation][pearson][medium]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  auto v1 = generate_random_vect(rnd_eng, 1'000'000, 0, 15'000);  // NOLINT
  auto v2 = generate_random_vect(rnd_eng, 1'000'000, 0, 15'000);  // NOLINT
  const auto pcc = compute_pearson(v1, v2);
  const auto pv = compute_pearson_significance(pcc, v1.size());
  const auto& [pcc_py, pv_py] = corr_scipy(v1, v2, "pearsonr", testdir());
  CHECK(Approx(pcc).margin(0) == pcc_py);
  CHECK(Approx(pv).margin(0) == pv_py);
}

TEST_CASE("Corr. test: Spearman wo ties", "[correlation][spearman][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  const auto rho = compute_spearman(v1, v2);
  const auto pv = compute_spearman_significance(rho, v1.size());
  CHECK(Approx(rho).margin(0) == -0.16363636363636364);  // NOLINT
  CHECK(Approx(pv).margin(0) == 0.6514773427962428);     // NOLINT
}

TEST_CASE("Corr. test: Spearman w ties", "[correlation][spearman][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 47, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  const auto rho = compute_spearman(v1, v2);
  const auto pv = compute_spearman_significance(rho, v1.size());
  CHECK(Approx(rho).margin(0) == 0.024316221747202587);  // NOLINT
  CHECK(Approx(pv).margin(0) == 0.9468397049085097);     // NOLINT
}

TEST_CASE("Corr. test: Spearman Scipy", "[correlation][spearman][short]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  auto v1 = generate_random_vect(rnd_eng, 1'000, 0, 15'000);  // NOLINT
  auto v2 = generate_random_vect(rnd_eng, 1'000, 0, 15'000);  // NOLINT
  const auto rho = compute_spearman(v1, v2);
  const auto pv = compute_spearman_significance(rho, v1.size());
  const auto& [rho_py, pv_py] = corr_scipy(v1, v2, "spearmanr", testdir());
  CHECK(Approx(rho).margin(0) == rho_py);
  CHECK(Approx(pv).margin(0) == pv_py);
}

TEST_CASE("Corr. test: Spearman Scipy long", "[correlation][spearman][long]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  for (auto i = 0UL; i < 100; ++i) {                            // NOLINT
    auto v1 = generate_random_vect(rnd_eng, 1'000, 0, 15'000);  // NOLINT
    auto v2 = generate_random_vect(rnd_eng, 1'000, 0, 15'000);  // NOLINT
    const auto rho = compute_spearman(v1, v2);
    const auto pv = compute_spearman_significance(rho, v1.size());
    const auto& [rho_py, pv_py] = corr_scipy(v1, v2, "spearmanr", testdir());
    CHECK(Approx(rho).margin(0) == rho_py);
    CHECK(Approx(pv).margin(0) == pv_py);
  }
}

TEST_CASE("Corr. test: Spearman Scipy long vect.", "[correlation][spearman][medium]") {
  std::mt19937 rnd_eng{std::random_device{}()};
  auto v1 = generate_random_vect(rnd_eng, 1'000'000, 0, 15'000);  // NOLINT
  auto v2 = generate_random_vect(rnd_eng, 1'000'000, 0, 15'000);  // NOLINT
  const auto rho = compute_spearman(v1, v2);
  const auto pv = compute_spearman_significance(rho, v1.size());
  const auto& [rho_py, pv_py] = corr_scipy(v1, v2, "spearmanr", testdir());
  CHECK(Approx(rho).margin(0) == rho_py);
  CHECK(Approx(pv).margin(0) == pv_py);
}

TEST_CASE("SED", "[correlation][sed][short]") {
  std::vector<uint32_t> v1{17, 86, 60, 77, 47, 3, 70, 87, 88, 92};   // NOLINT
  std::vector<uint32_t> v2{70, 29, 85, 61, 80, 34, 60, 31, 73, 66};  // NOLINT
  const auto sed = compute_sed(v1, v2);
  CHECK(Approx(sed).margin(0) == 13125.999999999998);  // NOLINT
}

}  // namespace modle::test::correlation
