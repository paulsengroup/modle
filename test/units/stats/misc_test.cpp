// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/stats/misc.hpp"

#include <absl/strings/str_split.h>
#include <fmt/format.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <vector>

#include "modle/common/numeric_utils.hpp"

namespace modle::stats::test {

[[maybe_unused]] static const std::filesystem::path& data_dir() {
  static const std::filesystem::path data_dir{"test/data/unit_tests"};
  return data_dir;
}

[[nodiscard]] static std::vector<double> import_reference_kernel(const usize radius,
                                                                 const double sigma) {
  const auto sbuff = [&]() {
    const auto path = data_dir() / std::filesystem::path{"reference_gaussian_kernels"} /
                      fmt::format(FMT_STRING("gaussian_kernel_{}_{:.1f}.csv"), radius, sigma);

    REQUIRE(std::filesystem::exists(path));
    std::ifstream f(path);

    std::stringstream buffer;
    buffer << f.rdbuf();
    return buffer.str();
  }();

  std::vector<double> buff;
  for (const auto& tok : absl::StrSplit(sbuff, absl::ByAnyChar("\t,\n"))) {
    if (!tok.empty()) {
      buff.push_back(utils::parse_numeric_or_throw<double>(tok));
    }
  }

  const auto kernel_size1d = radius * 2 + 1;
  REQUIRE(buff.size() == kernel_size1d * kernel_size1d);
  return buff;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Gaussian kernel", "[stats][short]") {
  for (usize radius = 1; radius < 15; ++radius) {
    for (const auto sigma : {0.5, 1.0, 1.5, 2.5, 6.3, 10.0}) {
      const auto ref_kernel = import_reference_kernel(radius, sigma);
      const auto kernel = compute_gauss_kernel2d(radius, sigma);

      REQUIRE(ref_kernel.size() == kernel.size());
      for (usize i = 0; i < kernel.size(); ++i) {
        CHECK(Catch::Approx(ref_kernel[i]) == kernel[i]);
      }
    }
  }
}

}  // namespace modle::stats::test
