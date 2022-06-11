// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/stats/misc.hpp"

#include <absl/strings/str_split.h>
#include <fmt/format.h>

#include <boost/process.hpp>
#include <catch2/catch.hpp>
#include <vector>

#include "modle/common/numeric_utils.hpp"

namespace modle::test::stats {
using namespace modle::stats;
using namespace std::string_view_literals;

// clang-format off
// https://stackoverflow.com/a/43346070
constexpr auto GENERATE_REFERENCE_KERNEL_CMD{
    "import numpy as np\n"
    "l = {:d}\n"
    "sig = {:.16e}\n"
    "ax = np.linspace(-(l - 1) / 2., (l - 1) / 2., l)\n"
    "gauss = np.exp(-0.5 * np.square(ax) / np.square(sig))\n"
    "kernel = np.outer(gauss, gauss)\n"
    "kernel /= np.sum(kernel)\n"
    "print(\",\".join([str(n) for n in kernel.flatten()]))\n"sv};
// clang-format on

[[nodiscard]] static std::vector<double> generate_reference_kernel(const usize size,
                                                                   const double sigma) {
  boost::process::ipstream stdout_stream;
  auto c = boost::process::child(
      boost::process::search_path("python3").string(), "-c",
      fmt::format(FMT_STRING(GENERATE_REFERENCE_KERNEL_CMD), size, sigma),
      boost::process::std_in.close(), boost::process::std_out > stdout_stream);
  assert(c.running());

  std::string sbuff;
  std::getline(stdout_stream, sbuff);

  std::vector<double> buff;
  for (const auto& tok : absl::StrSplit(sbuff, ',')) {
    buff.push_back(utils::parse_numeric_or_throw<double>(tok));
  }

  assert(buff.size() == size * size);
  assert(!std::getline(stdout_stream, sbuff));
  return buff;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Gaussian kernel (SciPy)", "[stats][short]") {
  for (usize size = 3; size < 25; size += 2) {
    for (const auto sigma : {0.5, 1.0, 1.5, 2.5, 6.3, 10.0}) {
      const auto ref_kernel = generate_reference_kernel(size, sigma);
      const auto kernel = compute_gauss_kernel(size, sigma);
      REQUIRE(ref_kernel.size() == kernel.size());
      for (usize i = 0; i < kernel.size(); ++i) {
        CHECK(Approx(ref_kernel[i]) == kernel[i]);
      }
    }
  }
}

}  // namespace modle::test::stats
