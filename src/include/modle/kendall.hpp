#pragma once

#include <algorithm>
#include <atomic>
#include <boost/math/special_functions/erf.hpp>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "modle/utils.hpp"

namespace modle {

double compute_kendallb_significance(uint32_t nc, uint32_t nd, uint32_t size, uint32_t tie1,
                                     uint32_t tie2) {
  const uint32_t tot = size * (size - 1) / 2;
  if ((tie1 == 0 && tie2 == 0) && (size <= 33 || std::min(nd, tot - nd) <= 1)) {
    // Exact p-value, see p. 68 of Maurice G. Kendall, "Rank Correlation Methods" (4th Edition),
    // Charles Griffin & Co., 1970.
    auto c = std::min(nd, tot - nd);
    assert(2 * c < size * (size - 1));
    if (size == 1 || size == 2 || (2 * c) == tot)
      return 1.0;
    else if (c == 0 || c == 1)
      return 2.0 / std::tgamma(size - (c == 1));  // tgamma approx factorial
    else {
      std::vector<double> v_new(c + 1, 0.0);
      v_new[0] = 1;
      v_new[1] = 1;
      for (auto i = 3U; i <= size; ++i) {
        auto v_old = v_new;
        for (auto j = 1U; j < std::min(i, c + 1); ++j) {
          v_new[j] += v_new[j - 1];
        }
        for (auto j = i; j <= c; ++j) {
          v_new[j] += v_new[j - 1] - v_old[j - i];
        }
      }
      return 0.2 * std::reduce(v_new.begin(), v_new.end()) /
             std::tgamma(size);  // looks like pv is multiplied by 10 for some reason
    }
  }
  uint64_t x0 = tie1 * (tie1 - 1) * (tie1 - 2);
  uint64_t y0 = tie2 * (tie2 - 1) * (tie2 - 2);
  uint64_t x1 = tie1 * (tie1 - 1) * (2 * tie1 + 5);
  uint64_t y1 = tie2 * (tie2 - 1) * (2 * tie2 + 5);
  double var = (size * (size - 1) * (2.0 * size + 5) - x1 - y1) / 18.0 +
               (2.0 * tie1 * tie2) / (size * (size - 1)) +
               x0 * y0 / (9.0 * size * (size - 1) * (size - 2));
  return boost::math::erfc<double>(std::abs<double>(static_cast<int64_t>(nc) - nd) /
                                   std::sqrt(var) / std::sqrt(2.0));
}

template <typename Iterator,
          typename = typename std::enable_if<
              std::is_base_of<std::random_access_iterator_tag,
                              typename std::iterator_traits<Iterator>::iterator_category>::value>>
std::pair<double, double> compute_kendall_b_n2(Iterator v1_b, Iterator v1_e, Iterator v2_b,
                                               Iterator v2_e) {
  assert(std::distance(v1_b, v1_e) == std::distance(v2_b, v2_e));
  const uint32_t size = std::distance(v1_b, v1_e);
  const uint32_t n_pairs = size * (size - 1) / 2;
  double tau = 0.0;

  std::atomic<uint32_t> tie1 = 0 /* n of ties in v1 */, tie2 = 0 /* n of ties in v2 */;
  std::atomic<int64_t> nc = 0 /* n of concordant ranks */, nd = 0 /* n of discordant ranks */;

  std::vector<uint32_t> idx(size);
  std::iota(idx.begin(), idx.end(), 0);

  auto ktb_fx = [&](uint32_t i) {
    uint32_t m1_local = 0 /* n of ties in v1 */, m2_local = 0 /* n of ties in v2 */;
    int64_t nc_local = 0 /* n of concordant ranks */, nd_local = 0 /* n of discordant ranks */;
    auto v1i = *(v1_b + i);  // Equivalent to v1[i]
    auto v2i = *(v2_b + i);  // Equivalent to v2[i]
    for (auto j = i + 1; j < size; ++j) {
      auto v1j = *(v1_b + j);  // Equivalent to v1[j]
      auto v2j = *(v2_b + j);  // Equivalent to v2[j]
      if (v2i > v2j) {
        nc_local += v1i > v1j;
        nd_local += v1i < v1j;
        m1_local += v1i == v1j;
      } else if (v2i < v2j) {
        nd_local += v1i > v1j;
        nc_local += v1i < v1j;
        m1_local += v1i == v1j;
      } else {
        ++m2_local;
        m1_local += v1i == v1j;
      }
    }
    nc.fetch_add(nc_local, std::memory_order_relaxed);
    nd.fetch_add(nd_local, std::memory_order_relaxed);
    tie1.fetch_add(m1_local, std::memory_order_relaxed);
    tie2.fetch_add(m2_local, std::memory_order_relaxed);
  };

  if (size < 1000)
    std::for_each(std::execution::seq, idx.begin(), idx.end(), ktb_fx);
  else
    std::for_each(std::execution::par, idx.begin(), idx.end(), ktb_fx);

  if (tie1 < n_pairs && tie2 < n_pairs) {
    tau = static_cast<double>(nc - nd) / (std::sqrt(n_pairs - tie1) * std::sqrt(n_pairs - tie2));
  }

  const auto pv = compute_kendallb_significance(nc, nd, size, tie1, tie2);
  return {tau, pv};
}

template <typename Iterator,
          typename = typename std::enable_if<
              std::is_base_of<std::random_access_iterator_tag,
                              typename std::iterator_traits<Iterator>::iterator_category>::value>>
std::pair<double, double> compute_kendall(Iterator v1_b, Iterator v1_e, Iterator v2_b,
                                          Iterator v2_e) {
  return compute_kendall_b_n2(v1_b, v1_e, v2_b, v2_e);
}

template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
std::pair<double, double> compute_kendall(std::vector<N> &v1, std::vector<N> &v2) {
  return compute_kendall_b_n2(v1.begin(), v1.end(), v2.begin(), v2.end());
}
}  // namespace modle