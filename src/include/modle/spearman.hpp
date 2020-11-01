#pragma once
#include <algorithm>
#include <boost/math/distributions/students_t.hpp>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <vector>

#include "modle/utils.hpp"

namespace modle {

double compute_spearman_significance(double r, uint32_t n) {
  auto dof = n - 2;
  double tscore = r * std::sqrt(dof / ((1.0 + r) * (1.0 - r)));
  boost::math::students_t_distribution<double> dist(dof);
  // See line 4229 at
  // https://github.com/scipy/scipy/blob/6703631bcd15750e86f4098b0421efabcac0f7c2/scipy/stats/stats.py
  // 2 * survival function
  return 2.0 * boost::math::cdf<double>(boost::math::complement(dist, std::fabs<double>(tscore)));
}

template <typename Iterator,
          typename = typename std::enable_if<
              std::is_base_of<std::random_access_iterator_tag,
                              typename std::iterator_traits<Iterator>::iterator_category>::value>>
std::pair<double, double> compute_spearman(Iterator v1_b, Iterator v1_e, Iterator v2_b,
                                           Iterator v2_e) {
  assert(std::distance(v1_b, v1_e) == std::distance(v2_b, v2_e));
  const uint32_t size = std::distance(v1_b, v1_e);
  auto v1r = compute_element_ranks(v1_b, v1_e);
  auto v2r = compute_element_ranks(v2_b, v2_e);

  auto v1r_avg = std::reduce(v1r.begin(), v1r.end(), 0.0) / size;
  auto v2r_avg = std::reduce(v2r.begin(), v2r.end(), 0.0) / size;
  double n = 0;           // numerator: sum (x - xm) * (y - ym)
  double d1 = 0, d2 = 0;  // denominator: sum (x - xm)^2 * sum (y - ym)^2

  for (auto i = 0U; i < size; ++i) {
    const auto &r1 = v1r[i];
    const auto &r2 = v2r[i];
    n += (r1 - v1r_avg) * (r2 - v2r_avg);
    d1 += std::pow(r1 - v1r_avg, 2);
    d2 += std::pow(r2 - v2r_avg, 2);
  }

  auto rho = n / std::sqrt(d1 * d2);
  return {rho, compute_spearman_significance(rho, size)};
}

template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
std::pair<double, double> compute_spearman(std::vector<N> &v1, std::vector<N> &v2) {
  return compute_spearman(v1.begin(), v1.end(), v2.begin(), v2.end());
}

}  // namespace modle