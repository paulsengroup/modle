#pragma once
#include <cstdint>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace modle {

template <typename N>
class CorrelationTest {
  static_assert(std::is_arithmetic<N>::value,
                "ContactMatrix requires a numeric type as template argument.");

 public:
  CorrelationTest(const std::vector<N>& v1, const std::vector<N>& v2);
  [[nodiscard]] std::pair<double, double> compute_spearman() const;
  [[nodiscard]] std::vector<std::pair<double, double>> compute_spearman(
      uint64_t window_size, uint64_t window_overlap) const;
  [[nodiscard]] std::pair<double, double> compute_kendall() const;
  [[nodiscard]] std::vector<std::pair<double, double>> compute_kendall(
      uint64_t window_size, uint64_t window_overlap) const;

 private:
  const std::vector<N>& _v1;
  const std::vector<N>& _v2;

  [[nodiscard]] static double compute_spearman_significance(double rho, uint32_t n);
  [[nodiscard]] static double compute_kendall_b_significance(uint32_t nc, uint32_t nd,
                                                             uint32_t size, uint32_t tie1,
                                                             uint32_t tie2);
  template <typename Iterator,
            typename = typename std::enable_if<
                std::is_base_of<std::random_access_iterator_tag,
                                typename std::iterator_traits<Iterator>::iterator_category>::value>>
  [[nodiscard]] static std::pair<double, double> compute_kendall_b_n2(Iterator v1_b, Iterator v1_e,
                                                                      Iterator v2_b);
  template <typename Iterator,
            typename = typename std::enable_if<
                std::is_base_of<std::random_access_iterator_tag,
                                typename std::iterator_traits<Iterator>::iterator_category>::value>>
  [[nodiscard]] static std::pair<double, double> compute_kendall(Iterator v1_b, Iterator v1_e,
                                                                 Iterator v2_b);
  template <typename Iterator,
            typename = typename std::enable_if<
                std::is_base_of<std::random_access_iterator_tag,
                                typename std::iterator_traits<Iterator>::iterator_category>::value>>
  [[nodiscard]] static std::pair<double, double> compute_spearman(Iterator v1_b, Iterator v1_e,
                                                                  Iterator v2_b);
};
}  // namespace modle

#include "modle/impl/correlation_test.hpp"
