// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span

#include <utility>  // for pair
#include <vector>   // for vector

#include "modle/common/common.hpp"  // for std::uint_fast8_t

namespace modle::correlation {

template <class FP>
class Spearman;

template <class FP = double>
class Pearson {
  friend Spearman<FP>;
  static_assert(std::is_floating_point_v<FP>,
                "Template argument FP should be a floating-point type.");

 public:
  inline Pearson() = default;

  struct Result {
    FP pcc{0};
    FP pvalue{1};
  };

  template <class It1, class It2>
  [[nodiscard]] inline Result operator()(It1 begin1, It1 end1, It2 begin2) const;
  template <class Range1, class Range2>
  [[nodiscard]] inline Result operator()(const Range1& r1, const Range2& r2) const;

 private:
  template <class It1, class It2>
  [[nodiscard]] static inline FP compute_pcc(It1 begin1, It1 end1, It2 begin2);

  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  [[nodiscard]] static inline FP compute_significance(FP pcc, I n);

#ifdef ENABLE_TESTING
 public:
  template <class It1, class It2>
  [[nodiscard]] static inline FP test_compute_pcc(It1 begin1, It1 end1, It2 begin2) {
    return compute_pcc(begin1, end1, begin2);
  };

  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  [[nodiscard]] static inline FP test_compute_significance(FP pcc, I n) {
    return compute_significance(pcc, n);
  }
#endif
};

template <class FP = double>
class Spearman {
  static_assert(std::is_floating_point_v<FP>,
                "Template argument FP should be a floating-point type.");

  std::vector<FP> _rank_buff1{};
  std::vector<FP> _rank_buff2{};
  std::vector<usize> _idx_buff{};

 public:
  inline Spearman() = default;
  inline explicit Spearman(usize n);

  struct Result {
    FP rho{0};
    FP pvalue{1};
  };

  template <class It1, class It2>
  [[nodiscard]] inline Result operator()(It1 begin1, It1 end1, It2 begin2);
  template <class Range1, class Range2>
  [[nodiscard]] inline Result operator()(const Range1& r1, const Range2& r2);

 private:
  // Return the indices corresponding to the sorted vector
  template <class It>
  static inline void sort_by_index(It it_begin, It it_end, std::vector<usize>& idx_buff);

  // This returns a vector of floating point numbers corresponding to the rank of the values from
  // the collection delimited by \p it_begin and \p it_end
  template <class It>
  static inline void compute_element_ranks(It it_begin, It it_end, std::vector<FP>& rank_buff,
                                           std::vector<usize>& idx_buff);

  template <class It1, class It2>
  [[nodiscard]] inline FP compute_rho(It1 begin1, It1 end1, It2 begin2);

  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  [[nodiscard]] static inline FP compute_significance(FP rho, I n);

#ifdef ENABLE_TESTING
 public:
  template <class It>
  static inline void test_sort_by_index(It it_begin, It it_end, std::vector<usize>& idx_buff) {
    return sort_by_index(it_begin, it_end, idx_buff);
  }

  // This returns a vector of floating point numbers corresponding to the rank of the values from
  // the collection delimited by \p it_begin and \p it_end
  template <class It>
  static inline void test_compute_element_ranks(It it_begin, It it_end, std::vector<FP>& rank_buff,
                                                std::vector<usize>& idx_buff) {
    return compute_element_ranks(it_begin, it_end, rank_buff, idx_buff);
  }

  template <class It1, class It2>
  [[nodiscard]] inline FP test_compute_rho(It1 begin1, It1 end1, It2 begin2) {
    return compute_rho(begin1, end1, begin2);
  }

  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  [[nodiscard]] static inline FP test_compute_significance(FP rho, I n) {
    return compute_significance(rho, n);
  }
#endif
};

template <class FP = double>
class SED {
  static_assert(std::is_floating_point_v<FP>,
                "Template argument FP should be a floating-point type.");

 public:
  inline SED() = default;
  template <class It1, class It2>
  [[nodiscard]] inline FP operator()(It1 begin1, It1 end1, It2 begin2) const;
  template <class Range1, class Range2>
  [[nodiscard]] inline FP operator()(const Range1& r1, const Range2& r2) const;
};

}  // namespace modle::correlation

#include "../../correlation_impl.hpp"  // IWYU pragma: export

// IWYU pragma: private, include "../../correlation_impl.hpp"
