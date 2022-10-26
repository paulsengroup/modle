// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>  // for enable_if, enable_if_t
#include <vector>       // for vector

#include "modle/common/common.hpp"  // for u8f
#include "modle/common/utils.hpp"   // for RepeatIterator

namespace modle::stats {

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

  template <class It1, class It2, class It3>
  [[nodiscard]] inline Result operator()(It1 first1, It1 last1, It2 first2,
                                         It3 weight_first = utils::RepeatIterator<FP>(1)) const;

  template <class Range1, class Range2>
  [[nodiscard]] inline Result operator()(const Range1& r1, const Range2& r2) const;
  template <class Range1, class Range2, class Range3>
  [[nodiscard]] inline Result operator()(const Range1& r1, const Range2& r2,
                                         const Range3& weights) const;

 private:
  template <class It1, class It2>
  [[nodiscard]] static inline FP compute_pcc(It1 first1, It1 last1, It2 first2);

  template <class It1, class It2, class It3>
  [[nodiscard]] static inline FP compute_weighted_pcc(It1 first1, It1 last1, It2 first2,
                                                      It3 weight_first);

  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  [[nodiscard]] static inline FP compute_significance(FP pcc, I n);

#ifdef MODLE_ENABLE_TESTING
 public:
  template <class It1, class It2>
  [[maybe_unused]] [[nodiscard]] static inline FP test_compute_pcc(It1 first1, It1 last1,
                                                                   It2 first2) {
    return compute_pcc(first1, last1, first2);
  }

  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  [[maybe_unused]] [[nodiscard]] static inline FP test_compute_significance(FP pcc, I n) {
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

  template <class It1, class It2, class It3>
  [[nodiscard]] inline Result operator()(It1 first1, It1 last1, It2 first2,
                                         It3 weight_first = utils::RepeatIterator<FP>(1));

  template <class Range1, class Range2>
  [[nodiscard]] inline Result operator()(const Range1& r1, const Range2& r2);
  template <class Range1, class Range2, class Range3>
  [[nodiscard]] inline Result operator()(const Range1& r1, const Range2& r2, const Range3& weights);

 private:
  template <class It1, class It2>
  [[nodiscard]] inline FP compute_rho(It1 first1, It1 last1, It2 first2);
  template <class It1, class It2, class It3>
  [[nodiscard]] inline FP compute_weighted_rho(It1 first1, It1 last1, It2 first2, It3 weight_first);

  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  [[nodiscard]] static inline FP compute_significance(FP rho, I n);

#ifdef MODLE_ENABLE_TESTING
 public:
  template <class It1, class It2>
  [[maybe_unused]] [[nodiscard]] inline FP test_compute_rho(It1 first1, It1 last1, It2 first2) {
    return compute_rho(first1, last1, first2);
  }

  template <class I, class = std::enable_if<std::is_integral_v<I>>>
  [[maybe_unused]] [[nodiscard]] static inline FP test_compute_significance(FP rho, I n) {
    return compute_significance(rho, n);
  }
#endif
};

namespace internal {

// Return the indices corresponding to the sorted vector
template <class It>
inline void sort_by_index(It first, It last, std::vector<usize>& idx_buff);

// This returns a vector of floating point numbers corresponding to the rank of the values from
// the collection delimited by \p first and \p last
template <class It, class FP, class = std::enable_if_t<std::is_floating_point_v<FP>>>
inline void compute_element_ranks(It first, It last, std::vector<FP>& rank_buff,
                                  std::vector<usize>& idx_buff);

// https://rdrr.io/cran/wCorr/f/inst/doc/wCorrFormulas.pdf
template <class It1, class It2, class FP, class = std::enable_if_t<std::is_floating_point_v<FP>>>
inline void compute_weighted_element_ranks(It1 first, It1 last, It2 weight_first,
                                           std::vector<FP>& rank_buff,
                                           std::vector<usize>& idx_buff);
}  // namespace internal

}  // namespace modle::stats

#include "../../../correlation_impl.hpp"  // IWYU pragma: export

// IWYU pragma: private, include "../../correlation_impl.hpp"
