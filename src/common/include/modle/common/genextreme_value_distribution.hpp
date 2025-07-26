// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <cmath>
#include <limits>

#include "modle/common/random.hpp"

namespace modle {
template <class RealType = double>
class genextreme_value_distribution {
 public:
  // types
  using result_type = RealType;

  class param_type {
    result_type _mu;
    result_type _sigma;
    result_type _xi;

   public:
    using distribution_type = genextreme_value_distribution;

    explicit param_type(result_type mu = 0, result_type sigma = 1, result_type xi = 0.1)
        : _mu(mu), _sigma(sigma), _xi(xi) {
      assert(sigma > 0);
    }

    [[nodiscard]] inline result_type mu() const { return _mu; }
    [[nodiscard]] inline result_type sigma() const { return _sigma; }
    [[nodiscard]] inline result_type xi() const { return _xi; }

    [[nodiscard]] friend inline bool operator==(const param_type& a, const param_type& b) {
      return a._mu == b._mu && a._sigma == b._sigma && a._xi == b._xi;
    }
    [[nodiscard]] friend inline bool operator!=(const param_type& a, const param_type& b) {
      return !(a == b);
    }
  };

 private:
  param_type _p;

 public:
  // constructors and reset functions
  inline explicit genextreme_value_distribution(result_type mu = 0, result_type sigma = 1,
                                                result_type xi = 0.1)
      : _p(param_type(mu, sigma, xi)) {}
  inline explicit genextreme_value_distribution(const param_type& p) : _p(p) {}
  inline void reset() {}

  // generating functions
  template <class URNG>
  [[nodiscard]] inline result_type operator()(URNG& g) {
    return (*this)(g, _p);
  }
  template <class URNG>
  inline result_type operator()(URNG& g, const param_type& p);

  // property functions
  [[nodiscard]] inline result_type mu() const { return _p.mu(); }
  [[nodiscard]] inline result_type sigma() const { return _p.sigma(); }
  [[nodiscard]] inline result_type xi() const { return _p.xi(); }
  [[nodiscard]] inline param_type param() const { return _p; }
  inline void param(const param_type& p) { _p = p; }

  [[nodiscard]] inline constexpr result_type min() const {
    return (std::numeric_limits<RealType>::min)();
  }
  [[nodiscard]] inline constexpr result_type max() const {
    return (std::numeric_limits<RealType>::max)();
  }

  friend inline bool operator==(const genextreme_value_distribution& a,
                                const genextreme_value_distribution& b) {
    return a._p == b._p;
  }
  friend inline bool operator!=(const genextreme_value_distribution& a,
                                const genextreme_value_distribution& b) {
    return a != b;
  }
};

template <class RealType>
template <class URNG>
[[nodiscard]] inline typename genextreme_value_distribution<RealType>::result_type
genextreme_value_distribution<RealType>::operator()(URNG& g, const param_type& p) {
  if (p.xi() == RealType(0)) {
    return (p.mu() - p.sigma()) *
           std::log(-std::log(
               random::generate_canonical<RealType, std::numeric_limits<RealType>::digits>(g)));
  }

  return p.mu() +
         (p.sigma() *
          (RealType(1) -
           std::pow(
               -std::log(
                   random::generate_canonical<RealType, std::numeric_limits<RealType>::digits>(g)),
               p.xi()))) /
             p.xi();
}
}  // namespace modle
