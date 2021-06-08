#pragma once

#include <boost/random/generate_canonical.hpp>
#include <cassert>
#include <cmath>

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
    typedef genextreme_value_distribution distribution_type;

    explicit param_type(result_type mu = 0, result_type sigma = 1, result_type xi = 0.1)  // NOLINT
        : _mu(mu), _sigma(sigma), _xi(xi) {
      assert(mu > 0);     // NOLINT
      assert(sigma > 0);  // NOLINT
    }

    [[nodiscard]] result_type mu() const { return _mu; }
    [[nodiscard]] result_type sigma() const { return _sigma; }
    [[nodiscard]] result_type xi() const { return _xi; }

    [[nodiscard]] friend bool operator==(const param_type& a, const param_type& b) {
      return a._mu == b._mu && a._sigma == b._sigma && a._xi == b._xi;
    }
    [[nodiscard]] friend bool operator!=(const param_type& a, const param_type& b) {
      return !(a == b);
    }
  };

 private:
  param_type _p;

 public:
  // constructors and reset functions
  explicit genextreme_value_distribution(result_type mu = 0, result_type sigma = 1,
                                         result_type xi = 0.1)  // NOLINT
      : _p(param_type(mu, sigma, xi)) {}
  explicit genextreme_value_distribution(const param_type& p) : _p(p) {}
  void reset() {}

  // generating functions
  template <class URNG>
  [[nodiscard]] result_type operator()(URNG& g) {
    return (*this)(g, this->_p);
  }
  template <class URNG>
  result_type operator()(URNG& g, const param_type& p);

  // property functions
  [[nodiscard]] result_type mu() const { return this->_p.mu(); }
  [[nodiscard]] result_type sigma() const { return this->_p.sigma(); }
  [[nodiscard]] result_type xi() const { return this->_p.xi(); }
  [[nodiscard]] param_type param() const { return this->_p; }
  [[nodiscard]] void param(const param_type& p) { this->_p = p; }

  [[nodiscard]] constexpr result_type min() const { return std::numeric_limits<RealType>::min(); }
  [[nodiscard]] constexpr result_type max() const { return std::numeric_limits<RealType>::max(); }

  friend bool operator==(const genextreme_value_distribution& a,
                         const genextreme_value_distribution& b) {
    return a._p == b._p;
  }
  friend bool operator!=(const genextreme_value_distribution& a,
                         const genextreme_value_distribution& b) {
    return a != b;
  }
};

template <class RealType>
template <class URNG>
[[nodiscard]] inline typename genextreme_value_distribution<RealType>::result_type
genextreme_value_distribution<RealType>::operator()(URNG& g, const param_type& p) {
  if (this->xi() == RealType(0)) {
    return (this->mu() - this->sigma()) *
           std::log(-std::log(
               boost::random::generate_canonical<RealType, std::numeric_limits<RealType>::digits>(
                   g)));
  }

  return this->mu() +
         (this->sigma() *
          (RealType(1) - std::pow(-std::log(boost::random::generate_canonical<
                                            RealType, std::numeric_limits<RealType>::digits>(g)),
                                  this->xi()))) /
             this->xi();
}
}  // namespace modle
