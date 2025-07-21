// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <xoshiro-cpp/XoshiroCpp.hpp>

#include "modle/common/common.hpp"

#ifdef MODLE_WITH_BOOST_RANDOM
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/generate_canonical.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#else
#include <random>
#endif

namespace modle::random {
[[nodiscard]] constexpr XoshiroCpp::Xoshiro256PlusPlus PRNG(u64 seed) noexcept(
    utils::ndebug_defined()) {
  auto seeder = XoshiroCpp::SplitMix64(seed);
  return XoshiroCpp::Xoshiro256PlusPlus(seeder.generateSeedSequence<4>());
}

using PRNG_t = decltype(PRNG(u64(1)));

#ifdef MODLE_WITH_BOOST_RANDOM
using bernoulli_trial = boost::random::bernoulli_distribution<double>;
template <class N>
using binomial_distribution = boost::random::binomial_distribution<N>;
template <class N>
using discrete_distribution = boost::random::discrete_distribution<N>;
template <class N, usize bits>
inline N generate_canonical(random::PRNG_t& rand_eng) {
  return boost::random::generate_canonical<N, bits>(rand_eng);
}
template <class N>
using normal_distribution = boost::random::normal_distribution<N>;
template <class N>
using poisson_distribution = boost::random::poisson_distribution<N>;
template <class N>
using uniform_int_distribution = boost::random::uniform_int_distribution<N>;
template <class N>
using uniform_real_distribution = boost::random::uniform_real_distribution<N>;

using random_device = boost::random_device;

#else
using bernoulli_trial = std::bernoulli_distribution;
template <class N>
using binomial_distribution = std::binomial_distribution<N>;
template <class N>
using discrete_distribution = std::discrete_distribution<N>;
template <class N, usize bits>
inline N generate_canonical(random::PRNG_t& rand_eng) {
  return std::generate_canonical<N, bits>(rand_eng);
}
template <class N>
using normal_distribution = std::normal_distribution<N>;
template <class N>
using poisson_distribution = std::poisson_distribution<N>;
template <class N>
using uniform_int_distribution = std::uniform_int_distribution<N>;
template <class N>
using uniform_real_distribution = std::uniform_real_distribution<N>;

using random_device = std::random_device;

#endif

}  // namespace modle::random
