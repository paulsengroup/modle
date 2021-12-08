// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <array>        // for array
#include <type_traits>  // for result_of

#include "modle/common/common.hpp"  // for u64
#include "modle/common/utils.hpp"   // for ndebug_defined

#if !defined(MODLE_USE_MERSENNE_TWISTER)
#define USE_XOSHIRO
#endif

#ifdef USE_XOSHIRO
#include <Xoshiro-cpp/XoshiroCpp.hpp>  // for SplitMix64, Xoshiro256PlusPlus
#endif

#if defined(MODLE_USE_MERSENNE_TWISTER) && defined(MODLE_WITH_BOOST_RANDOM)
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/seed_seq.hpp>
#endif

#ifdef MODLE_WITH_BOOST_RANDOM
#include <boost/random/bernoulli_distribution.hpp>     // for bernoulli_distribution
#include <boost/random/discrete_distribution.hpp>      // for discrete_distribution
#include <boost/random/generate_canonical.hpp>         // for generate_canonical
#include <boost/random/normal_distribution.hpp>        // for normal_distribution
#include <boost/random/poisson_distribution.hpp>       // for poisson_distribution
#include <boost/random/random_device.hpp>              // for random_device
#include <boost/random/uniform_int_distribution.hpp>   // for uniform_int_distribution
#include <boost/random/uniform_real_distribution.hpp>  // for uniform_real_distribution
#else
#include <random>
#endif

namespace modle::random {

#if defined(MODLE_USE_MERSENNE_TWISTER) && defined(MODLE_WITH_BOOST_RANDOM)
[[nodiscard]] inline boost::random::mt19937_64 PRNG(u64 seed) noexcept(utils::ndebug_defined()) {
  auto seeder = boost::random::seed_seq{seed};
  return boost::random::mt19937_64{seeder};
}
#endif

#if defined(MODLE_USE_MERSENNE_TWISTER) && !defined(MODLE_WITH_BOOST_RANDOM)
[[nodiscard]] inline std::mt19937_64 PRNG(u64 seed) noexcept(utils::ndebug_defined()) {
  auto seeder = std::seed_seq{seed};
  return std::mt19937_64{seeder};
}
#endif

#ifdef USE_XOSHIRO
[[nodiscard]] inline XoshiroCpp::Xoshiro256PlusPlus PRNG(u64 seed) noexcept(
    utils::ndebug_defined()) {
  auto seeder = XoshiroCpp::SplitMix64(seed);
  return XoshiroCpp::Xoshiro256PlusPlus(seeder.generateSeedSequence<4>());
}
#endif

#ifdef __APPLE__
using PRNG_t = std::invoke_result<decltype (&PRNG)(u64), u64>;
#else
using PRNG_t = decltype(std::function{PRNG})::result_type;
#endif

#ifdef MODLE_WITH_BOOST_RANDOM
using bernoulli_trial = boost::random::bernoulli_distribution<double>;
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
