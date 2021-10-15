// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <array>       // for array
#include <cstddef>     // for size_t
#include <cstdint>     // for uint64_t
#include <functional>  // for function<>::result_type
#include <limits>      // for numeric_limits

#include "modle/common/utils.hpp"  // for ndebug_defined

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
#include <boost/random/uniform_int_distribution.hpp>   // for uniform_int_distribution
#include <boost/random/uniform_real_distribution.hpp>  // for uniform_real_distribution
#else
#include <random>
#endif

namespace modle::random {

#if defined(MODLE_USE_MERSENNE_TWISTER) && defined(MODLE_WITH_BOOST_RANDOM)
[[nodiscard]] inline boost::random::mt19937_64 PRNG(uint64_t seed) noexcept(
    utils::ndebug_defined()) {
  auto seeder = boost::random::seed_seq{seed};
  return boost::random::mt19937_64{seeder};
}
#endif

#if defined(MODLE_USE_MERSENNE_TWISTER) && !defined(MODLE_WITH_BOOST_RANDOM)
[[nodiscard]] inline std::mt19937_64 PRNG(uint64_t seed) noexcept(utils::ndebug_defined()) {
  auto seeder = std::seed_seq{seed};
  return std::mt19937_64{seeder};
}
#endif

#ifdef USE_XOSHIRO
[[nodiscard]] inline XoshiroCpp::Xoshiro256PlusPlus PRNG(uint64_t seed) noexcept(
    utils::ndebug_defined()) {
  auto seeder = XoshiroCpp::SplitMix64(seed);
  return XoshiroCpp::Xoshiro256PlusPlus(seeder.generateSeedSequence<4>());
}
#endif

using PRNG_t = std::result_of<decltype (&PRNG)(uint64_t)>::type;

#ifdef MODLE_WITH_BOOST_RANDOM
using bernoulli_trial = boost::random::bernoulli_distribution<double>;
template <class N>
using discrete_distribution = boost::random::discrete_distribution<N>;
template <class N, size_t bits>
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

#else
using bernoulli_trial = std::bernoulli_distribution;
template <class N>
using discrete_distribution = std::discrete_distribution<N>;
template <class N, size_t bits>
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

#endif

}  // namespace modle::random
