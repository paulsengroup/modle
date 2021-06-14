#pragma once

#include <functional>  // for function<>::result_type
#include <limits>      // for numeric_limits

#include "modle/common/utils.hpp"

#if !defined(USE_MERSENNE_TWISTER)
#define USE_XOSHIRO
#endif

#ifdef USE_XOSHIRO
#include <Xoshiro-cpp/XoshiroCpp.hpp>  // for XoshiroCpp::Xoshiro256PlusPlus XoshiroCpp::SplitMix64
#endif

#if defined(USE_MERSENNE_TWISTER) && defined(USE_BOOST_RANDOM_LIB)
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/seed_seq.hpp>
#endif

#ifdef USE_BOOST_RANDOM_LIB
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/generate_canonical.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#else
#include <random>
#endif

namespace modle::random {

#if defined(USE_MERSENNE_TWISTER) && defined(USE_BOOST_RANDOM_LIB)
[[nodiscard]] inline boost::random::mt19937_64 PRNG(uint64_t seed) noexcept(
    utils::ndebug_defined()) {
  auto seeder = boost::random::seed_seq{seed};
  return boost::random::mt19937_64{seeder};
}
#endif

#if defined(USE_MERSENNE_TWISTER) && !defined(USE_BOOST_RANDOM_LIB)
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

using PRNG_t = decltype(std::function{PRNG})::result_type;

#ifdef USE_BOOST_RANDOM_LIB
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
