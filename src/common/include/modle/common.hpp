#pragma once

#include <array>                                    // for array
#include <boost/random/bernoulli_distribution.hpp>  // for nernoulli_distribution
#include <cstddef>                                  // IWYU pragma: keep for size_t
#include <cstdint>                                  // for uint32_t, uint_fast8_t
#include <functional>                               // for function<>::result_type

#include "modle/utils.hpp"  // for ndebug_defined

#ifdef USE_XOSHIRO
#include <Xoshiro-cpp/XoshiroCpp.hpp>  // for XoshiroCpp::Xoshiro256PlusPlus XoshiroCpp::SplitMix64
#else
#include <boost/random/random_number_generator.hpp>  // for mt19937_64, seed_seq
#endif

namespace modle {

#ifdef USE_XOSHIRO
[[nodiscard]] inline XoshiroCpp::Xoshiro256PlusPlus PRNG(uint64_t seed) noexcept(
    utils::ndebug_defined()) {
  auto seeder = XoshiroCpp::SplitMix64(seed);
  return XoshiroCpp::Xoshiro256PlusPlus(seeder.generateSeedSequence<4>());
}
#else
[[nodiscard]] inline std::mt19937_64 PRNG(uint64_t seed) noexcept(utils::ndebug_defined()) {
  auto seeder = std::seed_seq{seed};
  return std::mt19937_64{seeder};
}
#endif

using PRNG_t = decltype(std::function{PRNG})::result_type;

using bp_t = size_t;
using collision_t = uint_fast32_t;
using contacts_t = uint32_t;
using bernoulli_trial = boost::random::bernoulli_distribution<double>;

namespace dna {
enum Direction : uint_fast8_t { none = 0, fwd = 1, rev = 2, both = 3 };
}  // namespace dna

#if __cplusplus == 202002L
#define MODLE_LIKELY [[likely]]
#define MODLE_UNLIKELY [[unlikely]]
#else
#define MODLE_LIKELY
#define MODLE_UNLIKELY
#endif

}  // namespace modle
