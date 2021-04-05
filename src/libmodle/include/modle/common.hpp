#pragma once

#include <cstddef>  // IWYU pragma: keep for size_t
#include <cstdint>  // for uint32_t, uint_fast8_t

#ifdef USE_XOSHIRO
#include <Xoshiro-cpp/XoshiroCpp.hpp>  // for XoshiroCpp::Xoshiro256PlusPlus XoshiroCpp::SplitMix64
#else
#include <random>  // for mt19937_64, seed_seq
#endif

namespace modle {

#ifdef USE_XOSHIRO
[[nodiscard]] inline XoshiroCpp::Xoshiro256PlusPlus PRNG(uint64_t seed) {
  auto seeder = XoshiroCpp::SplitMix64(seed);
  return XoshiroCpp::Xoshiro256PlusPlus(seeder.generateSeedSequence<4>());
}
#else
[[nodiscard]] inline std::mt19937_64 PRNG(uint64_t seed) {
  auto seeder = std::seed_seq{seed};
  return std::mt19937_64{seeder};
}
#endif

using PRNG_t = decltype(std::function{PRNG})::result_type;

using bp_t = std::size_t;
using Contacts = uint32_t;

namespace dna {
enum Direction : uint_fast8_t { none = 0, fwd = 1, rev = 2, both = 3 };
}

}  // namespace modle
