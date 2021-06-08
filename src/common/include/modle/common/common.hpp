#pragma once

#include <array>    // for array
#include <cstddef>  // IWYU pragma: keep for size_t
#include <cstdint>  // for uint32_t, uint_fast8_t

#include "modle/common/utils.hpp"  // for ndebug_defined

namespace modle {

using bp_t = size_t;
using collision_t = uint_fast32_t;
using contacts_t = uint32_t;

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
