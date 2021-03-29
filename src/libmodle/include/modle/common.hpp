#pragma once

#include <cstddef>  // IWYU pragma: keep for size_t
#include <cstdint>  // for uint32_t, uint8_t

namespace modle {

using bp_t = std::size_t;
using Contacts = uint32_t;

namespace dna {
enum Direction : uint8_t { none = 0, fwd = 1, rev = 2, both = 3 };
}

}  // namespace modle
