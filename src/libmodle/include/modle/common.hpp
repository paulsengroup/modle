#pragma once

#include <cstdint>

namespace modle {

using Bp = std::size_t;
using Iter = uint32_t;
using Contacts = uint32_t;

namespace dna {
enum Direction : uint8_t { none = 0, fwd = 1, rev = 2, both = 3 };
}

}  // namespace modle
