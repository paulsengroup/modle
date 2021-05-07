#pragma once

#include <fmt/format.h>

#include <cuda/std/cstdint>  // for uint32_t, uint_fast8_t

namespace modle::cu {

#define CUDA_CALL(x)                                                                         \
  do {                                                                                       \
    if ((x) != cudaSuccess) {                                                                \
      throw std::runtime_error(fmt::format("Error at {}:{}: {}: {}\n", __FILE__, __LINE__,   \
                                           cudaGetErrorName((x)), cudaGetErrorString((x)))); \
    }                                                                                        \
  } while (0)

using bp_t = uint32_t;
using collision_t = uint32_t;
using contacts_t = uint32_t;

namespace dna {
enum Direction : uint_fast8_t { none = 0, fwd = 1, rev = 2, both = 3 };
}  // namespace dna

}  // namespace modle::cu
