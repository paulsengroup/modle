#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <device_launch_parameters.h>

#include <cassert>
#include <cstdint>  // for uint_fast8_t

#include "modle/common.cuh"  // for bp_t, Direction, PRNG_t

namespace modle::cu {

class ExtrusionBarrier {
 public:
  __host__ ExtrusionBarrier() = default;
  __host__ __device__ inline ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                                              double transition_prob_non_blocking_to_non_blocking,
                                              dna::Direction motif_direction);
  __host__ __device__ inline ExtrusionBarrier(bp_t pos, double transition_prob_blocking_to_blocking,
                                              double transition_prob_non_blocking_to_non_blocking,
                                              char motif_direction);

  [[nodiscard]] __host__ __device__ inline bp_t pos() const;
  [[nodiscard]] __device__ inline double prob_occupied_to_occupied() const;
  [[nodiscard]] __device__ inline double prob_occupied_to_not_occupied() const;
  [[nodiscard]] __device__ inline double prob_not_occupied_to_not_occupied() const;
  [[nodiscard]] __device__ inline double prob_not_occupied_to_occupied() const;
  [[nodiscard]] __host__ __device__ inline dna::Direction blocking_direction_major() const;
  [[nodiscard]] __host__ __device__ inline dna::Direction blocking_direction_minor() const;
  [[nodiscard]] __host__ __device__ inline bool operator<(const ExtrusionBarrier& other) const;
  [[nodiscard]] __device__ inline static double
  compute_blocking_to_blocking_transition_probabilities_from_pblock(
      double probability_of_barrier_block, double non_blocking_to_non_blocking_transition_prob);

 protected:
  bp_t _pos;                                             // NOLINT
  double _occupied_to_occupied_transition_prob;          // NOLINT
  double _non_occupied_to_not_occupied_transition_prob;  // NOLINT
  dna::Direction _blocking_direction;                    // NOLINT
};

namespace CTCF {
enum State : uint_fast8_t { NOT_OCCUPIED = 0, OCCUPIED = 1 };
__device__ inline State next_state(CTCF::State current_state, double occupied_self_transition_prob,
                                   double not_occupied_self_transition_prob, double prob);
//! Update CTCF states for the current iteration based on the states from the previous
//! iteration.

//! \param extr_barriers
//! \param mask bitset used to store CTCF states. States will be updated inplace. Bits set to 0
//!        represent CTCFs in NOT_OCCUPIED state, while bits set to 1 represent CTCFs in
//!        OCCUPIED state
//! \param rand_eng
__device__ inline void update_states(const ExtrusionBarrier* extr_barriers, size_t nbarriers,
                                     CTCF::State* mask, curandStatePhilox4_32_10_t* state);

[[nodiscard]] __host__ __device__ inline dna::Direction major_blocking_dir_to_motif_dir(
    dna::Direction d);

}  // namespace CTCF

__host__ __device__ ExtrusionBarrier::ExtrusionBarrier(
    bp_t pos, double transition_prob_blocking_to_blocking,
    double transition_prob_non_blocking_to_non_blocking, dna::Direction motif_direction)
    : _pos(pos),
      _occupied_to_occupied_transition_prob(transition_prob_blocking_to_blocking),
      _non_occupied_to_not_occupied_transition_prob(transition_prob_non_blocking_to_non_blocking),
      _blocking_direction(motif_direction == dna::Direction::fwd ? dna::Direction::rev
                                                                 : dna::Direction::fwd) {}

__host__ __device__ ExtrusionBarrier::ExtrusionBarrier(
    bp_t pos, double transition_prob_blocking_to_blocking,
    double transition_prob_non_blocking_to_non_blocking, char motif_direction)
    : _pos(pos),
      _occupied_to_occupied_transition_prob(transition_prob_blocking_to_blocking),
      _non_occupied_to_not_occupied_transition_prob(transition_prob_non_blocking_to_non_blocking),
      _blocking_direction(motif_direction == '+' ? dna::Direction::rev : dna::Direction::fwd) {
  // NOLINTNEXTLINE
  assert(motif_direction == '+' || motif_direction == '-');
  assert(transition_prob_blocking_to_blocking >= 0.0 &&  // NOLINT
         transition_prob_blocking_to_blocking <= 1.0);
  assert(transition_prob_non_blocking_to_non_blocking >= 0.0 &&  // NOLINT
         transition_prob_non_blocking_to_non_blocking <= 1.0);
}

__host__ __device__ bp_t ExtrusionBarrier::pos() const { return this->_pos; }

__device__ double ExtrusionBarrier::prob_occupied_to_occupied() const {
  return this->_occupied_to_occupied_transition_prob;
}

__device__ double ExtrusionBarrier::prob_occupied_to_not_occupied() const {
  return 1.0 - prob_occupied_to_occupied();
}

__device__ double ExtrusionBarrier::prob_not_occupied_to_not_occupied() const {
  return this->_non_occupied_to_not_occupied_transition_prob;
}

__device__ double ExtrusionBarrier::prob_not_occupied_to_occupied() const {
  return 1.0 - prob_not_occupied_to_not_occupied();
}

__device__ dna::Direction ExtrusionBarrier::blocking_direction_major() const {
  return this->_blocking_direction;
}

__device__ dna::Direction ExtrusionBarrier::blocking_direction_minor() const {
  return this->_blocking_direction == dna::fwd ? dna::rev : dna::fwd;
}

__host__ __device__ bool ExtrusionBarrier::operator<(const ExtrusionBarrier& other) const {
  return this->pos() < other.pos();
}

__device__ double
ExtrusionBarrier::compute_blocking_to_blocking_transition_probabilities_from_pblock(
    double probability_of_barrier_block, double non_blocking_to_non_blocking_transition_prob) {
  // pno = Transition prob. from non-occupied to occupied
  // pon = Transition prob. from occupied to non-occupied
  // occ = Occupancy
  const auto pno = 1.0 - non_blocking_to_non_blocking_transition_prob;
  const auto occ = probability_of_barrier_block;
  const auto pon = (pno - (occ * pno)) / occ;
  auto pblock = 1.0 - pon;

  if (pblock < 0.0) {
    pblock = 0.0;
  } else if (pblock > 1.0) {
    pblock = 1.0;
  }

  return pblock;
}

__device__ CTCF::State CTCF::next_state(CTCF::State current_state,
                                        double occupied_self_transition_prob,
                                        double not_occupied_self_transition_prob, double prob) {
  assert(occupied_self_transition_prob >= 0 && occupied_self_transition_prob <= 1);  // NOLINT
  assert(not_occupied_self_transition_prob >= 0 &&                                   // NOLINT
         not_occupied_self_transition_prob <= 1);
  assert(prob >= 0 && prob <= 1);  // NOLINT

  if (current_state == NOT_OCCUPIED && prob > not_occupied_self_transition_prob) {
    return OCCUPIED;
  }
  if (current_state == OCCUPIED && prob > occupied_self_transition_prob) {
    return NOT_OCCUPIED;
  }

  return current_state;
}

__device__ void CTCF::update_states(const ExtrusionBarrier* const extr_barriers, size_t nbarriers,
                                    CTCF::State* mask, curandStatePhilox4_32_10_t* state) {
  const auto id = threadIdx.x + blockIdx.x * blockDim.x;
  const auto chunk_size = (nbarriers + 1) / blockDim.x;
  // printf("id=%u; chunk_size=%lu;\n", id, chunk_size);

  auto local_state = state[id];

  const auto i0 = id * chunk_size;
  const auto i1 = i0 + chunk_size;
  auto i = i0;

  do {
    const auto buff = curand_uniform2_double(&local_state);
    mask[i] = CTCF::next_state(mask[i] ? CTCF::OCCUPIED : CTCF::NOT_OCCUPIED,
                               extr_barriers[i].prob_occupied_to_occupied(),
                               extr_barriers[i].prob_not_occupied_to_not_occupied(), buff.x);
    if (++i < i1) {
      mask[i] = CTCF::next_state(mask[i] ? CTCF::OCCUPIED : CTCF::NOT_OCCUPIED,
                                 extr_barriers[i].prob_occupied_to_occupied(),
                                 extr_barriers[i].prob_not_occupied_to_not_occupied(), buff.y);
      ++i;
    }
  } while (i < i1 && i < nbarriers);

  state[id] = local_state;
  // __syncthreads();
}

__host__ __device__ inline dna::Direction CTCF::major_blocking_dir_to_motif_dir(dna::Direction d) {
  return d == dna::fwd ? dna::rev : dna::fwd;
}

}  // namespace modle::cu
