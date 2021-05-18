#include <cuda_runtime.h>

#include <cassert>
#include <cstdint>
#include <cub/cub.cuh>
#include <cuda/runtime_api.hpp>

#include "modle/cu/common.hpp"
#include "modle/cu/simulation.hpp"
#include "modle/cu/simulation_internal.hpp"

namespace modle::cu::dev {

__device__ uint32_t detect_collisions_at_5prime_single_threaded(
    const bp_t* rev_unit_pos, const bp_t* fwd_unit_pos, const bp_t* rev_moves,
    collision_t* rev_collision_mask, uint32_t chrom_start_pos, uint32_t num_extr_units) {
  assert(threadIdx.x < 2);  // NOLINT
  auto nrev_units_at_5prime = 0U;

  const auto first_active_fwd_unit = fwd_unit_pos[0];
  for (auto i = 0U; i < num_extr_units; ++i) {
    assert(rev_unit_pos[i] != Simulation::EXTR_UNIT_IS_IDLE);   // NOLINT
    assert(fwd_unit_pos[i] != Simulation::EXTR_UNIT_IS_IDLE);   // NOLINT
    assert(chrom_start_pos + rev_moves[i] <= rev_unit_pos[i]);  // NOLINT

    if (rev_unit_pos[i] == chrom_start_pos) {  // Unit already at the 5'-end
      assert(rev_moves[i] == 0);               // NOLINT
      ++nrev_units_at_5prime;
      rev_collision_mask[i] = Simulation::EXTR_UNIT_AT_CHROM_BOUNDARY;

    } else if (rev_unit_pos[i] > first_active_fwd_unit) {
      // As the rev_unit is located downstream of a fwd unit, a LEF-LEF collision might take place
      // before the rev unit reaches the 5'-end
      break;

    } else if (rev_unit_pos[i] - rev_moves[i] == chrom_start_pos) {
      // Unit will reach the 5'-end by the end of the current epoch
      rev_collision_mask[i] = Simulation::EXTR_UNIT_AT_CHROM_BOUNDARY;
      ++nrev_units_at_5prime;
      break;
    }
  }
  return nrev_units_at_5prime;
}

__device__ uint32_t detect_collisions_at_3prime_single_threaded(
    const bp_t* rev_unit_pos, const bp_t* fwd_unit_pos, const bp_t* fwd_moves,
    collision_t* fwd_collision_mask, uint32_t chrom_end_pos, uint32_t num_extr_units) {
  assert(threadIdx.x < 2);  // NOLINT
  auto nfwd_units_at_3prime = 0U;

  const auto last_active_rev_unit = rev_unit_pos[num_extr_units - 1];
  for (auto i = 0U; i < num_extr_units; ++i) {
    assert(rev_unit_pos[i] != Simulation::EXTR_UNIT_IS_IDLE);  // NOLINT
    assert(fwd_unit_pos[i] != Simulation::EXTR_UNIT_IS_IDLE);  // NOLINT
    assert(fwd_unit_pos[i] + fwd_moves[i] < chrom_end_pos);    // NOLINT
    if (fwd_unit_pos[i] == chrom_end_pos - 1) {
      assert(fwd_moves[i] == 0);  // NOLINT
      ++nfwd_units_at_3prime;
      fwd_collision_mask[i] = Simulation::EXTR_UNIT_AT_CHROM_BOUNDARY;

    } else if (fwd_unit_pos[i] < last_active_rev_unit) {
      break;

    } else if (fwd_unit_pos[i] + fwd_moves[i] == chrom_end_pos - 1) {
      fwd_collision_mask[i] = Simulation::EXTR_UNIT_AT_CHROM_BOUNDARY;
      ++nfwd_units_at_3prime;
      break;
    }
  }

  return nfwd_units_at_3prime;
}

__device__ void extrude(GlobalStateDev* global_state) {
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  auto* rev_pos = global_state->block_states[bid].rev_unit_pos;
  auto* fwd_pos = global_state->block_states[bid].fwd_unit_pos;

  const auto* rev_moves = global_state->block_states[bid].rev_moves_buff;
  const auto* fwd_moves = global_state->block_states[bid].fwd_moves_buff;

  const auto num_active_lefs = global_state->block_states[bid].num_active_lefs;
  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  for (auto i = i0; i < i1; ++i) {
    assert(rev_pos[i] - global_state->tasks[bid].chrom_start >= rev_moves[i]);  // NOLINT
    assert(fwd_pos[i] + fwd_moves[i] < global_state->tasks[bid].chrom_end);     // NOLINT
    rev_pos[i] -= rev_moves[i];
    fwd_pos[i] += fwd_moves[i];
  }
}
}  // namespace modle::cu::dev
