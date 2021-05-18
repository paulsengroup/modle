#include <cuda_runtime_api.h>
#include <curand_kernel.h>
#include <thrust/find.h>
#include <thrust/pair.h>

#include <cassert>
#include <cstdint>

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

__device__ void extrude(bp_t* rev_units_pos, bp_t* fwd_units_pos, const bp_t* rev_moves,
                        const bp_t* fwd_moves, uint32_t num_active_lefs, uint32_t chrom_start_pos,
                        uint32_t chrom_end_pos) {
  const auto tid = threadIdx.x;
  (void)chrom_start_pos;
  (void)chrom_end_pos;

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  for (auto i = i0; i < i1; ++i) {
    assert(rev_units_pos[i] - chrom_start_pos >= rev_moves[i]);  // NOLINT
    assert(fwd_units_pos[i] + fwd_moves[i] < chrom_end_pos);     // NOLINT
    /*if (blockIdx.x == 0) {
      printf("i=%d; rev %d -> %d (%d) || fwd %d -> %d (%d);\n", i, rev_units_pos[i],
             rev_units_pos[i] - rev_moves[i], rev_moves[i], fwd_units_pos[i],
             fwd_units_pos[i] + fwd_moves[i], fwd_moves[i]);
    }*/
    rev_units_pos[i] -= rev_moves[i];
    fwd_units_pos[i] += fwd_moves[i];
  }
}

__device__ void detect_lef_bar_collisions(const bp_t* rev_units_pos, const bp_t* fwd_units_pos,
                                          const bp_t* rev_moves, const bp_t* fwd_moves,
                                          uint32_t num_active_lefs, const bp_t* barrier_pos,
                                          const dna::Direction* barrier_directions,
                                          const CTCF::State* barrier_states, uint32_t num_barriers,
                                          collision_t* rev_collisions, collision_t* fwd_collisions,
                                          curandStatePhilox4_32_10_t* rng_states,
                                          uint32_t num_rev_units_at_5prime,
                                          uint32_t num_fwd_units_at_3prime) {
  if (num_active_lefs == 0 || num_barriers == 0) {
    return;
  }

  const auto tid = threadIdx.x;
  auto local_rng_state = rng_states[tid];

  if (const auto num_active_rev_units = num_active_lefs - num_rev_units_at_5prime;
      num_active_rev_units > 0) {
    // Adjust the indices for the rev units such that they point to units whose position is
    // different than that of the next unit
    const auto [jr0, jr1] = [&]() {
      const auto chunk_size = (num_active_rev_units + blockDim.x - 1) / blockDim.x;
      auto jr0 = min(num_rev_units_at_5prime + (tid * chunk_size), num_active_rev_units);
      auto jr1 = min(jr0 + chunk_size, num_active_rev_units);

      while (jr0 > 0 && rev_units_pos[jr0] == rev_units_pos[jr0 - 1]) {
        --jr0;
      }
      while (jr1 > jr0 && rev_units_pos[jr1] == rev_units_pos[jr1 - 1]) {
        --jr1;
      }
      return thrust::make_pair(jr0, jr1);
    }();

    const auto [i0, i1] = [&, thresh = rev_units_pos[jr0] - rev_moves[jr0]]() {
      const auto chunk_size = (num_barriers + blockDim.x - 1) / blockDim.x;
      auto i0 = min(tid * chunk_size, num_barriers);
      auto i1 = min(i0 + chunk_size, num_barriers);

      while (i0 > 0 && barrier_pos[i0] >= thresh) {
        --i0;
      }

      return thrust::make_pair(i0, i1);
    }();

    __syncthreads();
    if ((jr0 == jr1 && jr0 == num_active_rev_units) || (i0 == i1 && i0 == num_barriers)) {
      goto process_fwd_units;
    }

    for (auto i = i0, j = jr0; i < i1 && j < jr1; ++i) {
      if (barrier_states[i] == CTCF::NOT_OCCUPIED) {
        continue;
      }
      const auto pblock = barrier_directions[i] == dna::rev ? config.lef_hard_collision_pblock
                                                            : config.lef_soft_collision_pblock;
      while (rev_units_pos[j] <= barrier_pos[i]) {
        if (++j >= jr1) {
          goto process_fwd_units;
        }
      }
      if (rev_units_pos[j] != Simulation::EXTR_UNIT_IS_IDLE) {  // Is this necessary?
        const auto delta = rev_units_pos[j] - barrier_pos[i];

        if (delta > 0 && delta <= rev_moves[j] &&
            (pblock == 1.0F || bernoulli_trial(&local_rng_state, pblock))) {
          rev_collisions[j] = i;
        }
      }
    }
  }

process_fwd_units:
  __syncthreads();
  // Same as above, but for fwd units
  if (const auto num_active_fwd_units = num_active_lefs - num_fwd_units_at_3prime;
      num_active_fwd_units > 0) {
    const auto [jf0, jf1] = [&]() {
      const auto chunk_size = (num_active_fwd_units + blockDim.x - 1) / blockDim.x;
      auto jf0 = min(num_fwd_units_at_3prime + (tid * chunk_size), num_active_fwd_units);
      auto jf1 = min(jf0 + chunk_size, num_active_fwd_units);

      while (jf0 > 0 && fwd_units_pos[jf0] == fwd_units_pos[jf0 - 1]) {
        --jf0;
      }
      while (jf1 > jf0 && fwd_units_pos[jf1] == fwd_units_pos[jf1 - 1]) {
        --jf1;
      }

      return thrust::make_pair(jf0, jf1);
    }();

    const auto [i0, i1] = [&, thresh = fwd_units_pos[jf0] + fwd_moves[jf0]]() {
      const auto chunk_size = (num_barriers + blockDim.x - 1) / blockDim.x;
      auto i0 = tid * chunk_size;
      auto i1 = min(i0 + chunk_size, num_barriers);

      while (i0 < i1 && barrier_pos[i0] <= thresh) {
        ++i0;
      }

      return thrust::make_pair(i0, i1);
    }();

    __syncthreads();
    if ((jf0 == jf1 && jf0 == num_active_fwd_units) || (i0 == i1 && i0 == num_barriers)) {
      goto end_of_loop;
    }

    for (auto i = i0, j = jf0; i < i1 && j < jf1; ++i) {
      if (barrier_states[i] == CTCF::NOT_OCCUPIED) {
        continue;
      }

      const auto pblock = barrier_directions[i] == dna::fwd ? config.lef_hard_collision_pblock
                                                            : config.lef_soft_collision_pblock;
      while (fwd_units_pos[j] < barrier_pos[i]) {
        if (++j >= jf1) {
          goto end_of_loop;
        }
      }
      j = j > 0 ? j - 1 : 0;

      if (fwd_units_pos[j] != Simulation::EXTR_UNIT_IS_IDLE) {  // Is this necessary?
        const auto delta = barrier_pos[i] - fwd_units_pos[j];

        if (delta > 0 && delta <= fwd_moves[j] &&
            (pblock == 1.0F || bernoulli_trial(&local_rng_state, pblock))) {
          fwd_collisions[j] = i;
        }
      }
    }
  }
end_of_loop:
  __syncthreads();
  rng_states[tid] = local_rng_state;
}

__device__ bool bernoulli_trial(curandStatePhilox4_32_10_t* const state, float prob_of_success) {
  assert(prob_of_success >= 0 && prob_of_success <= 1);  // NOLINT
  return curand_uniform(state) < prob_of_success;
}

__device__ void correct_moves_for_lef_bar_collisions(const bp_t* rev_units_pos,
                                                     const bp_t* fwd_units_pos, bp_t* rev_moves,
                                                     bp_t* fwd_moves, uint32_t num_active_lefs,
                                                     const bp_t* barr_pos, uint32_t num_barriers,
                                                     const collision_t* rev_collisions,
                                                     const collision_t* fwd_collisions) {
  if (num_active_lefs == 0 || num_barriers == 0) {
    return;
  }
  const auto tid = threadIdx.x;
  const auto& upper_bound = num_barriers;

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  auto i0 = tid * chunk_size;
  auto i1 = min(i0 + chunk_size, num_active_lefs);

  for (auto i = i0; i < i1; ++i) {
    if (rev_collisions[i] < upper_bound) {
      const auto barr_idx = rev_collisions[i];
      assert(rev_units_pos[i] > barr_pos[barr_idx]);  // NOLINT
      const auto distance = rev_units_pos[i] - barr_pos[barr_idx];
      rev_moves[i] = distance > 1U ? distance - 1 : 0U;
    }

    if (fwd_collisions[i] < upper_bound) {
      const auto barr_idx = fwd_collisions[i];
      assert(fwd_units_pos[i] < barr_pos[barr_idx]);  // NOLINT
      const auto distance = barr_pos[barr_idx] - fwd_units_pos[i];
      fwd_moves[i] = distance > 1U ? distance - 1 : 0U;
    }
  }
}
}  // namespace modle::cu::dev
