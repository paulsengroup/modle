#include <cuda_runtime_api.h>
#include <curand_kernel.h>
#include <thrust/binary_search.h>
#include <thrust/find.h>
#include <thrust/pair.h>

#include <cassert>
#include <cstdint>
#include <cub/cub.cuh>
#include <cub/device/device_scan.cuh>

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
  if (num_active_lefs == 0) {
    return;
  }
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

__device__ void detect_primary_lef_lef_collisions(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, const uint32_t* lef_rev_unit_idx,
    const bp_t* rev_moves, const bp_t* fwd_moves, uint32_t num_active_lefs, const bp_t* barrier_pos,
    uint32_t num_barriers, collision_t* rev_collisions, collision_t* fwd_collisions,
    curandStatePhilox4_32_10_t* rng_states, uint32_t num_rev_units_at_5prime,
    uint32_t num_fwd_units_at_3prime) {
  if (num_active_lefs - num_rev_units_at_5prime == 0 ||
      num_active_lefs - num_fwd_units_at_3prime == 0 || num_barriers == 0) {
    return;
  }

  const auto tid = threadIdx.x;

  const auto num_active_rev_units = num_active_lefs - num_rev_units_at_5prime;
  // Adjust the indices for the rev units such that they point to units whose position is
  // different than that of the next unit
  const auto [ir0, ir1] = [&]() {
    const auto chunk_size = (num_active_rev_units + blockDim.x - 1) / blockDim.x;
    auto ir0 = min(num_rev_units_at_5prime + (tid * chunk_size), num_active_rev_units);
    auto ir1 = min(ir0 + chunk_size, num_active_rev_units);

    while (ir0 > 0 && rev_units_pos[ir0] == rev_units_pos[ir0 - 1]) {
      --ir0;
    }
    while (ir1 > ir0 && rev_units_pos[ir1] == rev_units_pos[ir1 - 1]) {
      --ir1;
    }

    return thrust::make_pair(ir0, ir1);
  }();

  const auto [if0, if1] = [&, ir0 = ir0, ir1 = ir1]() {
    auto if0 = lef_rev_unit_idx[ir0];
    auto if1 = lef_rev_unit_idx[ir1];
    assert(rev_units_pos[ir0] <= fwd_units_pos[if0]);  // NOLINT
    assert(rev_units_pos[ir1] <= fwd_units_pos[if1]);  // NOLINT

    if0 -= if0 == 0 ? 0 : 1;
    while (if0 > 0 && fwd_units_pos[if0] >= rev_units_pos[ir0]) {
      --if0;
    }

    if1 -= if1 == 0 ? 0 : 1;
    while (if1 > if0 && fwd_units_pos[if1] >= rev_units_pos[ir1]) {
      --if1;
    }

    return thrust::make_pair(if0, if1);
  }();

  assert(ir0 <= num_active_lefs);  // NOLINT
  assert(ir1 <= num_active_lefs);  // NOLINT
  assert(if0 <= num_active_lefs);  // NOLINT
  assert(if1 <= num_active_lefs);  // NOLINT
  __syncthreads();
  if ((ir0 == ir1 && ir0 == num_active_rev_units) || (if0 == if1 && if0 == ir0)) {
    return;
  }

  auto local_rng_state = rng_states[tid];

  auto rev_idx = ir0;
  auto fwd_idx = if0;
  while (rev_idx < ir1 && fwd_idx < if1) {
    while (rev_units_pos[rev_idx] <= fwd_units_pos[fwd_idx]) {
      if (++rev_idx == ir1) {
        rng_states[tid] = local_rng_state;
        return;
      }
    }

    while (fwd_units_pos[fwd_idx] < rev_units_pos[rev_idx]) {
      if (++fwd_idx == if1) {
        rng_states[tid] = local_rng_state;
        return;
      }
    }
    const auto correct_fwd_idx = fwd_idx != 0;
    fwd_idx -= correct_fwd_idx ? 1 : 0;

    assert(rev_idx < num_active_lefs);                         // NOLINT
    assert(fwd_idx < num_active_lefs);                         // NOLINT
    assert(rev_units_pos[rev_idx] >= fwd_units_pos[fwd_idx]);  // NOLINT
    if (const auto delta = rev_units_pos[rev_idx] - fwd_units_pos[fwd_idx];
        delta > 0 && delta < rev_moves[rev_idx] + fwd_moves[fwd_idx] &&
        (config.probability_of_extrusion_unit_bypass == 0 ||
         bernoulli_trial(&local_rng_state, 1.0F - config.probability_of_extrusion_unit_bypass))) {
      auto& cause_of_collision_rev = rev_collisions[rev_idx];
      auto& cause_of_collision_fwd = fwd_collisions[fwd_idx];

      auto [collision_pos_rev, collision_pos_fwd] = compute_lef_lef_collision_pos(
          rev_units_pos[rev_idx], fwd_units_pos[fwd_idx], rev_moves[rev_idx], fwd_moves[fwd_idx]);

      if (cause_of_collision_rev == Simulation::NO_COLLISIONS &&
          cause_of_collision_fwd == Simulation::NO_COLLISIONS) {
        cause_of_collision_rev = num_barriers + fwd_idx;
        cause_of_collision_fwd = num_barriers + rev_idx;
      } else if (cause_of_collision_rev != Simulation::NO_COLLISIONS &&
                 cause_of_collision_fwd == Simulation::NO_COLLISIONS) {
        assert(cause_of_collision_rev < num_barriers);             // NOLINT
        assert(collision_pos_rev != 0 && collision_pos_fwd != 0);  // NOLINT
        const auto& barr_idx = cause_of_collision_rev;
        if (collision_pos_fwd > barrier_pos[barr_idx]) {
          cause_of_collision_rev = num_barriers + fwd_idx;
          cause_of_collision_fwd = num_barriers + rev_idx;
        } else {
          cause_of_collision_fwd = num_barriers + rev_idx;
        }
      } else if (cause_of_collision_rev == Simulation::NO_COLLISIONS &&
                 cause_of_collision_fwd != Simulation::NO_COLLISIONS) {
        assert(cause_of_collision_fwd < num_barriers);             // NOLINT
        assert(collision_pos_rev != 0 && collision_pos_fwd != 0);  // NOLINT

        const auto& barr_idx = cause_of_collision_fwd;
        if (collision_pos_rev < barrier_pos[barr_idx]) {
          cause_of_collision_rev = num_barriers + fwd_idx;
          cause_of_collision_fwd = num_barriers + rev_idx;
        } else {
          cause_of_collision_rev = num_barriers + fwd_idx;
        }
      }
    }
    fwd_idx += correct_fwd_idx ? 1 : 0;
  }
  rng_states[tid] = local_rng_state;
}

__device__ thrust::pair<bp_t, bp_t> compute_lef_lef_collision_pos(bp_t rev_unit_pos,
                                                                  bp_t fwd_unit_pos, bp_t rev_move,
                                                                  bp_t fwd_move) {
  const auto& rev_speed = rev_move;
  const auto& fwd_speed = fwd_move;

  const auto relative_speed = static_cast<float>(rev_speed + fwd_speed);
  assert(rev_unit_pos >= fwd_unit_pos);  // NOLINT
  const auto time_to_collision = static_cast<float>(rev_unit_pos - fwd_unit_pos) / relative_speed;

  const auto collision_pos =
      fwd_unit_pos + __float2uint_rn(static_cast<float>(fwd_speed) * time_to_collision);
  assert(collision_pos <= rev_unit_pos);  // NOLINT

  if (collision_pos == fwd_unit_pos) {
    assert(collision_pos >= fwd_unit_pos);      // NOLINT
    assert(collision_pos + 1 <= rev_unit_pos);  // NOLINT

    return thrust::make_pair(collision_pos + 1, collision_pos);
  }

  assert(collision_pos > 0);                  // NOLINT
  assert(collision_pos - 1 >= fwd_unit_pos);  // NOLINT

  return thrust::make_pair(collision_pos, collision_pos - 1);
}

__device__ void correct_moves_for_primary_lef_lef_collisions(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, bp_t* rev_moves, bp_t* fwd_moves,
    uint32_t num_active_lefs, uint32_t num_barriers, const collision_t* rev_collisions,
    const collision_t* fwd_collisions) {
  if (num_active_lefs == 0 || num_barriers == 0) {
    return;
  }
  const auto tid = threadIdx.x;
  const auto lower_bound = num_barriers;
  const auto upper_bound = lower_bound + num_active_lefs;

  auto is_primary_lef_lef_collision = [&](const auto i) constexpr {
    return i >= lower_bound && i < upper_bound;
  };

  auto is_lef_bar_collision = [&](const auto i) constexpr { return i < lower_bound; };

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  auto i0 = tid * chunk_size;
  auto i1 = min(i0 + chunk_size, num_active_lefs);

  for (auto rev_idx = i0; rev_idx < i1; ++rev_idx) {
    if (auto rev_collision = rev_collisions[rev_idx]; is_primary_lef_lef_collision(rev_collision)) {
      rev_collision -= lower_bound;
      const auto& rev_pos = rev_units_pos[rev_idx];
      const auto& fwd_idx = rev_collision;
      if (auto fwd_collision = fwd_collisions[fwd_idx];
          is_primary_lef_lef_collision(fwd_collision)) {
        const auto& fwd_pos = fwd_units_pos[fwd_idx];

        auto& rev_move = rev_moves[rev_idx];
        auto& fwd_move = fwd_moves[fwd_idx];

        const auto [p1, p2] = compute_lef_lef_collision_pos(rev_pos, fwd_pos, rev_move, fwd_move);
        assert(rev_pos >= p1);  // NOLINT
        assert(fwd_pos <= p2);  // NOLINT

        rev_move = rev_pos - p1;
        fwd_move = p2 - fwd_pos;
      } else if (is_lef_bar_collision(fwd_collision)) {
        const auto& fwd_pos = fwd_units_pos[fwd_idx];

        auto& rev_move = rev_moves[rev_idx];
        const auto& fwd_move = fwd_moves[fwd_idx];

        assert(rev_pos >= fwd_pos + fwd_move);  // NOLINT
        rev_move = rev_pos - (fwd_pos + fwd_move);
      }
    }
  }
  for (auto fwd_idx = i0; fwd_idx < i1; ++fwd_idx) {
    if (auto fwd_collision = fwd_collisions[fwd_idx]; is_primary_lef_lef_collision(fwd_collision)) {
      fwd_collision -= lower_bound;
      const auto& rev_idx = fwd_collision;

      if (const auto& rev_collision = rev_collisions[rev_idx];
          is_lef_bar_collision(rev_collision)) {
        const auto& rev_move = rev_moves[rev_idx];
        auto& fwd_move = fwd_moves[fwd_idx];

        assert(rev_units_pos[rev_idx] >= fwd_units_pos[fwd_idx] + rev_move);  // NOLINT
        fwd_move = (rev_units_pos[rev_idx] - rev_move) - fwd_units_pos[fwd_idx];
      }
    }
  }
}

__device__ void detect_secondary_lef_lef_collisions(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, const bp_t* rev_moves,
    const bp_t* fwd_moves, uint32_t num_active_lefs, uint32_t num_barriers,
    collision_t* rev_collisions, collision_t* fwd_collisions,
    curandStatePhilox4_32_10_t* rng_states, uint32_t num_rev_units_at_5prime,
    uint32_t num_fwd_units_at_3prime) {
  if (num_active_lefs - (num_rev_units_at_5prime + num_fwd_units_at_3prime) == 0 ||
      num_barriers == 0) {
    return;
  }

  const auto tid = threadIdx.x;
  const auto num_active_rev_units = num_active_lefs - num_rev_units_at_5prime;
  // Adjust the indices for the rev units such that they point to units whose position is
  // different than that of the next unit
  const auto [ir0, ir1] = [&]() {
    const auto chunk_size = (num_active_rev_units + blockDim.x - 1) / blockDim.x;
    auto ir0 = min(num_rev_units_at_5prime + (tid * chunk_size), num_active_rev_units);
    auto ir1 = min(ir0 + chunk_size, num_active_rev_units);

    while (ir0 > 0 && rev_collisions[ir0] != Simulation::NO_COLLISIONS) {
      --ir0;
    }
    while (ir1 > ir0 && rev_units_pos[ir1] != Simulation::NO_COLLISIONS) {
      --ir1;
    }

    return thrust::make_pair(ir0, ir1);
  }();

  const auto [if0, if1] = [&]() {
    const auto chunk_size = (num_active_rev_units + blockDim.x - 1) / blockDim.x;
    auto if0 = min(num_fwd_units_at_3prime + (tid * chunk_size), num_active_rev_units);
    auto if1 = min(if0 + chunk_size, num_active_rev_units);

    while (if0 > 0 && fwd_collisions[if0] != Simulation::NO_COLLISIONS) {
      --if0;
    }
    while (if1 > if0 && fwd_units_pos[if1] != Simulation::NO_COLLISIONS) {
      --if1;
    }

    return thrust::make_pair(if0, if1);
  }();

  assert(ir0 <= num_active_lefs);  // NOLINT
  assert(ir1 <= num_active_lefs);  // NOLINT

  assert(if0 <= num_active_lefs);  // NOLINT
  assert(if1 <= num_active_lefs);  // NOLINT
  __syncthreads();

  if (ir0 == ir1 && if0 == if1) {
    return;
  }

  auto local_rng_state = rng_states[tid];
  const auto offset = num_barriers + num_active_lefs;

  for (auto i = ir0; i < ir1 - 1; ++i) {
    if (rev_collisions[i] != Simulation::NO_COLLISIONS) {
      const auto& rev_pos1 = rev_units_pos[i];
      const auto& rev_pos2 = rev_units_pos[i + 1];

      const auto& rev_move1 = rev_moves[i];
      const auto& rev_move2 = rev_moves[i + 1];

      if (rev_pos2 - rev_move2 <= rev_pos1 - rev_move1 &&
          (config.probability_of_extrusion_unit_bypass == 0 ||
           bernoulli_trial(&local_rng_state, 1.0F - config.probability_of_extrusion_unit_bypass))) {
        rev_collisions[i + 1] = offset + i;
      }
    }
  }

  for (auto i = if0; i < if1 - 1; ++i) {
    if (fwd_collisions[i] != Simulation::NO_COLLISIONS) {
      const auto& fwd_pos1 = fwd_units_pos[i];
      const auto& fwd_pos2 = fwd_units_pos[i + 1];

      const auto& fwd_move1 = fwd_moves[i];
      const auto& fwd_move2 = fwd_moves[i + 1];

      if (fwd_pos1 + fwd_move1 >= fwd_pos2 + fwd_move2 &&
          (config.probability_of_extrusion_unit_bypass == 0 ||
           bernoulli_trial(&local_rng_state, 1.0F - config.probability_of_extrusion_unit_bypass))) {
        fwd_collisions[i] = offset + i + 1;
      }
    }
  }

  rng_states[tid] = local_rng_state;
}

__device__ void correct_moves_for_secondary_lef_lef_collisions(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, bp_t* rev_moves, bp_t* fwd_moves,
    uint32_t num_active_lefs, uint32_t num_barriers, const collision_t* rev_collisions,
    const collision_t* fwd_collisions) {
  if (num_active_lefs == 0 || num_barriers == 0) {
    return;
  }
  const auto tid = threadIdx.x;
  const auto lower_bound = num_barriers + num_active_lefs;
  const auto upper_bound = lower_bound + num_active_lefs;

  auto is_secondary_lef_lef_collision = [&](const auto i) constexpr {
    return i >= lower_bound && i < upper_bound;
  };

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  auto i0 = tid * chunk_size;
  auto i1 = min(i0 + chunk_size, num_active_lefs);

  for (auto rev_idx2 = i0; rev_idx2 < i1; ++rev_idx2) {
    if (auto rev_collision = rev_collisions[rev_idx2];
        is_secondary_lef_lef_collision(rev_collision)) {
      rev_collision -= lower_bound;

      const auto& rev_idx1 = rev_collision;
      const auto move = rev_units_pos[rev_idx2] - (rev_units_pos[rev_idx1] - rev_moves[rev_idx1]);
      rev_moves[rev_idx2] = move > 0U ? move - 1 : 0U;
    }
  }

  for (auto fwd_idx1 = i1; fwd_idx1 < i0; --fwd_idx1) {
    if (auto fwd_collision = fwd_collisions[fwd_idx1];
        is_secondary_lef_lef_collision(fwd_collision)) {
      fwd_collision -= lower_bound;

      const auto& fwd_idx2 = fwd_collision;
      const auto move = (fwd_units_pos[fwd_idx2] + fwd_moves[fwd_idx2]) - fwd_units_pos[fwd_idx1];
      fwd_moves[fwd_idx1] = move > 0U ? move - 1 : 0U;
    }
  }
}

__device__ void generate_lef_unloader_affinities(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, const dna::Direction* barrier_directions,
    const bp_t* lef_rev_idx, const collision_t* rev_collisions, const collision_t* fwd_collisions,
    uint32_t num_active_lefs, uint32_t num_barriers, float* lef_unloader_affinities) {
  (void)fwd_units_pos;
  if (num_active_lefs == 0) {
    return;
  }
  auto is_lef_bar_collision = [&](const auto i) { return i < num_barriers; };

  const auto tid = threadIdx.x;

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = min(tid * chunk_size, num_active_lefs);
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  auto local_total = 0.0F;
  __shared__ float block_total;
  if (tid == 0) {
    block_total = 0.0F;
  }
  __syncthreads();

  for (auto i = i0; i < i1; ++i) {
    const auto j = lef_rev_idx[i];
    if (!is_lef_bar_collision(rev_collisions[i]) || !is_lef_bar_collision(fwd_collisions[j])) {
      assert(fwd_units_pos[j] != Simulation::EXTR_UNIT_IS_IDLE);  // NOLINT
      lef_unloader_affinities[i] = 1.0F;
    } else if (rev_units_pos[i] != Simulation::EXTR_UNIT_IS_IDLE) {
      assert(fwd_units_pos[j] != Simulation::EXTR_UNIT_IS_IDLE);  // NOLINT
      const auto& barr1_idx = rev_collisions[i];
      const auto& barr2_idx = fwd_collisions[j];

      if (barrier_directions[barr1_idx] == dna::rev && barrier_directions[barr2_idx] == dna::fwd) {
        lef_unloader_affinities[i] = 1.0F / static_cast<float>(config.hard_stall_multiplier);
      } else {
        lef_unloader_affinities[i] = 1.0F;
      }
    } else {
      assert(rev_units_pos[i] == Simulation::EXTR_UNIT_IS_IDLE);  // NOLINT
      assert(fwd_units_pos[j] == Simulation::EXTR_UNIT_IS_IDLE);  // NOLINT
      lef_unloader_affinities[i] = 0.0F;
    }
    local_total += lef_unloader_affinities[i];
  }
  atomicAdd(&block_total, local_total);
  __syncthreads();
  if (block_total == 0.0F) {
    return;
  }
  for (auto i = i0; i < i1; ++i) {
    lef_unloader_affinities[i] /= block_total;
  }
}

__device__ void select_and_release_lefs(bp_t* rev_units_pos, bp_t* fwd_units_pos,
                                        const bp_t* lef_rev_idx, const bp_t* lef_fwd_idx,
                                        uint32_t num_active_lefs,
                                        const float* lef_unloader_affinities,
                                        float* lef_unloader_affinities_prefix_sum,
                                        curandStatePhilox4_32_10_t* rng_states) {
  if (num_active_lefs == 0) {
    return;
  }
  const auto tid = threadIdx.x;

  auto local_rng_state = rng_states[tid];
  __shared__ uint32_t num_lefs_to_release;
  if (tid == 0) {
    const auto avg_num_of_lefs_to_release =
        static_cast<double>((config.rev_extrusion_speed + config.fwd_extrusion_speed) *
                            num_active_lefs) /
        static_cast<double>(config.average_lef_lifetime);

    num_lefs_to_release =
        min(num_active_lefs, curand_poisson(&local_rng_state, avg_num_of_lefs_to_release));

    if (num_lefs_to_release > 0) {
      // TODO avoid allocating a new buffer each time?
      size_t tmp_storage_bytes_required = 0;
      auto status = cub::DeviceScan::InclusiveSum(
          nullptr, tmp_storage_bytes_required, lef_unloader_affinities,
          lef_unloader_affinities_prefix_sum, num_active_lefs - 1);
      assert(status == cudaSuccess);
      auto* d_temp_storage = malloc(tmp_storage_bytes_required);
      assert(d_temp_storage);  // NOLINT
      status = cub::DeviceScan::InclusiveSum(
          d_temp_storage, tmp_storage_bytes_required, lef_unloader_affinities,
          lef_unloader_affinities_prefix_sum, num_active_lefs - 1);
      // lef_unloader_affinities_prefix_sum[num_active_lefs + 1] = 1.0F;
      assert(status == cudaSuccess);
      free(d_temp_storage);
    }
  }
  __syncthreads();
  if (num_lefs_to_release == 0) {
    rng_states[tid] = local_rng_state;
    return;
  }

  const auto chunk_size = (num_lefs_to_release + blockDim.x - 1) / blockDim.x;
  const auto i0 = min(chunk_size * tid, num_lefs_to_release);
  const auto i1 = min(i0 + chunk_size, num_lefs_to_release);

  if (i0 == i1) {
    return;
  }

  for (auto i = i0; i < i1; ++i) {
    const auto rev_idx = static_cast<uint32_t>(
        thrust::upper_bound(
            thrust::seq, lef_unloader_affinities_prefix_sum,
            lef_unloader_affinities_prefix_sum + num_active_lefs - 1,
            curand_uniform(&local_rng_state) - lef_unloader_affinities_prefix_sum[0]) -
        lef_unloader_affinities_prefix_sum);
    const auto fwd_idx = lef_rev_idx[rev_idx];
    assert(lef_fwd_idx[fwd_idx] == rev_idx);  // NOLINT
    (void)lef_fwd_idx;

    rev_units_pos[rev_idx] = Simulation::EXTR_UNIT_IS_IDLE;
    fwd_units_pos[fwd_idx] = Simulation::EXTR_UNIT_IS_IDLE;
    // printf("releasing %d:%d (%d/%d)\n", rev_idx, fwd_idx, i, num_lefs_to_release);
  }
  rng_states[tid] = local_rng_state;
}

}  // namespace modle::cu::dev
