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
#include "modle/cu/simulation_internal.cuh"

namespace modle::cu::dev {

__device__ void generate_initial_loading_epochs(uint32_t* epoch_buff,
                                                curandStatePhilox4_32_10_t* rng_states,
                                                uint32_t num_lefs) {
  if (num_lefs == 0) {
    return;
  }

  const auto tid = threadIdx.x;

  const auto scaling_factor =
      static_cast<float>(4 * config.average_lef_lifetime) / static_cast<float>(config.bin_size);

  const auto chunk_size = (num_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = min(tid * chunk_size, num_lefs);
  const auto i1 = min(i0 + chunk_size, num_lefs);

  if (i0 != i1) {
    auto rng_state_local = rng_states[tid];
    thrust::generate(thrust::seq, epoch_buff + i0, epoch_buff + i1, [&]() {
      return static_cast<uint32_t>(roundf(curand_uniform(&rng_state_local) * scaling_factor));
    });

    rng_states[tid] = rng_state_local;
  }
}

__device__ void select_and_bind_lefs(
    bp_t* rev_units_pos, bp_t* fwd_units_pos, const uint32_t* lef_rev_unit_idx,
    uint32_t* lef_fwd_unit_idx, const uint32_t* loading_epochs,
    const uint32_t* lefs_to_load_per_epoch, uint32_t& epoch_idx, uint32_t num_unique_loading_epochs,
    uint32_t& num_active_lefs, uint32_t current_epoch, uint32_t tot_num_lefs, uint32_t chrom_start,
    uint32_t chrom_end, curandStatePhilox4_32_10_t* rng_states, bool& burnin_completed) {
  const auto tid = threadIdx.x;
  const auto chrom_simulated_size = static_cast<float>(chrom_end - chrom_start - 1);
  auto local_rng_state = rng_states[tid];

  // Determine if an how many new LEFs should be bound in the current iteration
  if (!burnin_completed && tid == 0) {
    assert(epoch_idx < num_unique_loading_epochs);  // NOLINT
    if (loading_epochs[epoch_idx] != current_epoch) {
      assert(current_epoch < loading_epochs[epoch_idx]);  // NOLINT
    } else {
      num_active_lefs += lefs_to_load_per_epoch[epoch_idx++];
      if (num_active_lefs == tot_num_lefs) {
        burnin_completed = true;
      }
    }
  }
  __syncthreads();

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  // Bind idle LEFs. LEFs that are supposed to bind for the first time during the current iteration,
  // are always idling
  for (auto i = i0; i < i1; ++i) {
    if (rev_units_pos[i] == Simulation::EXTR_UNIT_IS_IDLE) {
      const auto j = lef_rev_unit_idx[i];
      assert(i < num_active_lefs);  // NOLINT
      assert(j < num_active_lefs);  // NOLINT
      rev_units_pos[i] =
          static_cast<uint32_t>(roundf(curand_uniform(&local_rng_state) * chrom_simulated_size));

      assert(fwd_units_pos[j] == Simulation::EXTR_UNIT_IS_IDLE);  // NOLINT
      fwd_units_pos[j] = rev_units_pos[i];
      lef_fwd_unit_idx[j] = i;
    }
  }

  __syncthreads();
  rng_states[tid] = local_rng_state;
  __syncthreads();
}

__device__ void update_extr_unit_mappings(uint32_t* lef_rev_unit_idx, uint32_t* lef_fwd_unit_idx,
                                          uint32_t num_active_lefs, dna::Direction direction) {
  assert(direction == dna::fwd || direction == dna::rev);  // NOLINT

  const auto tid = threadIdx.x;

  const auto* source_idx = direction == dna::Direction::rev ? lef_rev_unit_idx : lef_fwd_unit_idx;
  auto* dest_idx = direction == dna::Direction::rev ? lef_fwd_unit_idx : lef_rev_unit_idx;

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  for (auto i = i0; i < i1; ++i) {
    const auto j = source_idx[i];
    assert(i < num_active_lefs);  // NOLINT
    assert(j < num_active_lefs);  // NOLINT
    dest_idx[j] = i;
  }
  __syncthreads();
}

__device__ void reset_collision_masks(collision_t* rev_collisions, collision_t* fwd_collisions,
                                      uint32_t num_active_lefs) {
  if (num_active_lefs == 0) {
    return;
  }
  const auto tid = threadIdx.x;

  constexpr auto fill_value = Simulation::NO_COLLISIONS;

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  if (i0 < i1) {
    thrust::fill(thrust::seq, rev_collisions + i0, rev_collisions + i1, fill_value);
    thrust::fill(thrust::seq, fwd_collisions + i0, fwd_collisions + i1, fill_value);
  }
  __syncthreads();
}

__device__ uint32_t detect_collisions_at_5prime_single_threaded(
    const bp_t* rev_unit_pos, const bp_t* fwd_unit_pos, const bp_t* rev_moves,
    collision_t* rev_collision_mask, uint32_t chrom_start_pos, uint32_t num_extr_units) {
  assert(threadIdx.x < 2);  // NOLINT
  if (num_extr_units == 0) {
    return 0;
  }
  auto nrev_units_at_5prime = 0U;

  const auto first_active_fwd_unit = fwd_unit_pos[0];
  for (auto i = 0U; i < num_extr_units; ++i) {
    assert(rev_unit_pos[i] != Simulation::EXTR_UNIT_IS_IDLE);   // NOLINT
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
  if (num_extr_units == 0) {
    return 0;
  }
  auto nfwd_units_at_3prime = 0U;

  const auto last_active_rev_unit = rev_unit_pos[num_extr_units - 1];
  for (auto i = num_extr_units - 1; i > 0; --i) {
    if (fwd_unit_pos[i] == Simulation::EXTR_UNIT_IS_IDLE) {
      printf("%d -> i=%u; pos=%u;\n", blockIdx.x, i, fwd_unit_pos[i]);
    }
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

    const auto [i0, i1] = [&, partition_pos = rev_units_pos[jr0] - rev_moves[jr0]]() {
      const auto i0 = static_cast<uint32_t>(
          thrust::lower_bound(thrust::seq, barrier_pos + num_rev_units_at_5prime,
                              barrier_pos + num_barriers, partition_pos) -
          barrier_pos);
      const auto chunk_size = (num_barriers + blockDim.x - 1) / blockDim.x;
      auto i1 = min(i0 + chunk_size, num_barriers);

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
      auto jf0 = min(tid * chunk_size, num_active_fwd_units);
      auto jf1 = min(jf0 + chunk_size, num_active_fwd_units);

      while (jf0 > 0 && fwd_units_pos[jf0] == fwd_units_pos[jf0 - 1]) {
        --jf0;
      }
      while (jf1 > jf0 && fwd_units_pos[jf1] == fwd_units_pos[jf1 - 1]) {
        --jf1;
      }

      return thrust::make_pair(jf0, jf1);
    }();

    const auto [i0, i1] = [&, partition_pos = fwd_units_pos[jf0] + fwd_moves[jf0]]() {
      const auto i0 = static_cast<uint32_t>(
          thrust::upper_bound(thrust::seq, barrier_pos, barrier_pos + num_barriers, partition_pos) -
          barrier_pos);

      const auto chunk_size = (num_barriers + blockDim.x - 1) / blockDim.x;
      auto i1 = min(i0 + chunk_size, num_barriers);

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
    if (blockIdx.x == 580 && tid == 0) {
      if (rev_units_pos[ir0] > fwd_units_pos[if0] || rev_units_pos[ir1] > fwd_units_pos[if1]) {
        printf("bid=%d;\nrev_pos=[%u", blockIdx.x, *rev_units_pos);
        for (auto k = 1U; k < num_active_lefs; ++k) {
          printf(", %u", rev_units_pos[k]);
        }
        printf("]\nlef_rev_idx=[%u", *lef_rev_unit_idx);
        for (auto k = 1U; k < num_active_lefs; ++k) {
          printf(", %u", lef_rev_unit_idx[k]);
        }
        printf("]\nfwd_pos=[%u", *fwd_units_pos);
        for (auto k = 1U; k < num_active_lefs; ++k) {
          printf(", %u", fwd_units_pos[k]);
        }
        printf("]\n");
      }
    }
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

__device__ void process_secondary_lef_lef_collisions(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, bp_t* rev_moves, bp_t* fwd_moves,
    uint32_t num_active_lefs, uint32_t num_barriers, collision_t* rev_collisions,
    collision_t* fwd_collisions, curandStatePhilox4_32_10_t* rng_states,
    uint32_t num_rev_units_at_5prime, uint32_t num_fwd_units_at_3prime) {
  if (num_active_lefs == 0 || num_barriers == 0) {
    return;
  }

  const auto tid = threadIdx.x;
  const auto num_active_rev_units = num_active_lefs - num_rev_units_at_5prime;
  const auto num_active_fwd_units = num_active_lefs - num_fwd_units_at_3prime;
  // Adjust the indices for the rev units such that they point to units whose position is
  // different than that of the next unit
  const auto [ir0, ir1] = [&]() {
    const auto chunk_size = (num_active_rev_units + blockDim.x - 1) / blockDim.x;
    auto ir0 = min(num_rev_units_at_5prime + (tid * chunk_size), num_active_rev_units);
    auto ir1 = min(ir0 + chunk_size, num_active_rev_units);

    if (ir0 < ir1) {
      while (ir0 > 0 && rev_units_pos[ir0] < rev_units_pos[ir0 + 1] - rev_units_pos[ir0 + 1]) {
        --ir0;
      }
      while (ir1 > ir0 && rev_units_pos[ir1 - 1] < rev_units_pos[ir1] - rev_units_pos[ir1]) {
        --ir1;
      }
    }

    return thrust::make_pair(ir0, ir1);
  }();

  const auto [if0, if1] = [&]() {
    const auto chunk_size = (num_active_rev_units + blockDim.x - 1) / blockDim.x;
    auto if0 = min(num_fwd_units_at_3prime + (tid * chunk_size), num_active_fwd_units);
    auto if1 = min(if0 + chunk_size, num_active_fwd_units);

    if (if0 < if1) {
      while (if0 > 0 && fwd_units_pos[if0] + fwd_moves[if0] < fwd_units_pos[if0 + 1]) {
        --if0;
      }
      while (if1 > if0 && fwd_units_pos[if1 - 1] + fwd_moves[if1 - 1] < fwd_units_pos[if1]) {
        --if1;
      }
    }

    return thrust::make_pair(if0, if1);
  }();

  assert(ir0 <= num_active_lefs);  // NOLINT
  assert(ir1 <= num_active_lefs);  // NOLINT

  assert(if0 <= num_active_lefs);  // NOLINT
  assert(if1 <= num_active_lefs);  // NOLINT
  __syncthreads();

  auto local_rng_state = rng_states[tid];
  const auto offset = num_barriers + num_active_lefs;

  if (ir0 < ir1) {
    for (auto i = ir0; i < ir1 - 1; ++i) {
      if (rev_collisions[i] != Simulation::NO_COLLISIONS &&
          rev_collisions[i + 1] == Simulation::NO_COLLISIONS) {
        const auto& rev_pos1 = rev_units_pos[i];
        const auto& rev_pos2 = rev_units_pos[i + 1];

        const auto& rev_move1 = rev_moves[i];
        auto& rev_move2 = rev_moves[i + 1];
        if (rev_pos2 - rev_move2 <= rev_pos1 - rev_move1 &&
            (config.probability_of_extrusion_unit_bypass == 0 ||
             bernoulli_trial(&local_rng_state,
                             1.0F - config.probability_of_extrusion_unit_bypass))) {
          rev_collisions[i + 1] = offset + i;
          const auto move = rev_pos2 - (rev_pos1 - rev_move1);
          rev_move2 = move > 0U ? move - 1 : 0U;
        }
      }
    }
  }

  if (if0 < if1) {
    for (auto i = if1 - 1; i > if0 + 1; --i) {
      if (fwd_collisions[i - 1] == Simulation::NO_COLLISIONS &&
          fwd_collisions[i] != Simulation::NO_COLLISIONS) {
        const auto& fwd_pos1 = fwd_units_pos[i - 1];
        const auto& fwd_pos2 = fwd_units_pos[i];

        auto& fwd_move1 = fwd_moves[i - 1];
        const auto& fwd_move2 = fwd_moves[i];

        if (fwd_pos1 + fwd_move1 >= fwd_pos2 + fwd_move2 &&
            (config.probability_of_extrusion_unit_bypass == 0 ||
             bernoulli_trial(&local_rng_state,
                             1.0F - config.probability_of_extrusion_unit_bypass))) {
          fwd_collisions[i - 1] = offset + i;
          const auto move = (fwd_pos2 + fwd_move2) - fwd_pos1;
          fwd_move1 = move > 0U ? move - 1 : 0U;
        }
      }
    }
  }

  rng_states[tid] = local_rng_state;
}

__device__ void generate_lef_unloader_affinities(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, const dna::Direction* barrier_directions,
    const uint32_t* lef_rev_idx, const collision_t* rev_collisions,
    const collision_t* fwd_collisions, uint32_t num_active_lefs, uint32_t num_barriers,
    float* lef_unloader_affinities) {
  (void)fwd_units_pos;
  if (num_active_lefs == 0) {
    return;
  }
  auto is_lef_bar_collision = [&](const auto i) { return i < num_barriers; };

  const auto tid = threadIdx.x;

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  auto local_total = 0.0F;
  __shared__ float block_total;
  if (tid == 0) {
    atomicExch(&block_total, 0.0F);
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
        lef_unloader_affinities[i] = 1.0F / config.hard_stall_multiplier;
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
  if (i0 < i1) {
    thrust::transform(thrust::seq, lef_unloader_affinities + i0, lef_unloader_affinities + i1,
                      lef_unloader_affinities + i0, [&](const auto n) { return n / block_total; });
  }
  __syncthreads();
}

__device__ void select_and_release_lefs(bp_t* rev_units_pos, bp_t* fwd_units_pos,
                                        const uint32_t* lef_rev_idx, const uint32_t* lef_fwd_idx,
                                        uint32_t num_active_lefs,
                                        const float* lef_unloader_affinities,
                                        float* lef_unloader_affinities_prefix_sum,
                                        void* tmp_storage, size_t tmp_storage_bytes,
                                        curandStatePhilox4_32_10_t* rng_states) {
  (void)lef_fwd_idx;
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
      cub::DeviceScan::InclusiveSum(tmp_storage, tmp_storage_bytes, lef_unloader_affinities,
                                    lef_unloader_affinities_prefix_sum, num_active_lefs - 1);
      cudaDeviceSynchronize();
      const auto status = cudaGetLastError();
      assert(status == cudaSuccess);
      // lef_unloader_affinities_prefix_sum[num_active_lefs + 1] = 1.0F;
    }
  }
  __syncthreads();
  if (num_lefs_to_release == 0) {
    rng_states[tid] = local_rng_state;
    return;
  }

  const auto chunk_size = (num_lefs_to_release + blockDim.x - 1) / blockDim.x;
  const auto i0 = chunk_size * tid;
  const auto i1 = min(i0 + chunk_size, num_lefs_to_release);

  for (auto i = i0; i < i1; ++i) {
    const auto rev_idx = static_cast<uint32_t>(
        thrust::upper_bound(
            thrust::seq, lef_unloader_affinities_prefix_sum,
            lef_unloader_affinities_prefix_sum + num_active_lefs - 1,
            curand_uniform(&local_rng_state) - lef_unloader_affinities_prefix_sum[0]) -
        lef_unloader_affinities_prefix_sum);
    const auto fwd_idx = lef_rev_idx[rev_idx];
    if (lef_fwd_idx[fwd_idx] != rev_idx) {
      printf("bid=%d;\nrev_pos=[%u", blockIdx.x, *rev_units_pos);
      for (auto k = 1U; k < num_active_lefs; ++k) {
        printf(", %u", rev_units_pos[k]);
      }
      printf("]\nlef_rev_idx=[%u", *lef_rev_idx);
      for (auto k = 1U; k < num_active_lefs; ++k) {
        printf(", %u", lef_rev_idx[k]);
      }
      printf("]\nfwd_pos=[%u", *fwd_units_pos);
      for (auto k = 1U; k < num_active_lefs; ++k) {
        printf(", %u", fwd_units_pos[k]);
      }
      printf("]\nlef_fwd_idx=[%u", *lef_fwd_idx);
      for (auto k = 1U; k < num_active_lefs; ++k) {
        printf(", %u", lef_fwd_idx[k]);
      }
      printf("]\n");
    }
    assert(lef_fwd_idx[fwd_idx] == rev_idx);  // NOLINT

    atomicExch(rev_units_pos + rev_idx, Simulation::EXTR_UNIT_IS_IDLE);
    atomicExch(fwd_units_pos + fwd_idx, Simulation::EXTR_UNIT_IS_IDLE);
    // printf("releasing %d:%d (%d/%d)\n", rev_idx, fwd_idx, i, num_lefs_to_release);
  }
  rng_states[tid] = local_rng_state;
}

__device__ void shuffle_lefs(uint32_t* shuffled_lef_idx_buff, uint32_t* keys_buff,
                             uint32_t* tmp_lef_buff1, uint32_t* tmp_lef_buff2, void* tmp_storage,
                             size_t tmp_storage_bytes, uint32_t num_active_lefs,
                             curandStatePhilox4_32_10_t* rng_states) {
  if (num_active_lefs == 0) {
    return;
  }
  const auto tid = threadIdx.x;

  const auto scaling_factor = static_cast<float>(num_active_lefs);
  const auto [i0, i1] = compute_chunk_size(num_active_lefs);

  if (i0 != i1) {
    auto local_rng_state = rng_states[tid];
    thrust::generate(thrust::seq, keys_buff + i0, keys_buff + i1, [&]() {
      return static_cast<uint32_t>(curand_uniform(&local_rng_state) * scaling_factor);
    });
    rng_states[tid] = local_rng_state;
    thrust::sequence(thrust::seq, shuffled_lef_idx_buff + i0, shuffled_lef_idx_buff + i1, i0);
  }
  __syncthreads();

  if (tid == 0) {
    cub::DoubleBuffer keys_tmp_buff(keys_buff, tmp_lef_buff1);
    cub::DoubleBuffer vals_tmp_buff(shuffled_lef_idx_buff, tmp_lef_buff2);
    cub::DeviceRadixSort::SortPairs(tmp_storage, tmp_storage_bytes, keys_tmp_buff, vals_tmp_buff,
                                    num_active_lefs);
    cudaDeviceSynchronize();
    const auto status = cudaGetLastError();
    assert(status == cudaSuccess);  // NOLINT TODO: Do proper error handling
    memcpy(shuffled_lef_idx_buff, vals_tmp_buff.Current(), num_active_lefs * sizeof(uint32_t));
  }
  __syncthreads();
}

__device__ thrust::pair<uint32_t, uint32_t> compute_chunk_size(uint32_t num_elements) {
  const auto tid = threadIdx.x;

  const auto chunk_size = (num_elements + blockDim.x - 1) / blockDim.x;
  const auto i0 = min(chunk_size * tid, num_elements);
  const auto i1 = min(i0 + chunk_size, num_elements);

  return thrust::make_pair(i0, i1);
}

__device__ void register_contacts(const bp_t* rev_unit_pos, const bp_t* fwd_unit_pos,
                                  const uint32_t* shuffled_idx, const uint32_t* lef_rev_unit_idx,
                                  uint2* contacts_buff, const uint32_t chrom_start,
                                  const uint32_t chrom_end, const uint32_t bin_size,
                                  const uint32_t sample_size, const uint32_t num_active_lefs,
                                  uint32_t* contacts_buff_size_shared,
                                  const uint32_t contatcts_buff_capacity) {
  const auto [i0, i1] = modle::cu::dev::compute_chunk_size(sample_size);

  const auto first_bin = chrom_start / config.bin_size;
  const auto last_bin = chrom_end / config.bin_size;

  for (auto i = i0; i < i1; ++i) {
    const auto rev_idx = shuffled_idx[i];
    const auto fwd_idx = lef_rev_unit_idx[rev_idx];
    assert(rev_idx < num_active_lefs);  // NOLINT
    assert(fwd_idx < num_active_lefs);  // NOLINT

    const auto rev_pos = rev_unit_pos[rev_idx] - chrom_start;
    const auto fwd_pos = fwd_unit_pos[fwd_idx] - chrom_start;

    const auto bin1 = static_cast<uint32_t>(roundf(static_cast<float>(rev_pos) / bin_size));
    const auto bin2 = static_cast<uint32_t>(roundf(static_cast<float>(fwd_pos) / bin_size));

    if (bin1 > first_bin && bin2 > first_bin && bin1 < last_bin && bin2 < last_bin) {
      if (const auto j = atomicAdd(contacts_buff_size_shared, 1); j < contatcts_buff_capacity) {
        assert(j < contatcts_buff_capacity);  // NOLINT
        contacts_buff[j].x = bin1;
        contacts_buff[j].y = bin2;
      } else {
        atomicExch(contacts_buff_size_shared, contatcts_buff_capacity);
        break;
      }
    }
  }
  __syncthreads();
}

}  // namespace modle::cu::dev
