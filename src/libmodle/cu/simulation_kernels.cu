#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>

#include <cstdio>

#include "modle/cu/simulation.hpp"
#include "modle/cu/simulation_internal.hpp"

namespace modle::cu::kernels {

__global__ void init_curand(GlobalStateDev* global_state) {
  const auto seed = global_state->tasks[blockIdx.x].seed;
  curandStatePhilox4_32_10_t rng_state;

  curand_init(seed, threadIdx.x, 0, &rng_state);

  global_state->block_states[blockIdx.x].rng_state[threadIdx.x] = rng_state;
}

__global__ void reset_buffers(GlobalStateDev* global_state) {
  const auto bid = blockIdx.x;
  const auto tid = threadIdx.x;

  const auto nlefs = global_state->tasks[bid].nlefs;

  auto& block_state = global_state->block_states[bid];

  constexpr auto fill_value = Simulation::EXTR_UNIT_IS_IDLE;

  if (threadIdx.x == 0) {
    thrust::fill_n(thrust::device, block_state.rev_unit_pos, nlefs, fill_value);
    thrust::fill_n(thrust::device, block_state.fwd_unit_pos, nlefs, fill_value);
    block_state.num_active_lefs = 0;
    block_state.burnin_completed = false;
    block_state.simulation_completed = false;
    block_state.contact_local_buff_size = 0;
  }
  __syncthreads();

  const auto chunk_size = (nlefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, nlefs);
  for (auto i = i0; i < i1; ++i) {
    block_state.lef_rev_unit_idx[i] = i;
    block_state.lef_fwd_unit_idx[i] = i;
  }
}

__global__ void generate_initial_loading_epochs(GlobalStateDev* global_state, uint32_t num_epochs) {
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;
  const auto id = tid + bid * blockDim.x;
  const auto nthreads = blockDim.x * gridDim.x;

  auto* block_state = global_state->block_states + bid;
  const auto scaling_factor =
      static_cast<float>(4 * config.average_lef_lifetime) / static_cast<float>(config.bin_size);
  auto rng_state_local = block_state->rng_state[tid];

  auto* global_buff = global_state->large_uint_buff1;

  const auto chunk_size = (num_epochs + nthreads - 1) / nthreads;
  const auto i0 = id * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_epochs);

  for (auto i = i0; i < i1; ++i) {
    global_buff[i] = __float2uint_rn(curand_uniform(&rng_state_local) * scaling_factor);
  }

  block_state->rng_state[tid] = rng_state_local;
}

__global__ void init_barrier_states(GlobalStateDev* global_state) {
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  auto local_rng_state = global_state->block_states[bid].rng_state[tid];

  const auto nbarriers = global_state->tasks[bid].nbarriers;
  auto* barrier_states = global_state->block_states[bid].barrier_mask;

  const auto chunk_size = (nbarriers + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, nbarriers);

  for (auto i = i0; i < i1; ++i) {
    const auto pno = 1.0F - global_state->barrier_probs_nocc_to_nocc[i];
    const auto pon = 1.0F - global_state->barrier_probs_occ_to_occ[i];
    const auto occ = pno / (pno + pon);

    barrier_states[i] =
        curand_uniform(&local_rng_state) > occ ? CTCF::NOT_OCCUPIED : CTCF::OCCUPIED;
  }

  global_state->block_states[bid].rng_state[tid] = local_rng_state;
}

__global__ void shift_and_scatter_lef_loading_epochs(const uint32_t* global_buff,
                                                     BlockState* block_states,
                                                     const uint32_t* begin_offsets,
                                                     const uint32_t* end_offsets) {
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  const auto* source_ptr = global_buff + begin_offsets[bid];
  auto* dest_ptr = block_states[bid].epoch_buff;
  const auto size = end_offsets[bid] - begin_offsets[bid];

  if (threadIdx.x == 0) {
    thrust::copy_n(thrust::device, source_ptr, size, dest_ptr);
  }

  __syncthreads();
  const auto offset = *source_ptr;
  const auto chunk_size = (size + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, size);

  for (auto i = i0; i < i1; ++i) {
    dest_ptr[i] -= offset;
  }
}

__global__ void select_and_bind_lefs(uint32_t current_epoch, GlobalStateDev* global_state) {
  if (global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }
  const auto bid = blockIdx.x;
  const auto tid = threadIdx.x;

  const auto chrom_start = global_state->tasks[bid].chrom_start;
  const auto chrom_end = global_state->tasks[bid].chrom_end;
  const auto chrom_simulated_size = __uint2float_rn(chrom_end - chrom_start - 1);

  auto local_rng_state = global_state->block_states[bid].rng_state[tid];

  if (!global_state->block_states[bid].burnin_completed && tid == 0) {
    const auto nlefs = global_state->tasks[bid].nlefs;
    auto num_active_lefs = global_state->block_states[bid].num_active_lefs;
    do {
      if (global_state->block_states[bid].epoch_buff[num_active_lefs] > current_epoch) {
        break;
      }
    } while (++num_active_lefs < nlefs);
    global_state->block_states[bid].num_active_lefs = num_active_lefs;
    if (nlefs == num_active_lefs) {
      global_state->block_states[bid].burnin_completed = true;
    }
  }
  __syncthreads();

  auto* rev_unit_pos = global_state->block_states[bid].rev_unit_pos;
  auto* fwd_unit_pos = global_state->block_states[bid].fwd_unit_pos;

  auto* rev_unit_idx = global_state->block_states[bid].lef_rev_unit_idx;

  const auto num_active_lefs = global_state->block_states[bid].num_active_lefs;

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  for (auto i = i0; i < i1; ++i) {
    if (rev_unit_pos[i] == Simulation::EXTR_UNIT_IS_IDLE) {
      rev_unit_pos[i] =
          __float2uint_rn(roundf(curand_uniform(&local_rng_state) * chrom_simulated_size));
      const auto j = rev_unit_idx[i];

      assert(fwd_unit_pos[j] == Simulation::EXTR_UNIT_IS_IDLE);  // NOLINT
      fwd_unit_pos[j] = rev_unit_pos[i];
    }
  }
  global_state->block_states[bid].rng_state[tid] = local_rng_state;
}

__global__ void prepare_extr_units_for_sorting(GlobalStateDev* global_state,
                                               dna::Direction direction,
                                               uint32_t* tot_num_units_to_sort) {
  assert(direction == dna::Direction::fwd || direction == dna::Direction::rev);  // NOLINT
  if (global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }

  if (threadIdx.x == 0) {
    const auto bid = blockIdx.x;

    const auto buff_alignment = global_state->large_uint_buff_chunk_alignment;
    auto* buff1 = global_state->large_uint_buff1 + (bid * buff_alignment);
    auto* buff2 = global_state->large_uint_buff3 + (bid * buff_alignment);
    auto* start_offsets = global_state->sorting_offset1_buff;
    auto* end_offsets = global_state->sorting_offset2_buff;

    const auto num_active_lefs = global_state->block_states[bid].num_active_lefs;

    if (direction == dna::Direction::rev) {
      thrust::copy_n(thrust::device, global_state->block_states[bid].rev_unit_pos, num_active_lefs,
                     buff1);
      thrust::copy_n(thrust::device, global_state->block_states[bid].lef_rev_unit_idx,
                     num_active_lefs, buff2);
    } else {
      thrust::copy_n(thrust::device, global_state->block_states[bid].fwd_unit_pos, num_active_lefs,
                     buff1);
      thrust::copy_n(thrust::device, global_state->block_states[bid].lef_fwd_unit_idx,
                     num_active_lefs, buff2);
    }

    start_offsets[bid] = buff_alignment * bid;
    end_offsets[bid] = (buff_alignment * bid) + num_active_lefs;
    if (tot_num_units_to_sort) {
      atomicAdd(tot_num_units_to_sort, num_active_lefs);
    }
  }
}

__global__ void update_unit_mappings_and_scatter_sorted_lefs(
    GlobalStateDev* global_state, dna::Direction direction, bool update_extr_unit_to_lef_mappings) {
  if (global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }

  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  auto* block_state = global_state->block_states + bid;

  const auto* source_ptr = global_state->large_uint_buff1 + global_state->sorting_offset1_buff[bid];
  auto* dest_ptr =
      direction == dna::Direction::rev ? block_state->rev_unit_pos : block_state->fwd_unit_pos;
  const auto size = block_state->num_active_lefs;
  assert(size == global_state->sorting_offset2_buff[bid] -
                     global_state->sorting_offset1_buff[bid]);  // NOLINT

  if (tid == 0) {
    thrust::copy_n(thrust::device, source_ptr, size, dest_ptr);
  }
  __syncthreads();

  if (update_extr_unit_to_lef_mappings) {
    const auto* source_idx = direction == dna::Direction::rev ? block_state->lef_rev_unit_idx
                                                              : block_state->lef_fwd_unit_idx;
    auto* dest_idx = direction == dna::Direction::rev ? block_state->lef_fwd_unit_idx
                                                      : block_state->lef_rev_unit_idx;

    const auto chunk_size = (size + blockDim.x - 1) / blockDim.x;
    const auto i0 = tid * chunk_size;
    const auto i1 = min(i0 + chunk_size, size);

    for (auto i = i0; i < i1; ++i) {
      const auto j = source_idx[i];
      dest_idx[j] = i;
    }
  }
}

__global__ void prepare_units_for_random_shuffling(GlobalStateDev* global_state,
                                                   uint32_t* tot_num_active_units) {
  if (global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }

  const auto bid = blockIdx.x;
  const auto buff_alignment = global_state->large_uint_buff_chunk_alignment;
  auto* start_offsets = global_state->sorting_offset1_buff;
  auto* end_offsets = global_state->sorting_offset2_buff;

  const auto block_offset = bid * buff_alignment;

  if (!global_state->block_states[bid].burnin_completed) {
    if (threadIdx.x == 0) {
      start_offsets[bid] = block_offset;
      end_offsets[bid] = block_offset;
    }
    return;
  }

  auto* idx_buff = global_state->large_uint_buff1 + block_offset;

  const auto num_active_lefs = global_state->block_states[bid].num_active_lefs;
  assert(num_active_lefs <= buff_alignment);  // NOLINT
  if (threadIdx.x == 0) {
    thrust::sequence(thrust::device, idx_buff, idx_buff + num_active_lefs);
    start_offsets[bid] = block_offset;
    end_offsets[bid] = block_offset + num_active_lefs;

    if (tot_num_active_units) {
      atomicAdd(tot_num_active_units, num_active_lefs);
    }
  }
  __syncthreads();

  const auto tid = threadIdx.x;
  auto local_rng_state = global_state->block_states[bid].rng_state[tid];
  constexpr auto scaling_factor = static_cast<float>(static_cast<uint32_t>(-1));
  auto* random_num_buff = global_state->large_uint_buff3 + block_offset;

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  for (auto i = i0; i < i1; ++i) {
    random_num_buff[i] = __float2uint_rn(curand_uniform(&local_rng_state) * scaling_factor);
  }

  global_state->block_states[bid].rng_state[tid] = local_rng_state;
}

__global__ void select_lefs_then_register_contacts(GlobalStateDev* global_state) {
  if (global_state->block_states[blockIdx.x].simulation_completed ||
      !global_state->block_states[blockIdx.x].burnin_completed) {
    return;
  }

  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  const auto* rev_unit_idx_buff = global_state->large_uint_buff1;
  const auto* lef_rev_unit_idx_buff = global_state->block_states[bid].lef_rev_unit_idx;

  const auto* rev_pos_buff = global_state->block_states[bid].rev_unit_pos;
  const auto* fwd_pos_buff = global_state->block_states[bid].fwd_unit_pos;

  const auto buff_alignment = global_state->large_uint_buff_chunk_alignment;

  const auto mean_sample_size =
      static_cast<double>(config.lef_fraction_contact_sampling) *
      static_cast<double>(global_state->block_states[bid].num_active_lefs);

  auto local_rng_state = global_state->block_states[bid].rng_state[0];
  __syncthreads();
  const auto sample_size = min(min(curand_poisson(&local_rng_state, mean_sample_size),
                                   global_state->block_states[bid].contact_local_buff_capacity -
                                       global_state->block_states[bid].contact_local_buff_size),
                               global_state->block_states[bid].num_active_lefs);
  if (tid == 0) {
    global_state->block_states[bid].rng_state[0] = local_rng_state;
  }
  __syncthreads();
  if (sample_size == 0) {
    return;
  }
  const auto bin_size = static_cast<float>(config.bin_size);

  const auto offset = global_state->block_states[bid].contact_local_buff_size;

  const auto chunk_size = (sample_size + blockDim.x - 1) / blockDim.x;
  const auto i0 = chunk_size * tid;
  const auto i1 = min(i0 + chunk_size, sample_size);

  for (auto i = i0; i < i1; ++i) {
    const auto rev_idx = rev_unit_idx_buff[(bid * buff_alignment) + i];
    const auto fwd_idx = lef_rev_unit_idx_buff[rev_idx];

    const auto rev_pos = rev_pos_buff[rev_idx];
    const auto fwd_pos = fwd_pos_buff[fwd_idx];

    global_state->block_states[bid].contact_local_buff[offset + i].x =
        __float2uint_rn(static_cast<float>(min(rev_pos, fwd_pos)) / bin_size);
    global_state->block_states[bid].contact_local_buff[offset + i].y =
        __float2uint_rn(static_cast<float>(max(rev_pos, fwd_pos)) / bin_size);
  }

  __syncthreads();
  if (tid == 0) {
    global_state->block_states[bid].contact_local_buff_size += sample_size;
    if (global_state->block_states[bid].contact_local_buff_size ==
        global_state->block_states[bid].contact_local_buff_capacity) {
      global_state->block_states[blockIdx.x].simulation_completed = true;
      atomicAdd(&global_state->ntasks_completed, 1);
    }
  }
}

__global__ void generate_moves(GlobalStateDev* global_state) {
  if (global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  auto local_rng_state = global_state->block_states[bid].rng_state[tid];
  const auto num_active_lefs = global_state->block_states[bid].num_active_lefs;

  const auto rev_mean = config.rev_extrusion_speed;
  const auto rev_stddev = config.rev_extrusion_speed_std;
  const auto fwd_mean = config.fwd_extrusion_speed;
  const auto fwd_stddev = config.fwd_extrusion_speed_std;

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  for (auto i = i0; i < i1; ++i) {
    const auto [n1, n2] = curand_normal2(&local_rng_state);
    global_state->block_states[bid].rev_moves_buff[i] =
        __float2uint_rn(static_cast<float>(rev_mean) + (n1 * rev_stddev));
    global_state->block_states[bid].fwd_moves_buff[i] =
        __float2uint_rn(static_cast<float>(fwd_mean) + (n2 * fwd_stddev));
  }

  global_state->block_states[bid].rng_state[tid] = local_rng_state;

  // Adjust moves
  __syncthreads();
  const auto chrom_start = global_state->tasks[bid].chrom_start;
  const auto chrom_end = global_state->tasks[bid].chrom_end;

  for (auto i = i1 - 1; i > 0; --i) {
    const auto pos1 = global_state->block_states[bid].rev_unit_pos[i - 1];
    const auto pos2 = global_state->block_states[bid].rev_unit_pos[i];
    assert(pos2 >= pos1);  // NOLINT
    if (pos1 == Simulation::EXTR_UNIT_IS_IDLE || pos2 == Simulation::EXTR_UNIT_IS_IDLE) {
      continue;
    }

    auto move1 = global_state->block_states[bid].rev_moves_buff[i - 1];
    auto move2 = global_state->block_states[bid].rev_moves_buff[i];

    if (pos1 < chrom_start + move1) {
      move1 = pos1 - chrom_start;
      atomicExch(global_state->block_states[bid].rev_moves_buff + i - 1, move1);
    }

    if (pos2 < chrom_start + move2) {
      move2 = pos2 - chrom_start;
      atomicExch(global_state->block_states[bid].rev_moves_buff + i, move2);
    }

    if (pos2 - move2 < pos1 - move1) {
      move1 = min(pos1 - chrom_start, (pos2 - pos1) + move2);
      atomicExch(global_state->block_states[bid].rev_moves_buff + i - 1, move1);
    } else if (i - 1 <= i0) {
      break;
    }
  }

  __syncthreads();
  for (auto i = i0; i < num_active_lefs - 1; ++i) {
    const auto pos1 = global_state->block_states[bid].fwd_unit_pos[i];
    const auto pos2 = global_state->block_states[bid].fwd_unit_pos[i + 1];
    assert(pos2 >= pos1);  // NOLINT
    if (pos1 == Simulation::EXTR_UNIT_IS_IDLE || pos2 == Simulation::EXTR_UNIT_IS_IDLE) {
      continue;
    }

    auto move1 = global_state->block_states[bid].fwd_moves_buff[i];
    auto move2 = global_state->block_states[bid].fwd_moves_buff[i + 1];

    if (pos1 + move1 >= chrom_end) {
      move1 = chrom_end - 1 - pos1;
      atomicExch(global_state->block_states[bid].fwd_moves_buff + i, move1);
    }

    if (pos2 + move2 >= chrom_end) {
      move2 = chrom_end - 1 - pos2;
      atomicExch(global_state->block_states[bid].fwd_moves_buff + i + 1, move2);
    }

    if (pos1 + move1 > pos2 + move2) {
      move2 = min(chrom_end - 1 - pos2, (pos1 + move1) - pos2);
      atomicExch(global_state->block_states[bid].fwd_moves_buff + i + 1, move2);
    } else if (i >= i1) {
      break;
    }
  }
}

__global__ void update_ctcf_states(GlobalStateDev* global_state) {
  if (global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  auto local_rng_state = global_state->block_states[bid].rng_state[tid];
  const auto nbarriers = global_state->tasks[bid].nbarriers;

  const auto chunk_size = (nbarriers + blockDim.x - 1) / blockDim.x;
  const auto i0 = tid * chunk_size;
  const auto i1 = min(i0 + chunk_size, nbarriers);

  for (auto i = i0; i < i1; ++i) {
    const auto current_state = global_state->block_states[bid].barrier_mask[i];
    const auto transition_prob =
        1.0F - (current_state == CTCF::OCCUPIED ? global_state->barrier_probs_occ_to_occ[i]
                                                : global_state->barrier_probs_nocc_to_nocc[i]);
    if (curand_uniform(&local_rng_state) < transition_prob) {
      global_state->block_states[bid].barrier_mask[i] =
          current_state == CTCF::OCCUPIED ? CTCF::NOT_OCCUPIED : CTCF::OCCUPIED;
    }
    __syncthreads();
  }
  global_state->block_states[bid].rng_state[tid] = local_rng_state;
}

__global__ void reset_collision_masks(GlobalStateDev* global_state) {
  if (global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  constexpr auto fill_value = Simulation::NO_COLLISIONS;

  if (tid == 0) {
    thrust::fill_n(thrust::device, global_state->block_states[bid].rev_collision_mask,
                   global_state->block_states[bid].num_active_lefs, fill_value);
    thrust::fill_n(thrust::device, global_state->block_states[bid].fwd_collision_mask,
                   global_state->block_states[bid].num_active_lefs, fill_value);
  }
}

__global__ void process_collisions(uint32_t current_epoch, GlobalStateDev* global_state) {
  if (global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  auto* block_state = global_state->block_states + bid;
  const auto* task = global_state->tasks + bid;

  if (tid == 0) {
    block_state->num_rev_units_at_5prime =
        modle::cu::dev::detect_collisions_at_5prime_single_threaded(
            block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->rev_moves_buff,
            block_state->rev_collision_mask, task->chrom_start, block_state->num_active_lefs);
  } else if (tid == 1) {
    block_state->num_fwd_units_at_3prime =
        modle::cu::dev::detect_collisions_at_3prime_single_threaded(
            block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->fwd_moves_buff,
            block_state->fwd_collision_mask, task->chrom_end, block_state->num_active_lefs);
  }
  __syncthreads();

  modle::cu::dev::detect_lef_bar_collisions(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->rev_moves_buff,
      block_state->fwd_moves_buff, block_state->num_active_lefs, global_state->barrier_pos,
      global_state->barrier_directions, block_state->barrier_mask,
      global_state->tasks[bid].nbarriers, block_state->rev_collision_mask,
      block_state->fwd_collision_mask, block_state->rng_state, block_state->num_rev_units_at_5prime,
      block_state->num_fwd_units_at_3prime);
  __syncthreads();

  modle::cu::dev::detect_primary_lef_lef_collisions(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->lef_rev_unit_idx,
      block_state->rev_moves_buff, block_state->fwd_moves_buff, block_state->num_active_lefs,
      global_state->barrier_pos, global_state->tasks[bid].nbarriers,
      block_state->rev_collision_mask, block_state->fwd_collision_mask, block_state->rng_state,
      block_state->num_rev_units_at_5prime, block_state->num_fwd_units_at_3prime);
  __syncthreads();

  modle::cu::dev::correct_moves_for_lef_bar_collisions(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->rev_moves_buff,
      block_state->fwd_moves_buff, block_state->num_active_lefs, global_state->barrier_pos,
      global_state->tasks[bid].nbarriers, block_state->rev_collision_mask,
      block_state->fwd_collision_mask);
  __syncthreads();

  modle::cu::dev::correct_moves_for_primary_lef_lef_collisions(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->rev_moves_buff,
      block_state->fwd_moves_buff, block_state->num_active_lefs, global_state->tasks[bid].nbarriers,
      block_state->rev_collision_mask, block_state->fwd_collision_mask);
  __syncthreads();

  modle::cu::dev::extrude(block_state->rev_unit_pos, block_state->fwd_unit_pos,
                          block_state->rev_moves_buff, block_state->fwd_moves_buff,
                          block_state->num_active_lefs, task->chrom_start, task->chrom_end);
  __syncthreads();
}

}  // namespace modle::cu::kernels
