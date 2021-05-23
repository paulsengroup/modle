#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>

#include <cstdio>
#include <cub/cub.cuh>

#include "modle/cu/simulation.hpp"
#include "modle/cu/simulation_internal.cuh"

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
    block_state.epoch_idx = 0;
  }
  __syncthreads();

  const auto chunk_size = (nlefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = min(tid * chunk_size, nlefs);
  const auto i1 = min(i0 + chunk_size, nlefs);

  if (i0 < i1) {
    thrust::sequence(thrust::seq, block_state.lef_rev_unit_idx + i0,
                     block_state.lef_rev_unit_idx + i1, i0);
    thrust::sequence(thrust::seq, block_state.lef_fwd_unit_idx + i0,
                     block_state.lef_fwd_unit_idx + i1, i0);
  }
}

__global__ void compute_initial_loading_epochs(GlobalStateDev* global_state) {
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  auto* block_state = global_state->block_states + bid;
  const auto* task = global_state->tasks + bid;

  const auto nlefs = task->nlefs;

  // Here we are using a temporary buffer to generate and sort epochs
  modle::cu::dev::generate_initial_loading_epochs(block_state->tmp_lef_buff1,
                                                  block_state->rng_state, nlefs);
  __syncthreads();

  // Sort epochs in ascending order
  if (tid == 0) {
    cub::DoubleBuffer keys_buff(block_state->tmp_lef_buff1, block_state->tmp_lef_buff2);
    cub::DeviceRadixSort::SortKeys(block_state->cub_tmp_storage, block_state->cub_tmp_storage_bytes,
                                   keys_buff, nlefs);
    cudaDeviceSynchronize();
    const auto status = cudaGetLastError();
    assert(status == cudaSuccess);  // NOLINT TODO: Do proper error handling
    block_state->tmp_lef_buff1 = keys_buff.Current();
    block_state->tmp_lef_buff2 = keys_buff.Alternate();
  }
  __syncthreads();
  const auto offset = *block_state->tmp_lef_buff1;
  __syncthreads();

  // Make epochs start from 0
  if (offset != 0) {
    const auto chunk_size = (nlefs + blockDim.x - 1) / blockDim.x;
    const auto i0 = min(tid * chunk_size, nlefs);
    const auto i1 = min(i0 + chunk_size, nlefs);

    if (i0 != i1) {
      thrust::transform(thrust::seq, block_state->tmp_lef_buff1 + i0,
                        block_state->tmp_lef_buff1 + i1, block_state->tmp_lef_buff1 + i0,
                        [&](const auto n) {
                          assert(offset <= n);  // NOLINT
                          return n - offset;
                        });
    }
  }
  __syncthreads();

  // Encode epochs using a run-length encoding scheme
  if (tid == 0) {
    cub::DeviceRunLengthEncode::Encode(
        block_state->cub_tmp_storage, block_state->cub_tmp_storage_bytes,
        block_state->tmp_lef_buff1, block_state->loading_epochs,
        block_state->lefs_to_load_per_epoch, &block_state->num_unique_loading_epochs, nlefs);
    cudaDeviceSynchronize();
    const auto status = cudaGetLastError();
    assert(status == cudaSuccess);  // NOLINT TODO: Do proper error handling
  }
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

__global__ void bind_and_sort_lefs(uint32_t current_epoch, GlobalStateDev* global_state) {
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;
  auto* block_state = global_state->block_states + bid;
  const auto* task = global_state->tasks + bid;
  const auto simulation_completed = block_state->simulation_completed;

  if (!simulation_completed) {
    modle::cu::dev::select_and_bind_lefs(
        block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->lef_rev_unit_idx,
        block_state->lef_fwd_unit_idx, block_state->loading_epochs,
        block_state->lefs_to_load_per_epoch, block_state->epoch_idx,
        block_state->num_unique_loading_epochs, block_state->num_active_lefs, current_epoch,
        task->nlefs, task->chrom_start, task->chrom_end, block_state->rng_state,
        block_state->burnin_completed);
  }
  __syncthreads();

  if (!simulation_completed && tid == 0) {  // NOLINT
    // Sort LEFs by pos. TODO: take into account binding order & other factors in case of a tie
    cub::DoubleBuffer keys_tmp_buff(block_state->rev_unit_pos, block_state->tmp_lef_buff1);
    cub::DoubleBuffer vals_tmp_buff(block_state->lef_rev_unit_idx, block_state->tmp_lef_buff2);

    cub::DeviceRadixSort::SortPairs(block_state->cub_tmp_storage,
                                    block_state->cub_tmp_storage_bytes, keys_tmp_buff,
                                    vals_tmp_buff, block_state->num_active_lefs);
    cudaDeviceSynchronize();
    const auto status = cudaGetLastError();
    assert(status == cudaSuccess);  // NOLINT TODO: Do proper error handling
    /* Assigning pointer in this way doesn't work, as if Current() happens to point to *_buff?, then
    the position of idle LEFs will be overwritten by whatever happens to be in the tmp buffers
    block_state->rev_unit_pos = keys_tmp_buff.Current();
    block_state->tmp_lef_buff1 = keys_tmp_buff.Alternate();
    block_state->lef_rev_unit_idx = vals_tmp_buff.Current();
    block_state->tmp_lef_buff2 = vals_tmp_buff.Alternate();
     */
    memcpy(block_state->rev_unit_pos, keys_tmp_buff.Current(),
           block_state->num_active_lefs * sizeof(bp_t));
    memcpy(block_state->lef_rev_unit_idx, vals_tmp_buff.Current(),
           block_state->num_active_lefs * sizeof(uint32_t));
  }

  __syncthreads();
  if (!simulation_completed) {
    modle::cu::dev::update_extr_unit_mappings(block_state->lef_rev_unit_idx,
                                              block_state->lef_fwd_unit_idx,
                                              block_state->num_active_lefs, dna::Direction::rev);
  }

  __syncthreads();
  if (!simulation_completed && tid == 0) {  // NOLINT
    cub::DoubleBuffer keys_tmp_buff(block_state->fwd_unit_pos, block_state->tmp_lef_buff1);
    cub::DoubleBuffer vals_tmp_buff(block_state->lef_fwd_unit_idx, block_state->tmp_lef_buff2);

    cub::DeviceRadixSort::SortPairs(block_state->cub_tmp_storage,
                                    block_state->cub_tmp_storage_bytes, keys_tmp_buff,
                                    vals_tmp_buff, block_state->num_active_lefs);
    cudaDeviceSynchronize();
    const auto status = cudaGetLastError();
    assert(status == cudaSuccess);  // NOLINT TODO: Do proper error handling
    /*
    block_state->fwd_unit_pos = keys_tmp_buff.Current();
    block_state->tmp_lef_buff1 = keys_tmp_buff.Alternate();
    block_state->lef_fwd_unit_idx = vals_tmp_buff.Current();
    block_state->tmp_lef_buff2 = vals_tmp_buff.Alternate();
     */
    memcpy(block_state->fwd_unit_pos, keys_tmp_buff.Current(),
           block_state->num_active_lefs * sizeof(bp_t));
    memcpy(block_state->lef_fwd_unit_idx, vals_tmp_buff.Current(),
           block_state->num_active_lefs * sizeof(uint32_t));
  }

  __syncthreads();
  if (!simulation_completed) {
    modle::cu::dev::update_extr_unit_mappings(block_state->lef_rev_unit_idx,
                                              block_state->lef_fwd_unit_idx,
                                              block_state->num_active_lefs, dna::Direction::fwd);
  }
}

__global__ void select_lefs_then_register_contacts(GlobalStateDev* global_state) {
  if (!global_state->block_states[blockIdx.x].burnin_completed ||
      global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  auto* block_state = global_state->block_states + bid;
  const auto* task = global_state->tasks + bid;

  auto* shuffled_idx = block_state->tmp_lef_buff1;
  if (block_state->burnin_completed) {
    modle::cu::dev::shuffle_lefs(shuffled_idx, block_state->tmp_lef_buff2,
                                 block_state->tmp_lef_buff3, block_state->tmp_lef_buff4,
                                 block_state->cub_tmp_storage, block_state->cub_tmp_storage_bytes,
                                 block_state->num_active_lefs, block_state->rng_state);
  }
  __syncthreads();

  __shared__ uint32_t sample_size, contact_buff_size_shared;
  if (tid == 0) {
    const auto mean_sample_size = static_cast<double>(config.lef_fraction_contact_sampling) *
                                  static_cast<double>(block_state->num_active_lefs);

    sample_size = min(curand_poisson(block_state->rng_state + tid, mean_sample_size),
                      block_state->num_active_lefs);
    atomicExch(&contact_buff_size_shared, block_state->contact_local_buff_size);
  }
  __syncthreads();
  if (sample_size == 0) {
    return;
  }

  const auto bin_size = static_cast<float>(config.bin_size);
  const auto chrom_start = task->chrom_start;
  const auto chrom_end = task->chrom_end;

  modle::cu::dev::register_contacts(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, shuffled_idx,
      block_state->lef_rev_unit_idx, block_state->contact_local_buff, chrom_start, chrom_end,
      bin_size, sample_size, block_state->num_active_lefs, &contact_buff_size_shared,
      block_state->contact_local_buff_capacity);

  __syncthreads();
  if (tid == 0) {
    global_state->block_states[bid].contact_local_buff_size = contact_buff_size_shared;
    if (global_state->block_states[bid].contact_local_buff_size ==
        global_state->block_states[bid].contact_local_buff_capacity) {
      global_state->block_states[blockIdx.x].simulation_completed = true;
      atomicAdd(&global_state->ntasks_completed, 1);
      __threadfence_block();
    }
  }
}

__global__ void generate_moves(GlobalStateDev* global_state) {
  const auto bid = blockIdx.x;
  const auto num_active_lefs = global_state->block_states[bid].num_active_lefs;
  if (global_state->block_states[bid].simulation_completed || num_active_lefs == 0) {
    return;
  }
  const auto tid = threadIdx.x;

  auto local_rng_state = global_state->block_states[bid].rng_state[tid];

  const auto rev_mean = config.rev_extrusion_speed;
  const auto rev_stddev = config.rev_extrusion_speed_std;
  const auto fwd_mean = config.fwd_extrusion_speed;
  const auto fwd_stddev = config.fwd_extrusion_speed_std;

  const auto chunk_size = (num_active_lefs + blockDim.x - 1) / blockDim.x;
  const auto i0 = min(tid * chunk_size, num_active_lefs);
  const auto i1 = min(i0 + chunk_size, num_active_lefs);

  for (auto i = i0; i < i1; ++i) {
    const auto [n1, n2] = curand_normal2(&local_rng_state);
    global_state->block_states[bid].rev_moves_buff[i] =
        static_cast<bp_t>(roundf(static_cast<float>(rev_mean) + (n1 * rev_stddev)));
    global_state->block_states[bid].fwd_moves_buff[i] =
        static_cast<bp_t>(roundf(static_cast<float>(fwd_mean) + (n2 * fwd_stddev)));
  }

  global_state->block_states[bid].rng_state[tid] = local_rng_state;

  // Adjust moves
  __syncthreads();
  const auto chrom_start = global_state->tasks[bid].chrom_start;
  const auto chrom_end = global_state->tasks[bid].chrom_end;
  if (i0 != i1) {
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
  }

  __syncthreads();
  if (i0 != i1) {
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

__global__ void process_collisions(GlobalStateDev* global_state) {
  if (global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  auto* block_state = global_state->block_states + bid;
  const auto* task = global_state->tasks + bid;

  modle::cu::dev::reset_collision_masks(block_state->rev_collision_mask,
                                        block_state->fwd_collision_mask,
                                        block_state->num_active_lefs);
  __syncthreads();
  if (tid == 0) {
    // printf("%d reset_collision_masks_completed\n", bid);
  }

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
  if (tid == 0) {
    // printf("%d detect_collisions_at_?prime_completed\n", bid);
  }

  modle::cu::dev::detect_lef_bar_collisions(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->rev_moves_buff,
      block_state->fwd_moves_buff, block_state->num_active_lefs, global_state->barrier_pos,
      global_state->barrier_directions, block_state->barrier_mask,
      global_state->tasks[bid].nbarriers, block_state->rev_collision_mask,
      block_state->fwd_collision_mask, block_state->rng_state, block_state->num_rev_units_at_5prime,
      block_state->num_fwd_units_at_3prime);
  __syncthreads();
  if (tid == 0) {
    // printf("%d detect_lef_bar_collisions_completed\n", bid);
  }

  modle::cu::dev::detect_primary_lef_lef_collisions(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->lef_rev_unit_idx,
      block_state->rev_moves_buff, block_state->fwd_moves_buff, block_state->num_active_lefs,
      global_state->barrier_pos, global_state->tasks[bid].nbarriers,
      block_state->rev_collision_mask, block_state->fwd_collision_mask, block_state->rng_state,
      block_state->num_rev_units_at_5prime, block_state->num_fwd_units_at_3prime);
  __syncthreads();
  if (tid == 0) {
    // printf("%d detect_primary_lef_lef_collisions_completed\n", bid);
  }

  modle::cu::dev::correct_moves_for_lef_bar_collisions(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->rev_moves_buff,
      block_state->fwd_moves_buff, block_state->num_active_lefs, global_state->barrier_pos,
      global_state->tasks[bid].nbarriers, block_state->rev_collision_mask,
      block_state->fwd_collision_mask);
  __syncthreads();
  if (tid == 0) {
    // printf("%d correct_moves_for_lef_bar_collisions_completed\n", bid);
  }

  modle::cu::dev::correct_moves_for_primary_lef_lef_collisions(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->rev_moves_buff,
      block_state->fwd_moves_buff, block_state->num_active_lefs, global_state->tasks[bid].nbarriers,
      block_state->rev_collision_mask, block_state->fwd_collision_mask);
  __syncthreads();
  if (tid == 0) {
    //  printf("%d correct_moves_for_primary_lef_lef_collisions_completed\n", bid);
  }

  modle::cu::dev::detect_secondary_lef_lef_collisions(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->rev_moves_buff,
      block_state->fwd_moves_buff, block_state->num_active_lefs, global_state->tasks[bid].nbarriers,
      block_state->rev_collision_mask, block_state->fwd_collision_mask, block_state->rng_state,
      block_state->num_rev_units_at_5prime, block_state->num_fwd_units_at_3prime);
  __syncthreads();
  if (tid == 0) {
    // printf("%d detect_secondary_lef_lef_collisions_completed\n", bid);
  }

  modle::cu::dev::correct_moves_for_secondary_lef_lef_collisions(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->rev_moves_buff,
      block_state->fwd_moves_buff, block_state->num_active_lefs, global_state->tasks[bid].nbarriers,
      block_state->rev_collision_mask, block_state->fwd_collision_mask);
  __syncthreads();
  if (tid == 0) {
    // printf("%d correct_moves_for_secondary_lef_lef_collisions_completed\n", bid);
  }
}

__global__ void extrude_and_release_lefs(GlobalStateDev* global_state) {
  if (global_state->block_states[blockIdx.x].simulation_completed) {
    return;
  }
  const auto bid = blockIdx.x;

  auto* block_state = global_state->block_states + bid;
  const auto* task = global_state->tasks + bid;

  modle::cu::dev::extrude(block_state->rev_unit_pos, block_state->fwd_unit_pos,
                          block_state->rev_moves_buff, block_state->fwd_moves_buff,
                          block_state->num_active_lefs, task->chrom_start, task->chrom_end);

  modle::cu::dev::generate_lef_unloader_affinities(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, global_state->barrier_directions,
      block_state->lef_rev_unit_idx, block_state->rev_collision_mask,
      block_state->fwd_collision_mask, block_state->num_active_lefs,
      global_state->tasks[bid].nbarriers, block_state->lef_unloader_affinities);
  __syncthreads();

  modle::cu::dev::select_and_release_lefs(
      block_state->rev_unit_pos, block_state->fwd_unit_pos, block_state->lef_rev_unit_idx,
      block_state->lef_fwd_unit_idx, block_state->num_active_lefs,
      block_state->lef_unloader_affinities, block_state->lef_unloader_affinities_prefix_sum,
      block_state->cub_tmp_storage, block_state->cub_tmp_storage_bytes, block_state->rng_state);
  __syncthreads();
}

}  // namespace modle::cu::kernels
