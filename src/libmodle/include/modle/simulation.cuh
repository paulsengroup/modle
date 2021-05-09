#pragma once

#include <cstdint>
#include <cuda/std/atomic>
#include <vector>

#include "modle/common.cuh"
#include "modle/config.cuh"
#include "modle/contacts.cuh"
#include "modle/extrusion_barriers.cuh"

using curandStatePhilox4_32_10_t = struct curandStatePhilox4_32_10;

namespace modle::cu::Simulation {

struct Task {
  __host__ __device__ inline Task() = default;
  uint32_t id{};

  uint32_t chrom_id{};
  bp_t chrom_start{};
  bp_t chrom_end{};

  uint32_t cell_id{};
  uint32_t n_target_epochs{};
  uint64_t n_target_contacts{};

  uint32_t nlefs{};
  uint32_t nbarriers{};
  uint64_t seed{};
};

struct BlockState {
  __host__ __device__ inline BlockState() = default;
  __host__ __device__ inline BlockState(const BlockState& other) = delete;
  __host__ __device__ inline BlockState(BlockState&& other) = delete;

  bool* barrier_mask{nullptr};

  bp_t* rev_unit_pos{nullptr};
  bp_t* fwd_unit_pos{nullptr};
  uint32_t* lef_rev_unit_idx{nullptr};
  uint32_t* lef_fwd_unit_idx{nullptr};

  uint32_t* rev_moves_buff{nullptr};
  uint32_t* fwd_moves_buff{nullptr};
  collision_t* rev_collision_mask{nullptr};
  collision_t* fwd_collision_mask{nullptr};

  float* lef_unloader_affinities{nullptr};
  uint2* contact_local_buff{nullptr};
  uint32_t contact_local_buff_size{};
  uint32_t contact_local_buff_capacity{};
  uint32_t* epoch_buff{nullptr};

  curandStatePhilox4_32_10_t* rng_state{nullptr};
};

struct GlobalState {
  __host__ __device__ GlobalState() = default;
  __host__ __device__ inline GlobalState(const GlobalState& other) = delete;
  __host__ __device__ inline GlobalState(GlobalState&& other) = delete;

  BlockState* block_states{nullptr};
  Task* tasks{nullptr};

  bp_t* barrier_pos{nullptr};
  dna::Direction* barrier_directions{nullptr};

  uint32_t nblock_states{};
  uint32_t ntasks{};

  uint32_t* large_uint_buff1{nullptr};
  uint32_t* large_uint_buff2{nullptr};
};

void init_simulation_params(const Config& c);

[[nodiscard]] GlobalState* call_allocate_global_state_kernel(size_t grid_size, size_t block_size,
                                                             uint32_t max_nlefs,
                                                             uint32_t max_nbarriers);

void call_free_global_state_kernel(size_t grid_size, size_t block_size, GlobalState* global_state);

void call_simulation_kernel(size_t grid_size, size_t block_size, GlobalState* global_state,
                            std::vector<Task>& tasks_host, const std::vector<uint32_t>& barrier_pos,
                            const std::vector<dna::Direction>& barrier_dir);

}  // namespace modle::cu::Simulation
