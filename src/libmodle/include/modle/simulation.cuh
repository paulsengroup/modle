#pragma once

#include <cstdint>
#include <cuda/std/atomic>
#include <vector>

#include "modle/common.hpp"
#include "modle/contacts.cuh"
#include "modle/extrusion_barriers.cuh"

using curandStatePhilox4_32_10_t = struct curandStatePhilox4_32_10;

namespace modle::cu::Simulation {
__global__ void mykernel();
__global__ void mykernel2(modle::cu::ContactMatrix<uint32_t>* m, cuda::std::atomic<uint32_t>* buff,
                          size_t nrows, size_t ncols, size_t& missed_updates, size_t& tot_contacts);

__global__ void mykernel3(const ExtrusionBarrier* barriers, size_t nbarriers, CTCF::State* mask,
                          curandStatePhilox4_32_10_t* rng_state);
[[nodiscard]] std::vector<uint32_t> run_mykernel2(size_t nrows, size_t ncols,
                                                  size_t& missed_updates, size_t& tot_contacts);

[[nodiscard]] std::vector<CTCF::State> run_mykernel3(
    const std::vector<ExtrusionBarrier>& host_barriers, uint64_t seed = 123456789);

struct Task {
  __host__ __device__ inline Task() = default;
  uint32_t id{};
  char* chrom_name{nullptr};
  bp_t chrom_start{};
  bp_t chrom_end{};
  uint32_t cell_id{};

  uint32_t n_target_epochs{};
  uint64_t n_target_contacts{};

  uint32_t nlefs{};
  uint32_t nbarriers{};
  uint64_t seed{};
};

struct State : Task {
  __host__ __device__ inline State() = default;
  __host__ __device__ inline State(const State& other) = delete;
  __host__ __device__ inline State(State&& other) = delete;
  __host__ __device__ inline ~State();

  bp_t* barrier_pos{nullptr};
  dna::Direction* barrier_directions{nullptr};
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

[[nodiscard]] State* call_init_global_buffers_kernel(size_t grid_size, size_t block_size,
                                                     uint32_t max_nlefs, uint32_t max_nbarriers);
void call_free_global_buffers_kernel(size_t grid_size, size_t block_size, State* states);

void call_simulation_kernel(size_t grid_size, size_t block_size, const std::vector<Task>& tasks,
                            const std::vector<uint32_t>& barrier_pos,
                            const std::vector<dna::Direction>& barrier_dir);

}  // namespace modle::cu::Simulation
