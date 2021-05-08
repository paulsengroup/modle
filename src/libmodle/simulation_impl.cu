#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cuda/std/atomic>
#include <string>
#include <vector>

#include "modle/contacts.cuh"
#include "modle/extrusion_barriers.cuh"
#include "modle/simulation.cuh"

namespace modle::cu::Simulation {

__global__ void allocate_global_state_kernel(GlobalState* global_state, uint32_t max_nlefs,
                                             uint32_t max_nbarriers,
                                             cuda::std::atomic<bool>* global_ok) {
  if (threadIdx.x == 0 && *global_ok) {
    bool local_ok = true;
    const auto i = blockIdx.x;
    // Figure out why aliasing global_state->block_states doesn't work
    local_ok &= !!(global_state->block_states[i].barrier_mask = new bool[max_nbarriers]);

    local_ok &= !!(global_state->block_states[i].rev_unit_pos = new bp_t[max_nlefs]);
    local_ok &= !!(global_state->block_states[i].fwd_unit_pos = new bp_t[max_nlefs]);
    local_ok &= !!(global_state->block_states[i].lef_rev_unit_idx = new bp_t[max_nlefs]);
    local_ok &= !!(global_state->block_states[i].lef_fwd_unit_idx = new bp_t[max_nlefs]);

    local_ok &= !!(global_state->block_states[i].rev_moves_buff = new bp_t[max_nlefs]);
    local_ok &= !!(global_state->block_states[i].fwd_moves_buff = new bp_t[max_nlefs]);
    local_ok &= !!(global_state->block_states[i].rev_collision_mask = new bp_t[max_nlefs]);
    local_ok &= !!(global_state->block_states[i].fwd_collision_mask = new bp_t[max_nlefs]);

    local_ok &= !!(global_state->block_states[i].lef_unloader_affinities = new float[max_nlefs]);
    local_ok &= !!(global_state->block_states[i].contact_local_buff =
                       new uint2[max_nlefs]);  // TODO CHANGEME
    local_ok &= !!(global_state->block_states[i].epoch_buff = new bp_t[max_nlefs]);

    local_ok &=
        !!(global_state->block_states[i].rng_state = new curandStatePhilox4_32_10_t[blockDim.x]);

    if (!local_ok) {
      *global_ok = false;
    }
  }
}

__global__ void free_global_state_kernel(GlobalState* global_state) {
  if (threadIdx.x == 0) {
    const auto i = blockIdx.x;

    auto& block_state = global_state->block_states[i];

    delete[] block_state.barrier_mask;

    delete[] block_state.rev_unit_pos;
    delete[] block_state.fwd_unit_pos;
    delete[] block_state.lef_rev_unit_idx;
    delete[] block_state.lef_fwd_unit_idx;

    delete[] block_state.rev_moves_buff;
    delete[] block_state.fwd_moves_buff;
    delete[] block_state.rev_collision_mask;
    delete[] block_state.fwd_collision_mask;

    delete[] block_state.lef_unloader_affinities;
    delete[] block_state.contact_local_buff;
    delete[] block_state.epoch_buff;

    delete[] block_state.rng_state;
  }
}

__global__ void init_curand_kernel(GlobalState* global_state) {
  const auto seed = global_state->tasks[blockIdx.x].seed;
  auto rng_state = global_state->block_states[blockIdx.x].rng_state[threadIdx.x];

  curand_init(seed, threadIdx.x, 0, &rng_state);

  global_state->block_states[blockIdx.x].rng_state[threadIdx.x] = rng_state;
}

__device__ void reset_buffers_kernel(GlobalState* global_state, uint32_t nlefs,
                                     uint32_t nbarriers) {
  const auto id = threadIdx.x + blockIdx.x * blockDim.x;
  const auto nthreads = blockDim.x * gridDim.x;

  const auto chunk_size = max(1, (global_state->nblock_states + 1) / nthreads);

  const int default_pos = -1;
  const int default_idx = -1;
  const int default_mask = -1;

  const auto i0 = id * chunk_size;
  const auto i1 = (id + 1) * chunk_size;
  for (auto i = i0; i < i1 && i < global_state->nblock_states; ++i) {
    auto& block_state = global_state->block_states[i];

    memset(block_state.rev_unit_pos, default_pos, nlefs);
    memset(block_state.fwd_unit_pos, default_pos, nlefs);
    memset(block_state.lef_rev_unit_idx, default_idx, nlefs);
    memset(block_state.lef_fwd_unit_idx, default_idx, nlefs);

    memset(block_state.rev_moves_buff, 0U, nlefs);
    memset(block_state.fwd_moves_buff, 0U, nlefs);
    memset(block_state.rev_collision_mask, default_mask, nlefs);
    memset(block_state.fwd_collision_mask, default_mask, nlefs);
  }
}

[[nodiscard]] GlobalState* call_allocate_global_state_kernel(size_t grid_size, size_t block_size,
                                                             uint32_t max_nlefs,
                                                             uint32_t max_nbarriers) {
  try {
    CUDA_CALL(cudaDeviceSetLimit(cudaLimitMallocHeapSize, 1024ULL * 1024ULL * 4096LL));  // 2GB
    assert(grid_size > 0);                                                               // NOLINT
    // We allocate one state per block

    GlobalState* global_state_dev{nullptr};
    cuda::std::atomic<bool>* status_dev{nullptr};
    cuda::std::atomic<bool> status_host{true};

    CUDA_CALL(cudaMalloc(&status_dev, sizeof(cuda::std::atomic<bool>)));
    CUDA_CALL(cudaMemcpy(status_dev, &status_host, sizeof(cuda::std::atomic<bool>),
                         cudaMemcpyHostToDevice));

    CUDA_CALL(cudaMalloc(&global_state_dev, sizeof(GlobalState)));
    GlobalState global_state_host;
    CUDA_CALL(cudaMemcpy(&global_state_host, global_state_dev, sizeof(GlobalState),
                         cudaMemcpyDeviceToHost));

    BlockState* block_states_dev{nullptr};
    Task* tasks_dev{nullptr};
    bp_t* barrier_pos_dev{nullptr};
    dna::Direction* barrirer_directions_dev{nullptr};

    CUDA_CALL(cudaMalloc(&block_states_dev, grid_size * sizeof(BlockState)));
    CUDA_CALL(cudaMalloc(&tasks_dev, grid_size * sizeof(Task)));
    CUDA_CALL(cudaMalloc(&barrier_pos_dev, max_nbarriers * sizeof(bp_t)));
    CUDA_CALL(cudaMalloc(&barrirer_directions_dev, max_nbarriers * sizeof(dna::Direction)));

    global_state_host.ntasks = grid_size;
    global_state_host.nblock_states = grid_size;
    global_state_host.block_states = block_states_dev;
    global_state_host.tasks = tasks_dev;
    global_state_host.barrier_pos = barrier_pos_dev;
    global_state_host.barrier_directions = barrirer_directions_dev;

    CUDA_CALL(cudaMemcpy(global_state_dev, &global_state_host, sizeof(GlobalState),
                         cudaMemcpyHostToDevice));

    allocate_global_state_kernel<<<grid_size, block_size>>>(global_state_dev, max_nlefs,
                                                            max_nbarriers, status_dev);
    CUDA_CALL(cudaDeviceSynchronize());
    CUDA_CALL(cudaMemcpy(&status_host, status_dev, sizeof(cuda::std::atomic<bool>),
                         cudaMemcpyDeviceToHost));
    if (!status_host) {
      throw std::runtime_error("Unable to allocate enough memory to initialize device buffers.");
    }

    return global_state_dev;
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("call_allocate_global_state_kernel failed: {}"), e.what()));
  }
}

void call_free_global_state_kernel(size_t grid_size, size_t block_size,
                                   GlobalState* global_state_dev) {
  free_global_state_kernel<<<grid_size, block_size>>>(global_state_dev);
  GlobalState global_state_host;
  CUDA_CALL(cudaMemcpy(&global_state_host, global_state_dev, sizeof(GlobalState),
                       cudaMemcpyDeviceToHost));

  CUDA_CALL(cudaFree(global_state_host.block_states));
  CUDA_CALL(cudaFree(global_state_host.tasks));
  CUDA_CALL(cudaFree(global_state_host.barrier_pos));
  CUDA_CALL(cudaFree(global_state_host.barrier_directions));
  CUDA_CALL(cudaFree(global_state_dev));
}

void init_global_state(GlobalState* global_state, std::vector<Task>& host_tasks,
                       const std::vector<uint32_t>& barrier_pos,
                       const std::vector<dna::Direction>& barrier_dir) {
  GlobalState tmp_global_state;
  CUDA_CALL(
      cudaMemcpy(&tmp_global_state, global_state, sizeof(GlobalState), cudaMemcpyDeviceToHost));
  tmp_global_state.nblock_states = static_cast<uint32_t>(host_tasks.size());
  tmp_global_state.ntasks = static_cast<uint32_t>(host_tasks.size());
  CUDA_CALL(
      cudaMemcpy(global_state, &tmp_global_state, sizeof(GlobalState), cudaMemcpyHostToDevice));
  auto* dev_tasks = tmp_global_state.tasks;

  CUDA_CALL(cudaMemcpy(dev_tasks, host_tasks.data(), host_tasks.size() * sizeof(Task),
                       cudaMemcpyHostToDevice));
}

void call_simulation_kernel(size_t grid_size, size_t block_size, GlobalState* global_state,
                            std::vector<Task>& host_tasks, const std::vector<uint32_t>& barrier_pos,
                            const std::vector<dna::Direction>& barrier_dir) {
  assert(barrier_pos.size() == barrier_dir.size());  // NOLINT
  assert(host_tasks.size() == grid_size);            // NOLINT
  // TODO handle case where nlefs and nbarriers can be different across states
  const auto nlefs = host_tasks.front().nlefs;
  const auto nbarriers = barrier_pos.size();

  init_global_state(global_state, host_tasks, barrier_pos, barrier_dir);

  init_curand_kernel<<<grid_size, block_size>>>(global_state);

  // run_simulation_kernel<<<grid_size, block_size>>>(states, grid_size * block_size, nlefs,
  //                                                 nbarriers);
}

}  // namespace modle::cu::Simulation
