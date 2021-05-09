#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <fmt/format.h>
//#include <thrust/execution_policy.h>
//#include <thrust/sort.h>
#include <thrust/device_vector.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cub/cub.cuh>
#include <cuda/std/atomic>
#include <numeric>
#include <string>
#include <vector>

#include "modle/config.cuh"
#include "modle/contacts.cuh"
#include "modle/extrusion_barriers.cuh"
#include "modle/simulation.cuh"

namespace modle::cu::Simulation {

__constant__ __device__ struct Config config;

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

  curand_init(seed, threadIdx.x + blockIdx.x * blockDim.x, 0, &rng_state);

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

__global__ void generate_initial_loading_epochs(GlobalState* global_state, uint32_t num_epochs) {
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;
  const auto id = tid + bid * blockDim.x;
  const auto nthreads = blockDim.x * gridDim.x;

  auto* block_state = global_state->block_states + bid;
  const auto scaling_factor = static_cast<float>(4 * config.average_lef_lifetime);
  auto rng_state_local = block_state->rng_state[tid];

  auto local_buff_size = 0U;
  const auto local_buff_capacity = 32U;
  uint32_t local_buff[local_buff_capacity];

  auto* global_buff = global_state->large_uint_buff1;

  const auto chunk_size = (num_epochs + nthreads - 1) / nthreads;
  const auto i0 = id * chunk_size;
  const auto i1 = min(i0 + chunk_size, num_epochs);

  for (auto i = i0; i < i1; ++i) {
    local_buff[local_buff_size++] =
        __float2uint_rn(curand_uniform(&rng_state_local) * scaling_factor);

    if (local_buff_size == local_buff_capacity) {
      memcpy(global_buff + i - local_buff_size, local_buff, local_buff_size * sizeof(uint32_t));
      local_buff_size = 0;
    }
  }

  memcpy(global_buff + i1 - local_buff_size, local_buff, local_buff_size * sizeof(uint32_t));
  block_state->rng_state[tid] = rng_state_local;
}

void init_simulation_params(const Config& c) {
  CUDA_CALL(cudaMemcpyToSymbol(config, &c, sizeof(modle::cu::Config)));
}

[[nodiscard]] GlobalState* call_allocate_global_state_kernel(size_t grid_size, size_t block_size,
                                                             uint32_t max_nlefs,
                                                             uint32_t max_nbarriers) {
  try {
    CUDA_CALL(cudaDeviceSetLimit(cudaLimitMallocHeapSize, 1024ULL * 1024ULL * 2048LL));  // 2GB
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

    CUDA_CALL(cudaMalloc(&global_state_host.block_states, grid_size * sizeof(BlockState)));
    CUDA_CALL(cudaMalloc(&global_state_host.tasks, grid_size * sizeof(Task)));
    CUDA_CALL(cudaMalloc(&global_state_host.barrier_pos, max_nbarriers * sizeof(bp_t)));
    CUDA_CALL(
        cudaMalloc(&global_state_host.barrier_directions, max_nbarriers * sizeof(dna::Direction)));
    CUDA_CALL(cudaMalloc(&global_state_host.large_uint_buff1,
                         std::max(max_nlefs, max_nbarriers) * grid_size * sizeof(uint32_t)));
    CUDA_CALL(cudaMalloc(&global_state_host.large_uint_buff2,
                         std::max(max_nlefs, max_nbarriers) * grid_size * sizeof(uint32_t)));

    global_state_host.ntasks = grid_size;
    global_state_host.nblock_states = grid_size;

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
  CUDA_CALL(cudaFree(global_state_host.large_uint_buff1));
  CUDA_CALL(cudaFree(global_state_host.large_uint_buff2));
  CUDA_CALL(cudaFree(global_state_dev));
}

void init_global_state(size_t grid_size, size_t block_size, GlobalState* global_state_dev,
                       std::vector<Task>& host_tasks, const std::vector<uint32_t>& barrier_pos,
                       const std::vector<dna::Direction>& barrier_dir) {
  GlobalState global_state_host;
  CUDA_CALL(cudaMemcpy(&global_state_host, global_state_dev, sizeof(GlobalState),
                       cudaMemcpyDeviceToHost));
  global_state_host.nblock_states = static_cast<uint32_t>(host_tasks.size());
  global_state_host.ntasks = static_cast<uint32_t>(host_tasks.size());
  CUDA_CALL(cudaMemcpy(global_state_dev, &global_state_host, sizeof(GlobalState),
                       cudaMemcpyHostToDevice));

  CUDA_CALL(cudaMemcpy(global_state_host.tasks, host_tasks.data(), host_tasks.size() * sizeof(Task),
                       cudaMemcpyHostToDevice));

  init_curand_kernel<<<grid_size, block_size>>>(global_state_dev);
}

template <typename T>
__global__ void shift_and_scatter_lef_loading_epochs(const T* global_buff, BlockState* block_states,
                                                     const uint32_t* begin_offsets,
                                                     const uint32_t* end_offsets) {
  const auto tid = threadIdx.x;
  const auto bid = blockIdx.x;

  const auto* source_ptr = global_buff + begin_offsets[bid];
  auto* dest_ptr = block_states[bid].epoch_buff;
  const auto size = end_offsets[bid] - begin_offsets[bid];

  if (threadIdx.x == 0) {
    memcpy(dest_ptr, source_ptr, size * sizeof(T));
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

void setup_burnin_phase(size_t grid_size, size_t block_size, GlobalState* global_state_dev,
                        const std::vector<Task>& tasks_host) {
  GlobalState global_state_host;
  CUDA_CALL(cudaMemcpy(&global_state_host, global_state_dev, sizeof(GlobalState),
                       cudaMemcpyDeviceToHost));

  const auto num_epochs_for_initial_loading =
      std::accumulate(tasks_host.begin(), tasks_host.end(), 0UL,
                      [](auto accumulator, const auto& task) { return accumulator + task.nlefs; });
  generate_initial_loading_epochs<<<grid_size, block_size>>>(global_state_dev,
                                                             num_epochs_for_initial_loading);

  // Generate segment offsets
  thrust::device_vector<uint32_t> offsets_dev;
  {
    thrust::host_vector<uint32_t> offsets_host(tasks_host.size() + 1);
    offsets_host.front() = 0;
    std::transform(tasks_host.begin(), tasks_host.end(), offsets_host.begin() + 1,
                   [offset = 0U](const auto& task) mutable { return (offset += task.nlefs); });
    offsets_dev = offsets_host;
  }
  auto* d_begin_offsets = thrust::device_pointer_cast(offsets_dev.data()).get();
  auto* d_end_offsets = d_begin_offsets + 1;

  auto* d_key_buf = global_state_host.large_uint_buff1;
  auto* d_key_alt_buf = global_state_host.large_uint_buff2;
  const auto num_items_to_sort =
      std::accumulate(tasks_host.begin(), tasks_host.end(), 0UL,
                      [](auto accumulator, const auto& task) { return accumulator + task.nlefs; });
  const auto num_segments_to_sort = tasks_host.size();
  size_t temp_storage_bytes = 0;

  cub::DoubleBuffer<uint32_t> d_keys(d_key_buf, d_key_alt_buf);  // Compute temp_storage size
  CUDA_CALL(cub::DeviceSegmentedRadixSort::SortKeys(nullptr, temp_storage_bytes, d_keys,
                                                    num_items_to_sort, num_segments_to_sort,
                                                    d_begin_offsets, d_end_offsets));
  thrust::device_vector<uint8_t> d_tmp_storage(temp_storage_bytes);  // Allocate temp_storage
  auto* d_tmp_storage_ptr = thrust::device_pointer_cast(d_tmp_storage.data()).get();

  CUDA_CALL(cub::DeviceSegmentedRadixSort::SortKeys(d_tmp_storage_ptr, temp_storage_bytes, d_keys,
                                                    num_items_to_sort, num_segments_to_sort,
                                                    d_begin_offsets, d_end_offsets));

  shift_and_scatter_lef_loading_epochs<<<grid_size, block_size>>>(
      d_key_buf, global_state_host.block_states, d_begin_offsets, d_end_offsets);
}

void call_simulation_kernel(size_t grid_size, size_t block_size, GlobalState* global_state,
                            std::vector<Task>& tasks_host, const std::vector<uint32_t>& barrier_pos,
                            const std::vector<dna::Direction>& barrier_dir) {
  assert(barrier_pos.size() == barrier_dir.size());  // NOLINT
  assert(tasks_host.size() == grid_size);            // NOLINT
  // TODO handle case where nlefs and nbarriers can be different across states
  const auto nlefs = tasks_host.front().nlefs;
  const auto nbarriers = barrier_pos.size();

  init_global_state(grid_size, block_size, global_state, tasks_host, barrier_pos, barrier_dir);

  setup_burnin_phase(grid_size, block_size, global_state, tasks_host);

  // simulate_kernel<<<grid_size, block_size>>>(global_state);
}

}  // namespace modle::cu::Simulation
