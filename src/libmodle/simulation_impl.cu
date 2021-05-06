#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "modle/contacts.cuh"
#include "modle/extrusion_barriers.cuh"
#include "modle/simulation.cuh"

namespace modle::cu::Simulation {

__global__ void setup_curand_kernel(curandStatePhilox4_32_10_t* state, uint64_t seed = 123456789) {
  const auto id = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, id, 0, &state[id]);
}

__global__ void init_barriers(ExtrusionBarrier* barriers, size_t nbarriers, const size_t* positions,
                              const dna::Direction* directions) {
  const auto id = threadIdx.x + blockIdx.x * blockDim.x;
  const auto chunk_size = nbarriers / blockDim.x;

  const auto i0 = id * chunk_size;
  const auto i1 = i0 + chunk_size;

  for (auto i = i0; i < i1 && i < nbarriers; ++i) {
    barriers[i] = ExtrusionBarrier{positions[i], 0.93, 0.7, directions[i]};
  }
  __syncthreads();
}

__global__ void init_ctcf_states_kernel(CTCF::State* mask, size_t nctcfs,
                                        curandStatePhilox4_32_10_t* rng_state,
                                        double pblock = 0.85) {
  const auto id = threadIdx.x + blockIdx.x * blockDim.x;
  const auto chunk_size = (nctcfs + 1) / blockDim.x;

  auto local_state = rng_state[id];

  const auto i0 = id * chunk_size;
  const auto i1 = i0 + chunk_size;
  auto i = i0;

  do {
    const auto buff = curand_uniform2_double(&local_state);
    mask[i++] = buff.x < pblock ? CTCF::OCCUPIED : CTCF::NOT_OCCUPIED;
    if (i < i1) {
      mask[i++] = buff.y < pblock ? CTCF::OCCUPIED : CTCF::NOT_OCCUPIED;
    }
  } while (i < i1 && i < nctcfs);

  rng_state[id] = local_state;
  __syncthreads();
}

__global__ void mykernel() {
  const auto id = threadIdx.x + blockIdx.x * blockDim.x;
  printf("id=%d\n", id);
}

__global__ void mykernel2(ContactMatrix<uint32_t>* m, cuda::std::atomic<uint32_t>* buff,
                          size_t nrows, size_t ncols) {
  const auto id = threadIdx.x + blockIdx.x * blockDim.x;
  if (id == 0) {
    m->reset(buff, nrows, ncols);
  }
  __syncthreads();
  m->add(id, id, static_cast<uint32_t>(id));
}

__global__ void mykernel3(const ExtrusionBarrier* barriers, size_t nbarriers, CTCF::State* mask,
                          curandStatePhilox4_32_10_t* rng_state) {
  const auto id = threadIdx.x + blockIdx.x * blockDim.x;
  for (auto i = 0UL; i < 100000; ++i) {
    if (id == 0 && i % 1000 == 0) printf("iter=%lu\n", i);
    CTCF::update_states(barriers, nbarriers, mask, rng_state);
  }
}

std::vector<uint32_t> run_mykernel2(size_t nrows, size_t ncols, size_t& missed_updates,
                                    size_t& tot_contacts) {
  cuda::std::atomic<uint32_t>* dev_buff{nullptr};
  modle::cu::ContactMatrix<uint32_t>* matrix{nullptr};

  if (const auto status = cudaMalloc(&matrix, sizeof(ContactMatrix<uint32_t>));
      status != cudaSuccess) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to allocate {} bytes of memory on the device"),
                    sizeof(ContactMatrix<uint32_t>)));
  }
  if (const auto status =
          cudaMalloc(&dev_buff, nrows * ncols * sizeof(cuda::std::atomic<uint32_t>));
      status != cudaSuccess) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to allocate {:.2f} MB of memory on the device"),
                    nrows * ncols * sizeof(cuda::std::atomic<uint32_t>)));
  }

  mykernel2<<<1, 10>>>(matrix, dev_buff, nrows, ncols);
  if (const auto status = cudaDeviceSynchronize(); status != cudaSuccess) {
    throw std::runtime_error("cudaDeviceSyncronize failed");
  }
  missed_updates = 0;
  tot_contacts = 0;

  if (const auto status = cudaFree(matrix); status != cudaSuccess) {
    throw std::runtime_error("cudaFree failed.");
  }

  std::vector<uint32_t> host_buff(nrows * ncols);
  if (const auto status =
          cudaMemcpy(host_buff.data(), dev_buff,
                     nrows * ncols * sizeof(cuda::std::atomic<uint32_t>), cudaMemcpyDeviceToHost);
      status != cudaSuccess) {
    throw std::runtime_error("cudaMemcpy failed.");
  }
  if (const auto status = cudaFree(dev_buff); status != cudaSuccess) {
    throw std::runtime_error("cudaFree failed.");
  }

  return host_buff;
}

std::vector<CTCF::State> run_mykernel3(const std::vector<ExtrusionBarrier>& host_barriers,
                                       uint64_t seed) {
  ExtrusionBarrier* dev_barriers{nullptr};
  CTCF::State* dev_barrier_states{nullptr};
  curandStatePhilox4_32_10_t* dev_rng_states{nullptr};
  std::vector<size_t> host_pos_buff(host_barriers.size());
  std::vector<dna::Direction> host_dir_buff(host_barriers.size());
  size_t* dev_pos_buff{nullptr};
  dna::Direction* dev_dir_buff{nullptr};

  const auto grid_size = 1UL;
  const auto block_size = 64UL;

  CUDA_CALL(cudaMalloc(&dev_barriers, sizeof(ExtrusionBarrier) * host_barriers.size()));
  CUDA_CALL(cudaMalloc(&dev_barrier_states, sizeof(CTCF::State) * host_barriers.size()));
  CUDA_CALL(cudaMalloc(&dev_rng_states, grid_size * block_size * sizeof(curandStatePhilox4_32_10_t)));
  CUDA_CALL(cudaMalloc(&dev_pos_buff, sizeof(size_t) * host_barriers.size()));
  CUDA_CALL(cudaMalloc(&dev_dir_buff, sizeof(CTCF::State) * host_barriers.size()));
  CUDA_CALL(cudaDeviceSynchronize());

  std::transform(host_barriers.begin(), host_barriers.end(), host_pos_buff.begin(),
                 [](const auto& b) { return b.pos(); });

  std::transform(host_barriers.begin(), host_barriers.end(), host_dir_buff.begin(),
                 [](const auto& b) {
                   return CTCF::major_blocking_dir_to_motif_dir(b.blocking_direction_major());
                 });

  CUDA_CALL(cudaMemcpy(dev_pos_buff, host_pos_buff.data(), sizeof(size_t) * host_barriers.size(),
                       cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(dev_dir_buff, host_dir_buff.data(),
                       sizeof(dna::Direction) * host_barriers.size(), cudaMemcpyHostToDevice));

  CUDA_CALL(cudaDeviceSynchronize());
  setup_curand_kernel<<<grid_size, block_size>>>(dev_rng_states, seed);
  init_barriers<<<grid_size, block_size>>>(dev_barriers, host_barriers.size(), dev_pos_buff, dev_dir_buff);
  CUDA_CALL(cudaDeviceSynchronize());
  init_ctcf_states_kernel<<<grid_size, block_size>>>(dev_barrier_states, host_barriers.size(), dev_rng_states);
  CUDA_CALL(cudaDeviceSynchronize());

  mykernel3<<<grid_size, block_size>>>(dev_barriers, host_barriers.size(), dev_barrier_states, dev_rng_states);
  CUDA_CALL(cudaDeviceSynchronize());
  CUDA_CALL(cudaFree(dev_barriers));
  CUDA_CALL(cudaFree(dev_rng_states));

  std::vector<CTCF::State> host_states(host_barriers.size());
  CUDA_CALL(cudaMemcpy(host_states.data(), dev_barrier_states, host_states.size(),
                       cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaDeviceSynchronize());
  CUDA_CALL(cudaFree(dev_barrier_states));
  CUDA_CALL(cudaDeviceSynchronize());

  return host_states;
}
}  // namespace modle::cu::Simulation
