#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include <cstdint>
#include <cstdio>
#include <vector>

#include "modle/contacts.cuh"
#include "modle/simulation.cuh"

namespace modle::cu::Simulation {

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

std::vector<uint32_t> run_mykernel2(size_t nrows, size_t ncols, size_t& missed_updates,
                                    size_t& tot_contacts) {
  cuda::std::atomic<uint32_t>* dev_buff{nullptr};
  modle::cu::ContactMatrix<uint32_t>* matrix{nullptr};

  if (const auto status = cudaMalloc(&matrix, sizeof(ContactMatrix<uint32_t>));
      status != cudaSuccess) {
    throw std::runtime_error("Unable to allocate enough memory on the device.");
  }
  if (const auto status =
          cudaMalloc(&dev_buff, nrows * ncols * sizeof(cuda::std::atomic<uint32_t>));
      status != cudaSuccess) {
    throw std::runtime_error("Unable to allocate enough memory on the device.");
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
}  // namespace modle::cu::Simulation
