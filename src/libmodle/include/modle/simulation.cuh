#pragma once
#include <thrust/device_vector.h>

#include <cstdint>
#include <cuda/std/atomic>
#include <vector>

#include "modle/contacts.cuh"
namespace modle::cu::Simulation {
__global__ void mykernel();
__global__ void mykernel2(modle::cu::ContactMatrix<uint32_t>* m, cuda::std::atomic<uint32_t>* buff,
                          size_t nrows, size_t ncols, size_t& missed_updates, size_t& tot_contacts);

[[nodiscard]] std::vector<uint32_t> run_mykernel2(size_t nrows, size_t ncols,
                                                  size_t& missed_updates, size_t& tot_contacts);

}  // namespace modle::cu::Simulation
