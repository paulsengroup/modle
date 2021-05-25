#include <absl/hash/hash.h>
#include <absl/time/clock.h>
#include <curand_kernel.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <thrust/sort.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cub/device/device_segmented_radix_sort.cuh>
#include <cuda/runtime_api.hpp>
#include <numeric>
#include <string>
#include <vector>

#include "modle/cu/common.hpp"
#include "modle/cu/config.hpp"
#include "modle/cu/simulation.hpp"
#include "modle/cu/simulation_internal.cuh"
namespace modle::cu {

void throw_on_cuda_error(const cudaError status, std::string_view exception_prefix_message) {
  if (status != cudaSuccess) {
    if (exception_prefix_message.empty()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("The device encountered the following error: {}: {}"),
                      cudaGetErrorName(status), cudaGetErrorString(status)));
    }
    throw std::runtime_error(fmt::format(FMT_STRING("{}: {}: {}"), exception_prefix_message,
                                         cudaGetErrorName(status), cudaGetErrorString(status)));
  }
}

void Simulation::sync() { this->_global_state_host._device.synchronize(); }

void Simulation::run_batch(const std::vector<Task>& new_tasks,
                           const std::vector<uint32_t>& barrier_pos,
                           const std::vector<dna::Direction>& barrier_dir,
                           const std::vector<float>& barrier_probs_occ_to_occ,
                           const std::vector<float>& barrier_probs_nocc_to_nocc,
                           size_t batch_number) {
  assert(!new_tasks.empty());    // NOLINT
  assert(!barrier_pos.empty());  // NOLINT
  const auto t0 = absl::Now();

  this->_global_state_host.shared_memory_bytes_per_block = this->compute_shared_memory_per_block();
  this->_grid_size = new_tasks.size();
  this->_global_state_host.init(this->_grid_size, this->_block_size, new_tasks, barrier_pos,
                                barrier_dir, barrier_probs_occ_to_occ, barrier_probs_nocc_to_nocc);

  const auto& gs = this->_grid_size;
  const auto& bs = this->_block_size;
  const auto& smem = this->_global_state_host.shared_memory_bytes_per_block;
  cudaStream_t stream;
  cudaStreamCreate(&stream);

  kernels::init_barrier_states<<<gs, bs, 0, stream>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  throw_on_cuda_error(
      cudaGetLastError(),
      fmt::format(
          FMT_STRING(
              "The following error occurred while initializing barrier states for batch #{}"),
          batch_number));

  kernels::compute_initial_loading_epochs<<<this->_grid_size, this->_block_size, 0, stream>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  throw_on_cuda_error(
      cudaGetLastError(),
      fmt::format(
          FMT_STRING(
              "The following error occurred while computing initial loading epochs for batch #{}"),
          batch_number));
  this->sync();

  for (this->_current_epoch = 0; true; ++this->_current_epoch) {
    if (this->_current_epoch % 50 == 0) {
      fmt::print(stderr,
                 FMT_STRING("Simulating epoch #{} of batch #{}. grid_size={}; block_size={};\n"),
                 this->_current_epoch, batch_number, this->_grid_size, this->_block_size);
    }
    kernels::bind_and_sort_lefs<<<gs, bs, smem, stream>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(
        cudaGetLastError(),
        fmt::format(
            FMT_STRING(
                "The following error occurred while binding LEFs during epoch #{} of batch #{}"),
            this->_current_epoch, batch_number));

    kernels::select_lefs_then_register_contacts<<<gs, bs, 0, stream>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(cudaGetLastError(),
                        fmt::format(FMT_STRING("The following error occurred while sampling "
                                               "contacts during epoch #{} of batch #{}"),
                                    this->_current_epoch, batch_number));

    kernels::generate_moves<<<gs, bs, 0, stream>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(cudaGetLastError(),
                        fmt::format(FMT_STRING("The following error occurred while generating "
                                               "moves during epoch #{} of batch #{}"),
                                    this->_current_epoch, batch_number));

    kernels::update_ctcf_states<<<gs, bs, 0, stream>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(
        cudaGetLastError(),
        fmt::format(FMT_STRING("The following error occurred while updating extrusion barrier "
                               "states during epoch #{} of batch #{}"),
                    this->_current_epoch, batch_number));

    kernels::process_collisions<<<gs, bs, 0, stream>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(cudaGetLastError(),
                        fmt::format(FMT_STRING("The following error occurred while processing LEF "
                                               "collisions during epoch #{} of batch #{}"),
                                    this->_current_epoch, batch_number));

    kernels::extrude_and_release_lefs<<<gs, bs, 0, stream>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(
        cudaGetLastError(),
        fmt::format(
            FMT_STRING(
                "The following error occurred while extruding during epoch #{} of batch #{}"),
            this->_current_epoch, batch_number));

    kernels::advance_epoch<<<1, 1, 0, stream>>>(this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(
        cudaGetLastError(),
        fmt::format(
            FMT_STRING(
                "The following error occurred while advancing epoch from {} to {} for batch #{}"),
            this->_current_epoch, this->_current_epoch + 1, batch_number));

    throw_on_cuda_error(
        cudaStreamSynchronize(stream),
        fmt::format(FMT_STRING("An error occurred while simulating epoch #{} for batch #{}"),
                    this->_current_epoch, batch_number));
    this->sync();
    const auto state = this->_global_state_host.get_copy_of_device_instance();
    if (state.ntasks_completed == state.ntasks) {
      throw_on_cuda_error(cudaStreamDestroy(stream), "Failed to destroy CUDA stream");
      break;
    }
  }
}

size_t Simulation::compute_shared_memory_per_block() const {
  const auto prop = this->_global_state_host._device.properties();
  //   const auto num_sm = prop.multiProcessorCount;
  const auto max_num_blocks_per_sm = static_cast<size_t>(prop.maxBlocksPerMultiProcessor);
  const auto max_num_threads_per_sm = static_cast<size_t>(prop.maxThreadsPerMultiProcessor);
  const auto shared_mem_per_sm = prop.sharedMemPerMultiprocessor;

  const auto num_resident_blocks_per_sm =
      std::min(max_num_threads_per_sm / this->_block_size, max_num_blocks_per_sm);
  const auto effective_shared_mem_per_block = shared_mem_per_sm / num_resident_blocks_per_sm;

  return effective_shared_mem_per_block;
}

}  // namespace modle::cu
