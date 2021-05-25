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
  this->_grid_size = new_tasks.size();
  this->_global_state_host.init(this->_grid_size, this->_block_size, new_tasks, barrier_pos,
                                barrier_dir, barrier_probs_occ_to_occ, barrier_probs_nocc_to_nocc);

  kernels::init_barrier_states<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  throw_on_cuda_error(
      cudaGetLastError(),
      fmt::format(
          FMT_STRING(
              "The following error occurred while initializing barrier states for batch #{}"),
          batch_number));

  kernels::compute_initial_loading_epochs<<<this->_grid_size, this->_block_size>>>(
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
    kernels::bind_and_sort_lefs<<<this->_grid_size, this->_block_size>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(
        cudaGetLastError(),
        fmt::format(
            FMT_STRING(
                "The following error occurred while binding LEFs during epoch #{} of batch #{}"),
            this->_current_epoch, batch_number));

    kernels::select_lefs_then_register_contacts<<<this->_grid_size, this->_block_size>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(cudaGetLastError(),
                        fmt::format(FMT_STRING("The following error occurred while sampling "
                                               "contacts during epoch #{} of batch #{}"),
                                    this->_current_epoch, batch_number));

    kernels::generate_moves<<<this->_grid_size, this->_block_size>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(cudaGetLastError(),
                        fmt::format(FMT_STRING("The following error occurred while generating "
                                               "moves during epoch #{} of batch #{}"),
                                    this->_current_epoch, batch_number));

    kernels::update_ctcf_states<<<this->_grid_size, this->_block_size>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(
        cudaGetLastError(),
        fmt::format(FMT_STRING("The following error occurred while updating extrusion barrier "
                               "states during epoch #{} of batch #{}"),
                    this->_current_epoch, batch_number));

    kernels::process_collisions<<<this->_grid_size, this->_block_size>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(cudaGetLastError(),
                        fmt::format(FMT_STRING("The following error occurred while processing LEF "
                                               "collisions during epoch #{} of batch #{}"),
                                    this->_current_epoch, batch_number));

    kernels::extrude_and_release_lefs<<<this->_grid_size, this->_block_size>>>(
        this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(
        cudaGetLastError(),
        fmt::format(
            FMT_STRING(
                "The following error occurred while extruding during epoch #{} of batch #{}"),
            this->_current_epoch, batch_number));

    kernels::advance_epoch<<<1, 1>>>(this->_global_state_host.get_ptr_to_dev_instance());
    throw_on_cuda_error(
        cudaGetLastError(),
        fmt::format(
            FMT_STRING(
                "The following error occurred while advancing epoch from {} to {} for batch #{}"),
            this->_current_epoch, this->_current_epoch + 1, batch_number));

    this->sync();
    const auto state = this->_global_state_host.get_copy_of_device_instance();
    if (state.ntasks_completed == state.ntasks) {
      break;
    }
  }
}
}  // namespace modle::cu
