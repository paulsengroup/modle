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
  this->setup_burnin_phase();
  this->_global_state_host._device.synchronize();

  this->_current_epoch = 0UL;
  while (this->simulate_next_epoch()) {
    if (this->_current_epoch++ % 50 == 0) {
      fmt::print(stderr,
                 FMT_STRING("Simulating epoch #{} of batch #{}. grid_size={}; block_size={};\n"),
                 this->_current_epoch - 1, batch_number, this->_grid_size, this->_block_size);
    }
  }
}

void Simulation::setup_burnin_phase() {
  kernels::compute_initial_loading_epochs<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  this->_global_state_host._device.synchronize();
}

void Simulation::bind_and_sort_lefs() {
  kernels::bind_and_sort_lefs<<<this->_grid_size, this->_block_size>>>(
      this->_current_epoch, this->_global_state_host.get_ptr_to_dev_instance());
  this->_global_state_host._device.synchronize();
}

void Simulation::select_lefs_and_register_contacts() {
  kernels::select_lefs_then_register_contacts<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  this->_global_state_host._device.synchronize();
}

void Simulation::generate_moves() {
  kernels::generate_moves<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  this->_global_state_host._device.synchronize();
}

void Simulation::update_ctcf_states() {
  kernels::update_ctcf_states<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  this->_global_state_host._device.synchronize();
}

void Simulation::process_collisions() {
  kernels::process_collisions<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  this->_global_state_host._device.synchronize();
}

void Simulation::extrude_and_release_lefs() {
  kernels::extrude_and_release_lefs<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  this->_global_state_host._device.synchronize();
}

bool Simulation::simulate_next_epoch() {
  this->bind_and_sort_lefs();

  this->select_lefs_and_register_contacts();

  this->generate_moves();

  this->update_ctcf_states();

  this->process_collisions();

  this->extrude_and_release_lefs();

  const auto state = this->_global_state_host.get_copy_of_device_instance();
  return state.ntasks_completed != state.ntasks;  // End of simulation
}

}  // namespace modle::cu
