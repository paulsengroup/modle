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
#include "modle/cu/simulation_internal.hpp"

namespace modle::cu {

void Simulation::run_batch(const std::vector<Task>& new_tasks,
                           const std::vector<uint32_t>& barrier_pos,
                           const std::vector<dna::Direction>& barrier_dir,
                           const std::vector<float>& barrier_probs_occ_to_occ,
                           const std::vector<float>& barrier_probs_nocc_to_nocc) {
  assert(!new_tasks.empty());    // NOLINT
  assert(!barrier_pos.empty());  // NOLINT
  this->_grid_size = new_tasks.size();
  this->_global_state_host.init(this->_grid_size, this->_block_size, new_tasks, barrier_pos,
                                barrier_dir, barrier_probs_occ_to_occ, barrier_probs_nocc_to_nocc);

  kernels::init_barrier_states<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  this->setup_burnin_phase();

  // const auto shared_memory_size = std::max(nlefs, nbarriers) * sizeof(uint32_t);

  this->_current_epoch = 0UL;
  do {
    if (this->_current_epoch++ % 50 == 0) {
      fmt::print(stderr, FMT_STRING("Simulating epoch #{}...\n"), this->_current_epoch - 1);
    }
  } while (this->simulate_next_epoch());

  std::vector<uint32_t> buff(new_tasks[0].nlefs);
  BlockState bstate;
  cuda::memory::copy_single(&bstate,
                            this->_global_state_host.get_copy_of_device_instance().block_states);
  cuda::memory::copy(buff.data(), bstate.rev_unit_pos, buff.size() * sizeof(bp_t));
  fmt::print(stderr, "rev_pos = [{}", buff.front());
  for (auto i = 1UL; i < buff.size(); ++i) {
    fmt::print(stderr, ", {}", buff[i]);
  }
  fmt::print(stderr, "]\n");
  cuda::memory::copy(buff.data(), bstate.fwd_unit_pos, buff.size() * sizeof(bp_t));
  fmt::print(stderr, "fwd_pos = [{}", buff.front());
  for (auto i = 1UL; i < buff.size(); ++i) {
    fmt::print(stderr, ", {}", buff[i]);
  }
  fmt::print(stderr, "]\n");

  std::vector<char> buff1(new_tasks[0].nbarriers);
  cuda::memory::copy(buff1.data(), bstate.barrier_mask, buff1.size() * sizeof(bool));
  fmt::print(stderr, "barrier_states = [{}", static_cast<int>(buff1.front()));
  for (auto i = 1UL; i < buff1.size(); ++i) {
    fmt::print(stderr, ", {}", static_cast<int>(buff1[i]));
  }
  fmt::print(stderr, "]\n");
}

void Simulation::setup_burnin_phase() {
  const auto num_epochs_for_initial_loading =
      std::accumulate(this->_tasks.begin(), this->_tasks.end(), 0UL,
                      [](auto accumulator, const auto& task) { return accumulator + task.nlefs; });

  kernels::generate_initial_loading_epochs<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance(), num_epochs_for_initial_loading);

  auto& offsets_buff_host = this->_global_state_host._buff1;
  auto* offsets_buff_dev = this->_global_state_host.sorting_offset1_buff.get();
  using OffsetT = std::remove_pointer_t<decltype(offsets_buff_dev)>;

  offsets_buff_host.resize(this->_tasks.size() + 1);
  offsets_buff_host.front() = 0;
  std::transform(this->_tasks.begin(), this->_tasks.end(), offsets_buff_host.begin() + 1,
                 [offset = 0U](const auto& task) mutable { return (offset += task.nlefs); });

  cuda::memory::copy(offsets_buff_dev, offsets_buff_host.data(),
                     offsets_buff_host.size() * sizeof(OffsetT));

  auto* begin_offsets_dev = offsets_buff_dev;
  auto* end_offsets_dev = offsets_buff_dev + 1;

  cub::DoubleBuffer<uint32_t> keys_dev(this->_global_state_host.large_uint_buff1.get(),
                                       this->_global_state_host.large_uint_buff2.get());

  const auto num_items_to_sort =
      std::accumulate(this->_tasks.begin(), this->_tasks.end(), 0,
                      [](auto accumulator, const auto& task) { return accumulator + task.nlefs; });
  const auto num_segments_to_sort = static_cast<int>(this->_grid_size);

  auto tmp_storage_bytes_required = 0UL;
  // This call just computes the required temp buffer size
  auto status = cub::DeviceSegmentedRadixSort::SortKeys(
      nullptr, tmp_storage_bytes_required, keys_dev, num_items_to_sort, num_segments_to_sort,
      begin_offsets_dev, end_offsets_dev);

  if (status != cudaSuccess) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to compute the memory requirements to sort LEF initial "
                               "binding epochs on the device ({}:{}): {}: {}"),
                    __FILE__, __LINE__, cudaGetErrorName(status), cudaGetErrorString(status)));
  }

  if (tmp_storage_bytes_required > this->_global_state_host.tmp_sorting_storage_bytes) {
    this->_global_state_host.tmp_sorting_storage_bytes = tmp_storage_bytes_required;
    this->_global_state_host.tmp_sorting_storage = cuda::memory::device::make_unique<uint32_t[]>(
        this->_global_state_host._device, this->_global_state_host.tmp_sorting_storage_bytes);
    this->_global_state_host.sync_state_with_device();
  }

  this->_global_state_host._device.synchronize();

  status = cub::DeviceSegmentedRadixSort::SortKeys(
      this->_global_state_host.tmp_sorting_storage.get(),
      this->_global_state_host.tmp_sorting_storage_bytes, keys_dev, num_items_to_sort,
      num_segments_to_sort, begin_offsets_dev, end_offsets_dev);
  this->_global_state_host._device.synchronize();

  if (status != cudaSuccess) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to sort LEF initial binding epochs on the device ({}:{}): {}: {}"),
        __FILE__, __LINE__, cudaGetErrorName(status), cudaGetErrorString(status)));
  }

  try {
    auto global_state_host = this->_global_state_host.get_copy_of_device_instance();
    kernels::shift_and_scatter_lef_loading_epochs<<<this->_grid_size, this->_block_size>>>(
        this->_global_state_host.large_uint_buff1.get(), global_state_host.block_states,
        begin_offsets_dev, end_offsets_dev);
    this->_global_state_host._device.synchronize();
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to shift and scatter LEF initial binding epochs: {}"), e.what()));
  }
}

void Simulation::sort_lefs() {
  auto* begin_offsets_dev = this->_global_state_host.sorting_offset1_buff.get();
  auto* end_offsets_dev = this->_global_state_host.sorting_offset2_buff.get();

  const auto num_segments_to_sort = this->_global_state_host.ntasks;

  auto tmp_storage_bytes_required = 0UL;
  auto* tmp_storage_dev = this->_global_state_host.tmp_sorting_storage.get();

  cub::DoubleBuffer<uint32_t> keys_dev(this->_global_state_host.large_uint_buff1.get(),
                                       this->_global_state_host.large_uint_buff2.get());
  cub::DoubleBuffer<uint32_t> vals_dev(this->_global_state_host.large_uint_buff3.get(),
                                       this->_global_state_host.large_uint_buff4.get());

  auto num_items_to_sort_dev =
      cuda::memory::device::allocate(this->_global_state_host._device, sizeof(uint32_t));
  cuda::memory::device::zero(num_items_to_sort_dev);

  auto direction = dna::Direction::rev;

  kernels::prepare_extr_units_for_sorting<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance(), direction,
      static_cast<uint32_t*>(num_items_to_sort_dev.get()));

  uint32_t num_items_to_sort_host;
  cuda::memory::copy_single(&num_items_to_sort_host,
                            static_cast<uint32_t*>(num_items_to_sort_dev.get()));
  cuda::memory::device::free(num_items_to_sort_dev);

  auto status = cub::DeviceSegmentedRadixSort::SortPairs(
      nullptr, tmp_storage_bytes_required, keys_dev, vals_dev,
      static_cast<int>(num_items_to_sort_host), static_cast<int>(num_segments_to_sort),
      begin_offsets_dev, end_offsets_dev);
  if (status != cudaSuccess) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "Failed to compute the memory requirements to sort LEFs on the device ({}:{}): {}: {}"),
        __FILE__, __LINE__, cudaGetErrorName(status), cudaGetErrorString(status)));
  }

  if (tmp_storage_bytes_required > this->_global_state_host.tmp_sorting_storage_bytes) {
    this->_global_state_host.tmp_sorting_storage_bytes = tmp_storage_bytes_required;
    this->_global_state_host.tmp_sorting_storage = cuda::memory::device::make_unique<uint32_t[]>(
        this->_global_state_host._device, tmp_storage_bytes_required);
    tmp_storage_dev = this->_global_state_host.tmp_sorting_storage.get();
    this->_global_state_host.sync_state_with_device();
  }

  status = cub::DeviceSegmentedRadixSort::SortPairs(
      tmp_storage_dev, tmp_storage_bytes_required, keys_dev, vals_dev,
      static_cast<int>(num_items_to_sort_host), static_cast<int>(num_segments_to_sort),
      begin_offsets_dev, end_offsets_dev);
  if (status != cudaSuccess) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to sort LEFs on the device ({}:{}): {}: {}"), __FILE__,
                    __LINE__, cudaGetErrorName(status), cudaGetErrorString(status)));
  }

  this->_global_state_host._device.synchronize();
  kernels::update_unit_mappings_and_scatter_sorted_lefs<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance(), direction);

  direction = dna::Direction::fwd;

  kernels::prepare_extr_units_for_sorting<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance(), direction);

  status = cub::DeviceSegmentedRadixSort::SortPairs(
      tmp_storage_dev, tmp_storage_bytes_required, keys_dev, vals_dev,
      static_cast<int>(num_items_to_sort_host), static_cast<int>(num_segments_to_sort),
      begin_offsets_dev, end_offsets_dev);
  if (status != cudaSuccess) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to sort LEFs on the device ({}:{}): {}: {}"), __FILE__,
                    __LINE__, cudaGetErrorName(status), cudaGetErrorString(status)));
  }

  this->_global_state_host._device.synchronize();

  kernels::update_unit_mappings_and_scatter_sorted_lefs<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance(), direction);
  this->_global_state_host._device.synchronize();
}

void Simulation::select_lefs_and_register_contacts() {
  auto num_items_to_sort_dev =
      cuda::memory::device::allocate(this->_global_state_host._device, sizeof(uint32_t));
  cuda::memory::device::zero(num_items_to_sort_dev);

  kernels::prepare_units_for_random_shuffling<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance(),
      static_cast<uint32_t*>(num_items_to_sort_dev.get()));
  this->_global_state_host._device.synchronize();
  uint32_t num_items_to_sort_host;
  cuda::memory::copy_single(&num_items_to_sort_host,
                            static_cast<uint32_t*>(num_items_to_sort_dev.get()));
  cuda::memory::device::free(num_items_to_sort_dev);

  auto* begin_offsets_dev = this->_global_state_host.sorting_offset1_buff.get();
  auto* end_offsets_dev = this->_global_state_host.sorting_offset2_buff.get();

  const auto num_segments_to_sort = this->_global_state_host.ntasks;

  auto tmp_storage_bytes_required = 0UL;
  auto* tmp_storage_dev = this->_global_state_host.tmp_sorting_storage.get();

  cub::DoubleBuffer<uint32_t> keys_dev(this->_global_state_host.large_uint_buff3.get(),
                                       this->_global_state_host.large_uint_buff4.get());
  cub::DoubleBuffer<uint32_t> vals_dev(this->_global_state_host.large_uint_buff1.get(),
                                       this->_global_state_host.large_uint_buff2.get());

  auto status = cub::DeviceSegmentedRadixSort::SortPairs(
      nullptr, tmp_storage_bytes_required, keys_dev, vals_dev,
      static_cast<int>(num_items_to_sort_host), static_cast<int>(num_segments_to_sort),
      begin_offsets_dev, end_offsets_dev);
  this->_global_state_host._device.synchronize();

  if (status != cudaSuccess) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to compute the memory requirements to sort LEF indexes needed to "
                   "perform random shuffling ({}:{}): {}: {}"),
        __FILE__, __LINE__, cudaGetErrorName(status), cudaGetErrorString(status)));
  }

  if (tmp_storage_bytes_required > this->_global_state_host.tmp_sorting_storage_bytes) {
    this->_global_state_host.tmp_sorting_storage_bytes = tmp_storage_bytes_required;
    this->_global_state_host.tmp_sorting_storage = cuda::memory::device::make_unique<uint32_t[]>(
        this->_global_state_host._device, tmp_storage_bytes_required);
    tmp_storage_dev = this->_global_state_host.tmp_sorting_storage.get();
    this->_global_state_host.sync_state_with_device();
  }

  status = cub::DeviceSegmentedRadixSort::SortPairs(
      tmp_storage_dev, tmp_storage_bytes_required, keys_dev, vals_dev,
      static_cast<int>(num_items_to_sort_host), static_cast<int>(num_segments_to_sort),
      begin_offsets_dev, end_offsets_dev);
  this->_global_state_host._device.synchronize();

  if (status != cudaSuccess) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to sort LEF indexes needed to perform a random shuffling on "
                               "the device ({}:{}): {}: {}"),
                    __FILE__, __LINE__, cudaGetErrorName(status), cudaGetErrorString(status)));
  }

  kernels::select_lefs_then_register_contacts<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());

  this->_global_state_host._device.synchronize();
}

void Simulation::generate_moves() {
  kernels::generate_moves<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
}

void Simulation::update_ctcf_states() {
  kernels::update_ctcf_states<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
}

void Simulation::process_collisions() {
  kernels::process_collisions<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
}

void Simulation::extrude_and_release_lefs() {
  kernels::extrude_and_release_lefs<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
}

bool Simulation::simulate_next_epoch() {
  kernels::select_and_bind_lefs<<<this->_grid_size, this->_block_size>>>(
      this->_current_epoch, this->_global_state_host.get_ptr_to_dev_instance());

  this->sort_lefs();
  this->_global_state_host._device.synchronize();
  this->select_lefs_and_register_contacts();
  this->_global_state_host._device.synchronize();
  this->generate_moves();
  this->_global_state_host._device.synchronize();
  this->update_ctcf_states();
  this->_global_state_host._device.synchronize();
  kernels::reset_collision_masks<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host.get_ptr_to_dev_instance());
  this->_global_state_host._device.synchronize();
  this->process_collisions();
  this->_global_state_host._device.synchronize();
  this->extrude_and_release_lefs();
  this->_global_state_host._device.synchronize();

  const auto state = this->_global_state_host.get_copy_of_device_instance();
  return state.ntasks_completed != state.ntasks;  // End of simulation
}

}  // namespace modle::cu
