// clang-format off
//#include <cuda_runtime_api.h>
#include "modle/cu/simulation.hpp"
#include <curand_kernel.h>
#include <cub/cub.cuh>
#include <cuda/runtime_api.hpp>
// clang-format on

#include <absl/hash/hash.h>
#include <absl/time/clock.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

#include "modle/config_cuda.hpp"
#include "modle/cu/common.hpp"
#include "modle/cu/simulation_internal.hpp"

namespace modle::cu {
/*
void Simulation::init_global_state(size_t grid_size, size_t block_size, size_t max_grid_size,
                                   size_t max_nlefs, size_t max_nbarriers, cuda::device_t dev) {
  this->_grid_size = grid_size;
  this->_block_size = block_size;
  this->_global_state_host = std::make_unique<GlobalStateHost>(grid_size, block_size, max_grid_size,
                                                               this->_config->to_cuda_config(),
                                                               max_nlefs, max_nbarriers, dev);
}
 */

void Simulation::update_tasks(const std::vector<Task>& new_tasks) noexcept {
  this->_tasks = new_tasks;
  this->_grid_size = this->_tasks.size();

  this->_global_state_host->write_tasks_to_device(this->_tasks);
}

void Simulation::update_barriers(const std::vector<uint32_t>& barrier_pos,
                                 const std::vector<dna::Direction>& barrier_dir) noexcept {
  assert(barrier_pos.size() == barrier_dir.size());  // NOLINT
  this->_barrier_positions = barrier_pos;
  this->_barrier_directions = barrier_dir;

  this->_global_state_host->write_barriers_to_device(this->_barrier_positions,
                                                     this->_barrier_directions);
}

void Simulation::run_batch(const std::vector<Task>& new_tasks,
                           const std::vector<uint32_t>& barrier_pos,
                           std::vector<dna::Direction>& barrier_dir) {
  this->update_tasks(new_tasks);
  this->update_barriers(barrier_pos, barrier_dir);
  assert(this->_grid_size == this->_tasks.size() && this->_grid_size != 0);  // NOLINT
  // TODO handle case where nlefs and nbarriers can be different across states

  this->_global_state_host->init(this->_grid_size, this->_block_size);
  this->setup_burnin_phase();

  // const auto shared_memory_size = std::max(nlefs, nbarriers) * sizeof(uint32_t);

  bool end_of_simulation{false};

  for (this->_current_epoch = 0UL; !end_of_simulation; ++this->_current_epoch) {
    end_of_simulation = this->simulate_one_epoch();
    if (this->_current_epoch == 10) {
      return;
    }
  }
}

void Simulation::setup_burnin_phase() {
  const auto num_epochs_for_initial_loading =
      std::accumulate(this->_tasks.begin(), this->_tasks.end(), 0UL,
                      [](auto accumulator, const auto& task) { return accumulator + task.nlefs; });

  kernels::generate_initial_loading_epochs<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host->get_ptr_to_dev_instance(), num_epochs_for_initial_loading);

  this->_global_state_host->_device.synchronize();

  auto& offsets_buff_host = this->_global_state_host->_buff1;
  auto* offsets_buff_dev = this->_global_state_host->sorting_offset1_buff.get();
  using OffsetT = std::remove_pointer_t<decltype(offsets_buff_dev)>;

  offsets_buff_host.resize(this->_tasks.size() + 1);
  offsets_buff_host.front() = 0;
  std::transform(this->_tasks.begin(), this->_tasks.end(), offsets_buff_host.begin() + 1,
                 [offset = 0U](const auto& task) mutable { return (offset += task.nlefs); });

  cuda::memory::copy(offsets_buff_dev, offsets_buff_dev,
                     offsets_buff_host.size() * sizeof(OffsetT));

  auto* begin_offsets_dev = offsets_buff_dev;
  auto* end_offsets_dev = offsets_buff_dev + 1;

  cub::DoubleBuffer<uint32_t> keys_dev(this->_global_state_host->large_uint_buff1.get(),
                                       this->_global_state_host->large_uint_buff2.get());

  const auto num_items_to_sort =
      std::accumulate(this->_tasks.begin(), this->_tasks.end(), 0,
                      [](auto accumulator, const auto& task) { return accumulator + task.nlefs; });
  const auto num_segments_to_sort = static_cast<int>(this->_grid_size);

  auto tmp_storage_bytes_required = 0UL;
  // This call just computes the required temp buffer size
  CUDA_CALL(cub::DeviceSegmentedRadixSort::SortKeys(nullptr, tmp_storage_bytes_required, keys_dev,
                                                    num_items_to_sort, num_segments_to_sort,
                                                    begin_offsets_dev, end_offsets_dev));

  if (tmp_storage_bytes_required > this->_global_state_host->tmp_sorting_storage_bytes) {
    this->_global_state_host->tmp_sorting_storage = cuda::memory::device::make_unique<uint32_t[]>(
        this->_global_state_host->_device,
        (tmp_storage_bytes_required + sizeof(uint32_t) - 1) / sizeof(uint32_t));
    this->_global_state_host->tmp_sorting_storage_bytes = tmp_storage_bytes_required;
    this->_global_state_host->sync_state_with_device();
  }

  const auto status1 = cub::DeviceSegmentedRadixSort::SortKeys(
      this->_global_state_host->tmp_sorting_storage.get(),
      this->_global_state_host->tmp_sorting_storage_bytes, keys_dev, num_items_to_sort,
      num_segments_to_sort, begin_offsets_dev, end_offsets_dev);
  this->_global_state_host->_device.synchronize();

  if (status1 != cudaSuccess) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to sort LEF initial binding epochs ({}:{}): {}: {}"),
                    __FILE__, __LINE__, cudaGetErrorName(status1), cudaGetErrorString(status1)));
  }

  try {
    auto global_state_host = this->_global_state_host->get_copy_of_device_instance();
    kernels::shift_and_scatter_lef_loading_epochs<<<this->_grid_size, this->_block_size>>>(
        this->_global_state_host->large_uint_buff1.get(), global_state_host.block_states,
        begin_offsets_dev, end_offsets_dev);
    this->_global_state_host->_device.synchronize();
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to shift and scatter LEF initial binding epochs: {}"), e.what()));
  }
}

void Simulation::sort_lefs() {
  auto& offsets_buff_host = this->_global_state_host->_buff1;
  auto* begin_offsets_dev = this->_global_state_host->sorting_offset1_buff.get();
  auto* end_offsets_dev = this->_global_state_host->sorting_offset2_buff.get();

  const auto num_segments_to_sort = this->_global_state_host->ntasks;

  auto tmp_storage_bytes_required = 0UL;
  auto* tmp_storage_dev = this->_global_state_host->tmp_sorting_storage.get();

  cub::DoubleBuffer<uint32_t> keys_dev(this->_global_state_host->large_uint_buff1.get(),
                                       this->_global_state_host->large_uint_buff2.get());
  cub::DoubleBuffer<uint32_t> vals_dev(this->_global_state_host->large_uint_buff3.get(),
                                       this->_global_state_host->large_uint_buff4.get());

  auto num_items_to_sort_dev =
      cuda::memory::device::allocate(this->_global_state_host->_device, sizeof(uint32_t));
  cuda::memory::device::zero(num_items_to_sort_dev);

  auto direction = dna::Direction::rev;

  kernels::prepare_extr_units_for_sorting<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host->get_ptr_to_dev_instance(), direction,
      static_cast<uint32_t*>(num_items_to_sort_dev.get()));

  uint32_t num_items_to_sort_host;
  cuda::memory::copy(&num_items_to_sort_host, num_items_to_sort_dev.get(),
                     num_items_to_sort_dev.size());
  cuda::memory::device::free(num_items_to_sort_dev);

  CUDA_CALL(cub::DeviceSegmentedRadixSort::SortPairs(
      nullptr, tmp_storage_bytes_required, keys_dev, vals_dev,
      static_cast<int>(num_items_to_sort_host), static_cast<int>(num_segments_to_sort),
      begin_offsets_dev, end_offsets_dev));

  if (tmp_storage_bytes_required > this->_global_state_host->tmp_sorting_storage_bytes) {
    this->_global_state_host->tmp_sorting_storage = cuda::memory::device::make_unique<uint32_t[]>(
        this->_global_state_host->_device, tmp_storage_bytes_required);
    tmp_storage_dev = this->_global_state_host->tmp_sorting_storage.get();
    this->_global_state_host->tmp_sorting_storage_bytes = tmp_storage_bytes_required;
    this->_global_state_host->sync_state_with_device();
  }

  CUDA_CALL(cub::DeviceSegmentedRadixSort::SortPairs(
      tmp_storage_dev, tmp_storage_bytes_required, keys_dev, vals_dev,
      static_cast<int>(num_items_to_sort_host), static_cast<int>(num_segments_to_sort),
      begin_offsets_dev, end_offsets_dev));

  this->_global_state_host->_device.synchronize();

  direction = dna::Direction::fwd;

  kernels::prepare_extr_units_for_sorting<<<this->_grid_size, this->_block_size>>>(
      this->_global_state_host->get_ptr_to_dev_instance(), direction,
      static_cast<uint32_t*>(num_items_to_sort_dev.get()));

  CUDA_CALL(cub::DeviceSegmentedRadixSort::SortPairs(
      tmp_storage_dev, tmp_storage_bytes_required, keys_dev, vals_dev,
      static_cast<int>(num_items_to_sort_host), static_cast<int>(num_segments_to_sort),
      begin_offsets_dev, end_offsets_dev));

  this->_global_state_host->_device.synchronize();

  kernels::scatter_sorted_lefs<<<this->_block_size, this->_grid_size>>>(
      this->_global_state_host->get_ptr_to_dev_instance(), direction,
      static_cast<uint32_t*>(num_items_to_sort_dev.get()));
  this->_global_state_host->_device.synchronize();

  /*

  GlobalStateDev global_state_host;
  CUDA_CALL(cudaMemcpy(&global_state_host, global_state_dev, sizeof(GlobalStateDev),
                       cudaMemcpyDeviceToHost));

  auto* d_begin_offsets = global_state_host.sorting_offset1_buff;
  auto* d_end_offsets = global_state_host.sorting_offset2_buff;

  const auto num_segments_to_sort = global_state_host.ntasks;

  auto temp_storage_bytes = global_state_host.tmp_sorting_storage_bytes;
  auto* d_tmp_storage = global_state_host.tmp_sorting_storage;

  cub::DoubleBuffer<uint32_t> d_keys(global_state_host.large_uint_buff1,
                                     global_state_host.large_uint_buff2);
  cub::DoubleBuffer<uint32_t> d_vals(global_state_host.large_uint_buff3,
                                     global_state_host.large_uint_buff4);

  uint32_t* num_items_to_sort_dev{nullptr};
  uint32_t num_items_to_sort_host{0};
  CUDA_CALL(cudaMalloc(&num_items_to_sort_dev, sizeof(uint32_t)));
  CUDA_CALL(cudaMemcpy(num_items_to_sort_dev, &num_items_to_sort_host, sizeof(uint32_t),
                       cudaMemcpyHostToDevice));

  prepare_extr_units_for_sorting<<<grid_size, block_size>>>(global_state_dev, dna::Direction::rev,
                                                            num_items_to_sort_dev);
  CUDA_CALL(cudaMemcpy(&num_items_to_sort_host, num_items_to_sort_dev, sizeof(uint32_t),
                       cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaFree(num_items_to_sort_dev));

  CUDA_CALL(cub::DeviceSegmentedRadixSort::SortPairs(
      nullptr, temp_storage_bytes, d_keys, d_vals, static_cast<int>(num_items_to_sort_host),
      static_cast<int>(num_segments_to_sort), d_begin_offsets, d_end_offsets));
  if (temp_storage_bytes > global_state_host.tmp_sorting_storage_bytes) {
    CUDA_CALL(cudaFree(d_tmp_storage));
    CUDA_CALL(cudaMalloc(&d_tmp_storage, temp_storage_bytes));
    global_state_host.tmp_sorting_storage_bytes = temp_storage_bytes;
    global_state_host.tmp_sorting_storage = d_tmp_storage;
    CUDA_CALL(cudaMemcpy(global_state_dev, &global_state_host, sizeof(GlobalStateDev),
                         cudaMemcpyHostToDevice));
  }

  CUDA_CALL(cub::DeviceSegmentedRadixSort::SortPairs(
      d_tmp_storage, temp_storage_bytes, d_keys, d_vals, static_cast<int>(num_items_to_sort_host),
      static_cast<int>(num_segments_to_sort), d_begin_offsets, d_end_offsets));

  cudaDeviceSynchronize();
  scatter_sorted_lefs<<<block_size, grid_size>>>(global_state_dev, dna::Direction::rev);

  prepare_extr_units_for_sorting<<<grid_size, block_size>>>(global_state_dev, dna::Direction::fwd,
                                                            num_items_to_sort_dev);

  CUDA_CALL(cub::DeviceSegmentedRadixSort::SortPairs(
      d_tmp_storage, temp_storage_bytes, d_keys, d_vals, static_cast<int>(num_items_to_sort_host),
      static_cast<int>(num_segments_to_sort), d_begin_offsets, d_end_offsets));
  cudaDeviceSynchronize();

  scatter_sorted_lefs<<<block_size, grid_size>>>(global_state_dev, dna::Direction::fwd);
   */
}

bool Simulation::simulate_one_epoch() {
  kernels::select_and_bind_lefs<<<this->_grid_size, this->_block_size>>>(
      this->_current_epoch, this->_global_state_host->get_ptr_to_dev_instance());

  this->sort_lefs();

  std::vector<uint32_t> buff1(10000);
  std::vector<uint32_t> buff2(10000);
  GlobalStateDev global_state_host;
  cudaMemcpy(&global_state_host, this->_global_state_host.get(), sizeof(GlobalStateDev),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(buff1.data(), global_state_host.large_uint_buff1, buff1.size() * sizeof(uint32_t),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(buff1.data(), global_state_host.large_uint_buff3, buff2.size() * sizeof(uint32_t),
             cudaMemcpyDeviceToHost);

  printf("pos = [%d", buff1[0]);
  for (auto i = 1UL; i < buff1.size(); ++i) {
    printf(", %d", buff1[i]);
  }

  printf("]\nidx = [%d", buff2[0]);
  for (auto i = 1UL; i < buff2.size(); ++i) {
    printf(", %d", buff2[i]);
  }
  printf("]\n");
  return false;
}

}  // namespace modle::cu
