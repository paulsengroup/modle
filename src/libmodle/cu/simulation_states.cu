// clang-format off
#include <cuda/runtime_api.hpp>
//#include <cuda_runtime_api.h>

#include <curand_kernel.h>
#include <cub/cub.cuh>
// clang-format on

#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

#include "modle/config_cuda.hpp"
#include "modle/cu/common.hpp"
#include "modle/cu/simulation.hpp"
#include "modle/cu/simulation_internal.hpp"

namespace modle::cu {

__constant__ struct Config config;  // NOLINT

size_t GlobalStateHost::get_grid_size() const { return this->grid_size; };
size_t GlobalStateHost::get_block_size() const { return this->block_size; };

GlobalStateDev* GlobalStateHost::get_ptr_to_dev_instance() {
  if (!this->_global_state_dev) {
    this->_global_state_dev = static_cast<GlobalStateDev*>(
        cuda::memory::device::allocate(this->_device, sizeof(GlobalStateDev)).get());
    GlobalStateDev gs{};

    gs.config = this->config;

    gs.grid_size = this->grid_size;
    gs.block_size = this->block_size;

    gs.block_states = this->block_states;
    gs.tasks = this->tasks;

    gs.nblock_states = this->nblock_states;
    gs.ntasks = this->ntasks;

    gs.barrier_pos = this->barrier_pos.get();
    gs.barrier_directions = this->barrier_directions.get();

    gs.large_uint_buff1 = this->large_uint_buff1.get();
    gs.large_uint_buff2 = this->large_uint_buff2.get();
    gs.large_uint_buff3 = this->large_uint_buff3.get();
    gs.large_uint_buff4 = this->large_uint_buff4.get();

    gs.large_uint_buff_chunk_alignment = this->large_uint_buff_chunk_alignment;

    gs.sorting_offset1_buff = this->sorting_offset1_buff.get();
    gs.sorting_offset2_buff = this->sorting_offset2_buff.get();

    gs.tmp_sorting_storage = this->tmp_sorting_storage.get();
    gs.tmp_sorting_storage_bytes = this->tmp_sorting_storage_bytes;

    cuda::memory::copy(this->_global_state_dev, &gs, sizeof(GlobalStateDev));
  }
  return this->_global_state_dev;
}

GlobalStateDev GlobalStateHost::get_copy_of_device_instance() {
  GlobalStateDev state_host;
  cuda::memory::copy_single(&state_host, this->get_ptr_to_dev_instance());

  return state_host;
}

GlobalStateHost::GlobalStateHost(size_t grid_size_, size_t block_size_, size_t max_grid_size_,
                                 size_t device_heap_size, const Config& c, size_t max_nlefs,
                                 size_t max_nbarriers, cuda::device_t dev)
    : config(write_config_to_device(c)),
      grid_size(grid_size_),
      block_size(block_size_),
      max_grid_size(max_grid_size_),
      block_states(static_cast<BlockState*>(
          cuda::memory::device::allocate(dev, max_grid_size_ * sizeof(BlockState)).get())),
      tasks(static_cast<Task*>(
          cuda::memory::device::allocate(dev, max_grid_size_ * sizeof(Task)).get())),
      barrier_pos(cuda::memory::device::make_unique<bp_t[]>(dev, max_nbarriers * sizeof(bp_t))),
      barrier_directions(cuda::memory::device::make_unique<dna::Direction[]>(
          dev, max_nbarriers * sizeof(dna::Direction))),
      nblock_states(static_cast<uint32_t>(max_grid_size_)),
      ntasks(static_cast<uint32_t>(max_grid_size_)),
      large_uint_buff1(cuda::memory::device::make_unique<uint32_t[]>(
          dev, max_grid_size_ * std::max(max_nlefs, max_nbarriers))),
      large_uint_buff2(cuda::memory::device::make_unique<uint32_t[]>(
          dev, max_grid_size_ * std::max(max_nlefs, max_nbarriers))),
      large_uint_buff3(cuda::memory::device::make_unique<uint32_t[]>(
          dev, max_grid_size_ * std::max(max_nlefs, max_nbarriers))),
      large_uint_buff4(cuda::memory::device::make_unique<uint32_t[]>(
          dev, max_grid_size_ * std::max(max_nlefs, max_nbarriers))),
      sorting_offset1_buff(cuda::memory::device::make_unique<uint32_t[]>(dev, max_grid_size_ + 1)),
      sorting_offset2_buff(cuda::memory::device::make_unique<uint32_t[]>(dev, max_grid_size_ + 1)),
      large_uint_buff_chunk_alignment(static_cast<uint32_t>(std::max(max_nlefs, max_nbarriers))),
      _device(dev) {
  assert(this->grid_size > 0);   // NOLINT
  assert(this->block_size > 0);  // NOLINT
  this->_device.set_resource_limit(cudaLimitMallocHeapSize, device_heap_size);
  this->sync_state_with_device();
  auto ok_buff = cuda::memory::device::allocate(this->_device, sizeof(int));
  cuda::memory::device::set(ok_buff, true);

  const auto global_state_dev = this->get_copy_of_device_instance();

  kernels::allocate_block_states<<<this->max_grid_size, this->block_size>>>(
      global_state_dev.block_states, static_cast<uint32_t>(max_nlefs),
      static_cast<uint32_t>(max_nbarriers), static_cast<int*>(ok_buff.get()));
  this->_device.synchronize();
  int ok;
  cuda::memory::copy(&ok, ok_buff.get(), sizeof(int));
  cuda::memory::device::free(ok_buff);
  if (!ok) {
    throw std::runtime_error("Failed to initialize block states on the device");
  }
  this->sync_state_with_device();
}

GlobalStateHost::~GlobalStateHost() {
  kernels::free_block_states<<<this->max_grid_size, this->block_size>>>(this->block_states);
  this->_device.synchronize();

  cuda::memory::device::free(this->block_states);
  cuda::memory::device::free(this->tasks);
  cuda::memory::device::free(this->_global_state_dev);
}

void GlobalStateHost::sync_state_with_device(const GlobalStateDev& state) {
  cuda::memory::copy_single(this->get_ptr_to_dev_instance(), &state);
}

void GlobalStateHost::sync_state_with_device() {
  this->sync_state_with_device(this->get_copy_of_device_instance());
}

void GlobalStateHost::write_tasks_to_device(const std::vector<Task>& new_tasks) {
  assert(!new_tasks.empty());  // NOLINT
  try {
    auto dev_state = this->get_copy_of_device_instance();
    cuda::memory::copy(dev_state.tasks, new_tasks.data(), new_tasks.size() * sizeof(Task));
    dev_state.ntasks = static_cast<uint32_t>(new_tasks.size());

    this->sync_state_with_device(dev_state);
  } catch (const cuda::runtime_error& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Falied to write tasks to device global memory: {}"), e.what()));
  }
}

void GlobalStateHost::write_barriers_to_device(
    const std::vector<uint32_t>& new_barrier_positions,
    const std::vector<dna::Direction>& new_barrier_directions) {
  assert(new_barrier_positions.size() == new_barrier_directions.size());  // NOLINT
  try {
    auto dev_state = this->get_copy_of_device_instance();
    cuda::memory::copy(dev_state.barrier_pos, new_barrier_positions.data(),
                       new_barrier_positions.size() * sizeof(uint32_t));
    cuda::memory::copy(dev_state.barrier_directions, new_barrier_directions.data(),
                       new_barrier_directions.size() * sizeof(dna::Direction));

    this->sync_state_with_device(dev_state);
  } catch (const cuda::runtime_error& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Falied to write extrusion barriers to device global memory: {}"), e.what()));
  }
}

void GlobalStateHost::init(size_t grid_size_, size_t block_size_,
                           const std::vector<Task>& tasks_host,
                           const std::vector<uint32_t>& barrier_pos_host,
                           const std::vector<dna::Direction>& barrier_dir_host) {
  assert(barrier_pos_host.size() == barrier_dir_host.size());  // NOLINT
  this->grid_size = grid_size_;
  this->block_size = block_size_;
  if (!tasks_host.empty()) {
    this->write_tasks_to_device(tasks_host);
  }

  if (!barrier_pos_host.empty()) {
    this->write_barriers_to_device(barrier_pos_host, barrier_dir_host);
  }

  this->sync_state_with_device();

  kernels::init_curand<<<this->grid_size, this->block_size>>>(this->get_ptr_to_dev_instance());
  this->_device.synchronize();

  kernels::reset_buffers<<<this->grid_size, this->block_size>>>(this->get_ptr_to_dev_instance());
  this->_device.synchronize();

  /*
  BlockState foo;
  Simulation::Task foo2;
  cuda::memory::copy_single(&foo2, this->get_copy_of_device_instance().tasks);
  std::vector<uint32_t> foo1(foo2.nlefs);

  cuda::memory::copy_single(&foo, this->get_copy_of_device_instance().block_states);
  cuda::memory::copy(foo1.data(), foo.lef_rev_unit_idx, 1 * sizeof(uint32_t));

  for (auto& n : foo1) {
    fmt::print(stderr, "{}\t", n);
  }
  fmt::print(stderr, "\n");
   */

  /*
  GlobalStateDev global_state_host;
  CUDA_CALL(cudaMemcpy(&global_state_host, global_state_dev, sizeof(GlobalStateDev),
                       cudaMemcpyDeviceToHost));
  global_state_host.nblock_states = static_cast<uint32_t>(host_tasks.size());
  global_state_host.ntasks = static_cast<uint32_t>(host_tasks.size());
  CUDA_CALL(cudaMemcpy(global_state_host.barrier_pos, barrier_pos.data(),
                       barrier_pos.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(global_state_host.barrier_directions, barrier_dir.data(),
                       barrier_dir.size() * sizeof(dna::Direction), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(global_state_host.tasks, host_tasks.data(), host_tasks.size() * sizeof(Task),
                       cudaMemcpyHostToDevice));

  CUDA_CALL(cudaMemcpy(global_state_dev, &global_state_host, sizeof(GlobalStateDev),
                       cudaMemcpyHostToDevice));

  init_curand_kernel<<<grid_size, block_size>>>(global_state_dev);

  reset_buffers_kernel<<<grid_size, block_size>>>(global_state_dev);
  CUDA_CALL(cudaDeviceSynchronize());
   */
}

Config* write_config_to_device(const Config& c) {
  // TODO Rewrite using cuda::memory::locate when it is released in a stable release
  // https://github.com/eyalroz/cuda-api-wrappers/commit/d07c19f8488c5ea028f716b1fc219c84ee4bbaf6

  if (cudaMemcpyToSymbol(config, &c, sizeof(modle::cu::Config)) != cudaSuccess) {
    throw std::runtime_error("Failed to write simulation parameters to device constant memory");
  }

  cudaDeviceSynchronize();
  void* symbol_ptr{};
  if (cudaGetSymbolAddress(&symbol_ptr, config) != cudaSuccess) {
    throw std::runtime_error(
        "Failed to retrieve the address of the simulation parameters from the device");
  }

  return static_cast<Config*>(symbol_ptr);
}

}  // namespace modle::cu
