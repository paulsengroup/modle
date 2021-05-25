#include <curand_kernel.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cub/cub.cuh>
#include <cuda/runtime_api.hpp>
#include <numeric>
#include <string>
#include <vector>

#include "modle/cu/common.hpp"
#include "modle/cu/config.hpp"
#include "modle/cu/simulation.hpp"
#include "modle/cu/simulation_internal.cuh"

namespace modle::cu {

__constant__ struct Config config;  // NOLINT

size_t GlobalStateHost::get_grid_size() const { return this->grid_size; };
size_t GlobalStateHost::get_block_size() const { return this->block_size; };

GlobalStateDev* GlobalStateHost::get_ptr_to_dev_instance() {
  if (!this->_global_state_dev.get()) {
    this->_global_state_dev = cuda::memory::device::allocate(this->_device, sizeof(GlobalStateDev));
    GlobalStateDev gs{};

    gs._config = this->_config;

    gs.grid_size = this->grid_size;
    gs.block_size = this->block_size;
    gs.shared_memory_bytes_per_block = this->shared_memory_bytes_per_block;

    gs.block_states = this->block_states;
    gs.tasks = this->tasks;

    gs.current_epoch = this->current_epoch;
    gs.nblock_states = this->nblock_states;
    gs.ntasks = this->ntasks;
    gs.ntasks_completed = this->ntasks_completed;

    gs.barrier_pos = this->barrier_pos.get();
    gs.barrier_directions = this->barrier_directions.get();
    gs.barrier_probs_occ_to_occ = this->barrier_probs_occ_to_occ.get();
    gs.barrier_probs_nocc_to_nocc = this->barrier_probs_nocc_to_nocc.get();

    cuda::memory::copy(this->_global_state_dev.get(), &gs, sizeof(GlobalStateDev));
  }
  return static_cast<GlobalStateDev*>(this->_global_state_dev.get());
}

GlobalStateDev GlobalStateHost::get_copy_of_device_instance() {
  GlobalStateDev state_host;
  cuda::memory::copy_single(&state_host, this->get_ptr_to_dev_instance());

  return state_host;
}

GlobalStateHost::GlobalStateHost(size_t grid_size_, size_t block_size_, size_t max_grid_size_,
                                 size_t device_heap_size, const Config& c, size_t max_nlefs,
                                 size_t max_nbarriers, size_t max_ncontacts_per_block,
                                 cuda::device_t dev)
    : _config(write_config_to_device(c)),
      grid_size(grid_size_),
      block_size(block_size_),
      max_grid_size(max_grid_size_),
      block_states(allocate_block_states(dev, max_nbarriers, max_nlefs, max_ncontacts_per_block,
                                         max_grid_size, block_size)),
      tasks(static_cast<Task*>(
          cuda::memory::device::allocate(dev, max_grid_size_ * sizeof(Task)).get())),
      barrier_pos(cuda::memory::device::make_unique<bp_t[]>(dev, max_nbarriers * sizeof(bp_t))),
      barrier_directions(cuda::memory::device::make_unique<dna::Direction[]>(
          dev, max_nbarriers * sizeof(dna::Direction))),
      barrier_probs_occ_to_occ(
          cuda::memory::device::make_unique<float[]>(dev, max_nbarriers * sizeof(float))),
      barrier_probs_nocc_to_nocc(
          cuda::memory::device::make_unique<float[]>(dev, max_nbarriers * sizeof(float))),
      nblock_states(static_cast<uint32_t>(max_grid_size)),
      ntasks(static_cast<uint32_t>(max_grid_size)),
      _device(dev),
      _max_ncontacts_per_block(max_ncontacts_per_block) {
  assert(this->grid_size > 0);   // NOLINT
  assert(this->block_size > 0);  // NOLINT
  this->_device.set_resource_limit(cudaLimitMallocHeapSize, device_heap_size);
  this->sync_state_with_device();
}

GlobalStateHost::~GlobalStateHost() {
  deallocate_block_states(this->block_states, this->max_grid_size);
  cuda::memory::device::free(this->block_states);
  cuda::memory::device::free(this->tasks);
  cuda::memory::device::free(this->_global_state_dev);
}

void GlobalStateHost::sync_state_with_device(const GlobalStateDev& state) {
  cuda::memory::copy_single(this->get_ptr_to_dev_instance(), &state);
}

void GlobalStateHost::sync_state_with_device() {
  GlobalStateDev gs;

  gs._config = this->_config;

  gs.grid_size = this->grid_size;
  gs.block_size = this->block_size;
  gs.shared_memory_bytes_per_block = this->shared_memory_bytes_per_block;

  gs.block_states = this->block_states;
  gs.tasks = this->tasks;

  gs.current_epoch = this->current_epoch;
  gs.nblock_states = this->nblock_states;
  gs.ntasks = this->ntasks;
  gs.ntasks_completed = this->ntasks_completed;

  gs.barrier_pos = this->barrier_pos.get();
  gs.barrier_directions = this->barrier_directions.get();
  gs.barrier_probs_occ_to_occ = this->barrier_probs_occ_to_occ.get();
  gs.barrier_probs_nocc_to_nocc = this->barrier_probs_nocc_to_nocc.get();

  this->sync_state_with_device(gs);
}

void GlobalStateHost::write_tasks_to_device(const std::vector<Task>& new_tasks) {
  assert(!new_tasks.empty());  // NOLINT
  try {
    auto dev_state = this->get_copy_of_device_instance();
    cuda::memory::copy(dev_state.tasks, new_tasks.data(), new_tasks.size() * sizeof(Task));

    this->ntasks = static_cast<uint32_t>(new_tasks.size());
    dev_state.ntasks = this->ntasks;

    this->sync_state_with_device(dev_state);
  } catch (const cuda::runtime_error& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Falied to write tasks to device global memory: {}"), e.what()));
  }
}

void GlobalStateHost::write_barriers_to_device(
    const std::vector<uint32_t>& new_barrier_positions,
    const std::vector<dna::Direction>& new_barrier_directions,
    const std::vector<float>& new_barrier_probs_occ_to_occ_host,
    const std::vector<float>& new_barrier_probs_nocc_to_nocc_host) {
  const auto nbarriers = new_barrier_positions.size();
  assert(new_barrier_directions.size() == nbarriers);               // NOLINT
  assert(new_barrier_probs_occ_to_occ_host.size() == nbarriers);    // NOLINT
  assert(new_barrier_probs_nocc_to_nocc_host.size() == nbarriers);  // NOLINT

  try {
    auto dev_state = this->get_copy_of_device_instance();

    cuda::memory::copy(dev_state.barrier_pos, new_barrier_positions.data(),
                       nbarriers * sizeof(uint32_t));
    cuda::memory::copy(dev_state.barrier_directions, new_barrier_directions.data(),
                       nbarriers * sizeof(dna::Direction));
    cuda::memory::copy(dev_state.barrier_probs_occ_to_occ, new_barrier_probs_occ_to_occ_host.data(),
                       nbarriers * sizeof(float));
    cuda::memory::copy(dev_state.barrier_probs_nocc_to_nocc,
                       new_barrier_probs_nocc_to_nocc_host.data(), nbarriers * sizeof(float));

    this->sync_state_with_device(dev_state);
  } catch (const cuda::runtime_error& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Falied to write extrusion barriers to device global memory: {}"), e.what()));
  }
}

void GlobalStateHost::init(size_t grid_size_, size_t block_size_,
                           const std::vector<Task>& tasks_host,
                           const std::vector<uint32_t>& barrier_pos_host,
                           const std::vector<dna::Direction>& barrier_dir_host,
                           const std::vector<float>& barrier_probs_occ_to_occ_host,
                           const std::vector<float>& barrier_probs_nocc_to_nocc_host) {
  assert(barrier_pos_host.size() == barrier_dir_host.size());  // NOLINT
  this->grid_size = grid_size_;
  this->block_size = block_size_;
  this->current_epoch = 0;
  this->ntasks_completed = 0;
  this->sync_state_with_device();

  if (!tasks_host.empty()) {
    this->write_tasks_to_device(tasks_host);
  }

  if (!barrier_pos_host.empty()) {
    this->write_barriers_to_device(barrier_pos_host, barrier_dir_host,
                                   barrier_probs_occ_to_occ_host, barrier_probs_nocc_to_nocc_host);
  }

  this->sync_state_with_device();

  kernels::init_curand<<<this->grid_size, this->block_size>>>(this->get_ptr_to_dev_instance());
  throw_on_cuda_error(cudaGetLastError(),
                      "An error occurred while initializing PRNG states on the device");

  kernels::reset_buffers<<<this->grid_size, this->block_size>>>(this->get_ptr_to_dev_instance());
  throw_on_cuda_error(cudaGetLastError(), "An error occurred while resetting device buffers");
  this->_device.synchronize();
}

BlockState* GlobalStateHost::allocate_block_states(const cuda::device_t& dev, size_t max_nbarriers,
                                                   size_t max_nlefs, size_t max_ncontacts,
                                                   size_t max_grid_size, size_t max_block_size) {
  auto allocate_or_throw = [&dev](auto& buff, uint32_t size) {
    using T = std::remove_pointer_t<std::decay_t<decltype(buff)>>;
    buff = static_cast<T*>(cuda::memory::device::allocate(dev, size * sizeof(T)).get());
  };

  auto* block_states_dev = static_cast<BlockState*>(
      cuda::memory::device::allocate(dev, max_grid_size * sizeof(BlockState)).get());
  std::vector<BlockState> block_states_host(max_grid_size);

  for (auto i = 0UL; i < max_grid_size; ++i) {
    allocate_or_throw(block_states_host[i].barrier_mask, max_nbarriers);

    allocate_or_throw(block_states_host[i].rev_unit_pos, max_nlefs);
    allocate_or_throw(block_states_host[i].fwd_unit_pos, max_nlefs);
    allocate_or_throw(block_states_host[i].lef_rev_unit_idx, max_nlefs);
    allocate_or_throw(block_states_host[i].lef_fwd_unit_idx, max_nlefs);

    allocate_or_throw(block_states_host[i].rev_moves_buff, max_nlefs);
    allocate_or_throw(block_states_host[i].fwd_moves_buff, max_nlefs);
    allocate_or_throw(block_states_host[i].rev_collision_mask, max_nlefs);
    allocate_or_throw(block_states_host[i].fwd_collision_mask, max_nlefs);

    allocate_or_throw(block_states_host[i].lef_unloader_affinities, max_nlefs);
    allocate_or_throw(block_states_host[i].lef_unloader_affinities_prefix_sum, max_nlefs + 1);
    allocate_or_throw(block_states_host[i].contact_local_buff, max_ncontacts);
    allocate_or_throw(block_states_host[i].loading_epochs, max_nlefs);
    allocate_or_throw(block_states_host[i].lefs_to_load_per_epoch, max_nlefs);

    allocate_or_throw(block_states_host[i].rng_state, max_block_size);

    block_states_host[i].contact_local_buff_capacity = max_ncontacts;

    size_t cub_tmp_storage_bytes_required = 0;
    auto status = cub::DeviceRunLengthEncode::Encode(
        nullptr, cub_tmp_storage_bytes_required, block_states_host[i].loading_epochs,
        block_states_host[i].loading_epochs, block_states_host[i].loading_epochs,
        block_states_host[i].loading_epochs, static_cast<int>(max_nlefs));
    assert(status == cudaSuccess);  // NOLINT
    block_states_host[i].cub_tmp_storage_bytes = cub_tmp_storage_bytes_required;

    cub::DoubleBuffer tmp_buff(block_states_host[i].rev_unit_pos,
                               block_states_host[i].lef_rev_unit_idx);

    status = cub::DeviceRadixSort::SortPairs(nullptr, block_states_host[i].cub_tmp_storage_bytes,
                                             tmp_buff, tmp_buff, static_cast<int>(max_nlefs));
    assert(status == cudaSuccess);  // NOLINT
    if (cub_tmp_storage_bytes_required > block_states_host[i].cub_tmp_storage_bytes) {
      block_states_host[i].cub_tmp_storage_bytes = cub_tmp_storage_bytes_required;
    }

    block_states_host[i].cub_tmp_storage =
        cuda::memory::device::allocate(dev, block_states_host[i].cub_tmp_storage_bytes).get();

    allocate_or_throw(block_states_host[i].tmp_lef_buff1, max_nlefs);
    allocate_or_throw(block_states_host[i].tmp_lef_buff2, max_nlefs);
    allocate_or_throw(block_states_host[i].tmp_lef_buff3, max_nlefs);
    allocate_or_throw(block_states_host[i].tmp_lef_buff4, max_nlefs);
  }

  cuda::memory::copy(block_states_dev, block_states_host.data(),
                     block_states_host.size() * sizeof(BlockState));

  return block_states_dev;
}

void GlobalStateHost::deallocate_block_states(BlockState* block_states_dev, size_t nstates) {
  BlockState block_state_host;
  for (auto i = 0UL; i < nstates; ++i) {
    cuda::memory::copy_single(&block_state_host, block_states_dev + i);
    cuda::memory::device::free(block_state_host.barrier_mask);

    cuda::memory::device::free(block_state_host.rev_unit_pos);
    cuda::memory::device::free(block_state_host.fwd_unit_pos);
    cuda::memory::device::free(block_state_host.lef_rev_unit_idx);
    cuda::memory::device::free(block_state_host.lef_fwd_unit_idx);

    cuda::memory::device::free(block_state_host.rev_moves_buff);
    cuda::memory::device::free(block_state_host.fwd_moves_buff);
    cuda::memory::device::free(block_state_host.rev_collision_mask);
    cuda::memory::device::free(block_state_host.fwd_collision_mask);

    cuda::memory::device::free(block_state_host.lef_unloader_affinities);
    cuda::memory::device::free(block_state_host.lef_unloader_affinities_prefix_sum);
    cuda::memory::device::free(block_state_host.contact_local_buff);
    cuda::memory::device::free(block_state_host.loading_epochs);

    cuda::memory::device::free(block_state_host.rng_state);

    cuda::memory::device::free(block_state_host.cub_tmp_storage);
    cuda::memory::device::free(block_state_host.tmp_lef_buff1);
    cuda::memory::device::free(block_state_host.tmp_lef_buff2);
    cuda::memory::device::free(block_state_host.tmp_lef_buff3);
    cuda::memory::device::free(block_state_host.tmp_lef_buff4);
  }
}

Config* write_config_to_device(const Config& c) {
  // TODO Rewrite using cuda::memory::locate when it is released in a
  // stable release
  // https://github.com/eyalroz/cuda-api-wrappers/commit/d07c19f8488c5ea028f716b1fc219c84ee4bbaf6

  if (cudaMemcpyToSymbol(config, &c, sizeof(modle::cu::Config)) != cudaSuccess) {
    throw std::runtime_error(
        "Failed to write simulation parameters to device constant "
        "memory");
  }

  cudaDeviceSynchronize();
  void* symbol_ptr{};
  if (cudaGetSymbolAddress(&symbol_ptr, config) != cudaSuccess) {
    throw std::runtime_error(
        "Failed to retrieve the address of the simulation parameters "
        "from the device");
  }

  return static_cast<Config*>(symbol_ptr);
}

}  // namespace modle::cu
