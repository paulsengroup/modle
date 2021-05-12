#pragma once
#include <cstdint>

#include "modle/cu/common.hpp"
#include "modle/cu/simulation.hpp"

namespace modle::cu::kernels {

// Memory management
__global__ void allocate_block_states(BlockState* block_state, uint32_t max_nlefs,
                                      uint32_t max_nbarriers, int* ok);

__global__ void free_block_states(BlockState* block_state);

__global__ void init_curand(GlobalStateDev* global_state);

__global__ void reset_buffers(GlobalStateDev* global_state);

// One off
__global__ void generate_initial_loading_epochs(GlobalStateDev* global_state, uint32_t num_epochs);

__global__ void init_barrier_states(GlobalStateDev* global_state);

__global__ void shift_and_scatter_lef_loading_epochs(const uint32_t* global_buff,
                                                     BlockState* block_states,
                                                     const uint32_t* begin_offsets,
                                                     const uint32_t* end_offsets);

// Simulation body
__global__ void select_and_bind_lefs(uint32_t current_epoch, GlobalStateDev* global_state);

__global__ void prepare_extr_units_for_sorting(GlobalStateDev* global_state,
                                               dna::Direction direction,
                                               uint32_t* tot_num_items_to_sort);
__global__ void scatter_sorted_lefs(GlobalStateDev* global_state, dna::Direction direction,
                                    bool update_extr_unit_to_lef_mappings = true);
}  // namespace modle::cu::kernels
