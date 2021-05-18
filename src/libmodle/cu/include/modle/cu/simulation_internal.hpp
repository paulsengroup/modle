#pragma once
#include <curand_kernel.h>

#include <cstdint>

#include "modle/cu/common.hpp"
#include "modle/cu/simulation.hpp"

namespace modle::cu {

namespace kernels {

// Memory management
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
                                               uint32_t* tot_num_units_to_sort = nullptr);
__global__ void update_unit_mappings_and_scatter_sorted_lefs(
    GlobalStateDev* global_state, dna::Direction direction,
    bool update_extr_unit_to_lef_mappings = true);

__global__ void prepare_units_for_random_shuffling(GlobalStateDev* global_state,
                                                   uint32_t* tot_num_active_units);

__global__ void select_lefs_then_register_contacts(GlobalStateDev* global_state);

__global__ void generate_moves(GlobalStateDev* global_state);

__global__ void update_ctcf_states(GlobalStateDev* global_state);

__global__ void reset_collision_masks(GlobalStateDev* global_state);

__global__ void process_collisions(GlobalStateDev* global_state);
}  // namespace kernels

namespace dev {
__device__ uint32_t detect_collisions_at_5prime_single_threaded(
    const bp_t* rev_unit_pos, const bp_t* fwd_unit_pos, const bp_t* rev_moves,
    collision_t* rev_collision_mask, uint32_t chrom_start_pos, uint32_t num_extr_units);

__device__ uint32_t detect_collisions_at_3prime_single_threaded(
    const bp_t* rev_unit_pos, const bp_t* fwd_unit_pos, const bp_t* fwd_moves,
    collision_t* fwd_collision_mask, uint32_t chrom_end_pos, uint32_t num_extr_units);

__device__ void extrude(GlobalStateDev* global_state);
}  // namespace dev

}  // namespace modle::cu
