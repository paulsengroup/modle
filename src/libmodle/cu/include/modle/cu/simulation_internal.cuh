#pragma once
#include <cuda_runtime_api.h>
#include <curand_kernel.h>
#include <thrust/pair.h>

#include <cstdint>

#include "modle/cu/common.hpp"
#include "modle/cu/simulation.hpp"

namespace modle::cu {

namespace kernels {

// Memory management
__global__ void init_curand(GlobalStateDev* global_state);

__global__ void reset_buffers(GlobalStateDev* global_state);

// One off
__global__ void compute_initial_loading_epochs(GlobalStateDev* global_state);

__global__ void init_barrier_states(GlobalStateDev* global_state);

// Simulation body
__global__ void bind_and_sort_lefs(GlobalStateDev* global_state);

__global__ void select_lefs_then_register_contacts(GlobalStateDev* global_state);

__global__ void generate_moves(GlobalStateDev* global_state);

__global__ void update_ctcf_states(GlobalStateDev* global_state);

__global__ void reset_collision_masks(GlobalStateDev* global_state);

__global__ void process_collisions(GlobalStateDev* global_state);

__global__ void extrude_and_release_lefs(GlobalStateDev* global_state);

__global__ void advance_epoch(GlobalStateDev* global_state);
}  // namespace kernels

namespace dev {
__device__ uint32_t detect_collisions_at_5prime_single_threaded(
    const bp_t* rev_unit_pos, const bp_t* fwd_unit_pos, const bp_t* rev_moves,
    collision_t* rev_collision_mask, uint32_t chrom_start_pos, uint32_t num_extr_units);

__device__ uint32_t detect_collisions_at_3prime_single_threaded(
    const bp_t* rev_unit_pos, const bp_t* fwd_unit_pos, const bp_t* fwd_moves,
    collision_t* fwd_collision_mask, uint32_t chrom_end_pos, uint32_t num_extr_units);

__device__ void extrude(bp_t* rev_units_pos, bp_t* fwd_units_pos, const bp_t* rev_moves,
                        const bp_t* fwd_moves, uint32_t num_active_lefs, uint32_t chrom_start_pos,
                        uint32_t chrom_end_pos);

__device__ void detect_lef_bar_collisions(const bp_t* rev_units_pos, const bp_t* fwd_units_pos,
                                          const bp_t* rev_moves, const bp_t* fwd_moves,
                                          uint32_t num_active_lefs, const bp_t* barrier_pos,
                                          const dna::Direction* barrier_directions,
                                          const CTCF::State* barrier_states, uint32_t num_barriers,
                                          collision_t* rev_collisions, collision_t* fwd_collisions,
                                          curandStatePhilox4_32_10_t* rng_states,
                                          uint32_t num_rev_units_at_5prime,
                                          uint32_t num_fwd_units_at_3prime);

[[nodiscard]] __device__ bool bernoulli_trial(curandStatePhilox4_32_10_t* state,
                                              float prob_of_success);

__device__ void correct_moves_for_lef_bar_collisions(const bp_t* rev_units_pos,
                                                     const bp_t* fwd_units_pos, bp_t* rev_moves,
                                                     bp_t* fwd_moves, uint32_t num_active_lefs,
                                                     const bp_t* barr_pos, uint32_t num_barriers,
                                                     const collision_t* rev_collisions,
                                                     const collision_t* fwd_collisions);

__device__ void detect_primary_lef_lef_collisions(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, const uint32_t* lef_rev_unit_idx,
    const bp_t* rev_moves, const bp_t* fwd_moves, uint32_t num_active_lefs, const bp_t* barrier_pos,
    uint32_t num_barriers, collision_t* rev_collisions, collision_t* fwd_collisions,
    curandStatePhilox4_32_10_t* rng_states, uint32_t num_rev_units_at_5prime,
    uint32_t num_fwd_units_at_3prime);

__device__ void process_secondary_lef_lef_collisions(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, bp_t* rev_moves, bp_t* fwd_moves,
    uint32_t num_active_lefs, uint32_t num_barriers, collision_t* rev_collisions,
    collision_t* fwd_collisions, curandStatePhilox4_32_10_t* rng_states,
    uint32_t num_rev_units_at_5prime, uint32_t num_fwd_units_at_3prime);

__device__ void correct_moves_for_primary_lef_lef_collisions(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, bp_t* rev_moves, bp_t* fwd_moves,
    uint32_t num_active_lefs, uint32_t num_barriers, const collision_t* rev_collisions,
    const collision_t* fwd_collisions);

__device__ thrust::pair<bp_t, bp_t> compute_lef_lef_collision_pos(bp_t rev_unit_pos,
                                                                  bp_t fwd_unit_pos, bp_t rev_move,
                                                                  bp_t fwd_move);
__device__ void generate_lef_unloader_affinities(
    const bp_t* rev_units_pos, const bp_t* fwd_units_pos, const dna::Direction* barrier_directions,
    const uint32_t* lef_rev_idx, const collision_t* rev_collisions,
    const collision_t* fwd_collisions, uint32_t num_active_lefs, uint32_t num_barriers,
    float* lef_unloader_affinities);

__device__ void select_and_release_lefs(bp_t* rev_units_pos, bp_t* fwd_units_pos,
                                        const uint32_t* lef_rev_idx, const uint32_t* lef_fwd_idx,
                                        uint32_t num_active_lefs,
                                        const float* lef_unloader_affinities,
                                        float* lef_unloader_affinities_prefix_sum,
                                        void* tmp_storage, size_t tmp_storage_bytes,
                                        curandStatePhilox4_32_10_t* rng_states);

__device__ void generate_initial_loading_epochs(uint32_t* epoch_buff,
                                                curandStatePhilox4_32_10_t* rng_states,
                                                uint32_t num_lefs);

__device__ void select_and_bind_lefs(
    bp_t* rev_units_pos, bp_t* fwd_units_pos, const uint32_t* lef_rev_unit_idx,
    uint32_t* lef_fwd_unit_idx, const uint32_t* loading_epochs,
    const uint32_t* lefs_to_load_per_epoch, void* local_buff, uint32_t local_buff_bytes_per_block,
    uint32_t& epoch_idx, uint32_t num_unique_loading_epochs, uint32_t& num_active_lefs,
    uint32_t current_epoch, uint32_t tot_num_lefs, uint32_t chrom_start, uint32_t chrom_end,
    curandStatePhilox4_32_10_t* rng_states, bool& burnin_completed);

__device__ void update_extr_unit_mappings(uint32_t* lef_rev_unit_idx, uint32_t* lef_fwd_unit_idx,
                                          void* local_buff, uint32_t local_buff_bytes_per_block,
                                          uint32_t num_active_lefs, dna::Direction direction);

__device__ void reset_collision_masks(collision_t* rev_collisions, collision_t* fwd_collisions,
                                      uint32_t num_active_lefs);

__device__ thrust::pair<uint32_t, uint32_t> compute_chunk_size(uint32_t num_elements);

__device__ void shuffle_lefs(uint32_t* shuffled_lef_idx_buff, uint32_t* keys_buff,
                             uint32_t* tmp_lef_buff1, uint32_t* tmp_lef_buff2, void* tmp_storage,
                             size_t tmp_storage_bytes, uint32_t num_active_lefs,
                             curandStatePhilox4_32_10_t* rng_states);

__device__ void register_contacts(const bp_t* rev_unit_pos, const bp_t* fwd_unit_pos,
                                  const uint32_t* shuffled_idx, const uint32_t* lef_rev_unit_idx,
                                  uint2* contacts_buff, uint32_t chrom_start, uint32_t chrom_end,
                                  uint32_t bin_size, uint32_t sample_size, uint32_t num_active_lefs,
                                  uint32_t* contacts_buff_size_shared,
                                  uint32_t contacts_buff_capacity);

__device__ thrust::pair<int, int> compute_bit_width(uint32_t num);
}  // namespace dev

}  // namespace modle::cu
