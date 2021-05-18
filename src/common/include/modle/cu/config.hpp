#pragma once

#include "modle/cu/common.hpp"  // for bp_t

namespace modle::cu {
struct Config {
  // clang-format off
  // General settings
  bp_t bin_size{};
  bp_t fwd_extrusion_speed{};
  bp_t rev_extrusion_speed{};
  float fwd_extrusion_speed_std{};
  float rev_extrusion_speed_std{};
  uint32_t ncells{};
  bp_t diagonal_width{};
  bp_t simulation_iterations{};
  float target_contact_density{};
  float number_of_lefs_per_mbp{};
  bp_t average_lef_lifetime{};
  float hard_stall_multiplier{};
  float soft_stall_multiplier{};
  bool skip_burnin{};
  uint64_t seed{};

  // Misc probabilities
  float probability_of_extrusion_unit_bypass{};
  float probability_of_extrusion_barrier_block{};
  float ctcf_occupied_self_prob{};
  float ctcf_not_occupied_self_prob{};
  float lef_hard_collision_pblock{};
  float lef_soft_collision_pblock{};
  float lef_fraction_contact_sampling{};
  // clang-format on
};
}  // namespace modle::cu
