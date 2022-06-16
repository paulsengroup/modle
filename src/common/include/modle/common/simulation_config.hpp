// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span

#include <bitflags/bitflags.hpp>
#include <cmath>       // for round
#include <filesystem>  // for path
#include <limits>      // for numeric_limits
#include <string>      // for string
#include <thread>      // for thread
#include <vector>      // for vector

#include "modle/common/common.hpp"  // for bp_t, i64, u64

namespace modle {
class Cli;
class Simulation;

// NOLINTNEXTLINE(altera-struct-pack-align,clang-analyzer-optin.performance.Padding)
struct Config {
  friend Cli;
  friend Simulation;

  enum class StoppingCriterion : u8f { contact_density, simulation_epochs };

  // Even though we don't have any use for a none flag, it is required in order
  // for the automatically generated enum values make sense.
  BEGIN_BITFLAGS(ContactSamplingStrategy)
  FLAG(none)
  FLAG(noisify)
  FLAG(tad)
  FLAG(loop)
  END_BITFLAGS(ContactSamplingStrategy)

  std::filesystem::path path_to_chrom_sizes;
  std::filesystem::path path_to_chrom_subranges;
  std::filesystem::path path_to_output_prefix;
  std::filesystem::path path_to_output_file_cool{};
  std::filesystem::path path_to_config_file{};
  std::filesystem::path path_to_log_file;
  std::filesystem::path path_to_model_state_log_file;
  std::filesystem::path path_to_extr_barriers;
  bool force{false};
  bool quiet{false};
  std::filesystem::path path_to_reference_contacts{};
  std::vector<std::filesystem::path> path_to_feature_bed_files{};
  std::filesystem::path path_to_output_file_bedpe{};
  std::filesystem::path path_to_deletion_bed{};
  std::filesystem::path path_to_task_file;
  std::filesystem::path path_to_task_filter_file{};
  bool write_header{true};
  bool skip_output{false};
  bool log_model_internal_state{false};

  // Stopping criteria
  usize target_simulation_epochs{2000};
  double target_contact_density{1.0};
  StoppingCriterion stopping_criterion{StoppingCriterion::contact_density};

  // Contact matrix and sampling params
  bp_t bin_size{5'000};
  bp_t diagonal_width{3'000'000 /* 3 Mbp */};
  ContactSamplingStrategy contact_sampling_strategy{ContactSamplingStrategy::tad |
                                                    ContactSamplingStrategy::loop |
                                                    ContactSamplingStrategy::noisify};
  double tad_to_loop_contact_ratio{5.0};
  double genextreme_mu{0};
  double genextreme_sigma{5'000};
  double genextreme_xi{0.001};

  // LEFs params
  bp_t fwd_extrusion_speed{bin_size * 8 / 10};  // 80% of the bin size
  bp_t rev_extrusion_speed{fwd_extrusion_speed};
  double fwd_extrusion_speed_std{0.05};
  double rev_extrusion_speed_std{fwd_extrusion_speed_std};
  double number_of_lefs_per_mbp{20};
  double prob_of_lef_release{};
  double prob_of_lef_release_burnin{};
  bp_t avg_lef_processivity{300'000};
  bp_t contact_sampling_interval{50'000};

  // Extrusion barrier params
  double extrusion_barrier_occupancy{0.825};
  bool override_extrusion_barrier_occupancy{false};
  double barrier_occupied_stp{0.0};
  double barrier_not_occupied_stp{0.70};

  // Collision/stall params
  double hard_stall_lef_stability_multiplier{5.0};
  double soft_stall_lef_stability_multiplier{1.0};
  double probability_of_extrusion_unit_bypass{0.1};
  double lef_bar_major_collision_pblock{1.0};
  double lef_bar_minor_collision_pblock{0.0};

  // Miscellaneous
  bool simulate_chromosomes_wo_barriers{false};
  usize num_cells{512};
  usize nthreads{std::thread::hardware_concurrency()};
  u64 seed{0};
  bp_t probability_normalization_factor{rev_extrusion_speed + fwd_extrusion_speed};
  bool normalize_probabilities{true};

  // MoDLE perturbate
  bp_t deletion_size{10'000};
  bool compute_reference_matrix{false};
  usize block_size{9};

  // Burn-in
  bool skip_burnin{false};
  usize burnin_history_length{100};
  usize burnin_smoothing_window_size{5};
  usize min_burnin_epochs{0};
  usize max_burnin_epochs{(std::numeric_limits<usize>::max)()};
  usize burnin_target_epochs_for_lef_activation{320};
  double burnin_speed_coefficient{1.0};
  bp_t fwd_extrusion_speed_burnin{bp_t(double(fwd_extrusion_speed) * burnin_speed_coefficient)};
  bp_t rev_extrusion_speed_burnin{fwd_extrusion_speed_burnin};

  static constexpr std::string_view model_internal_state_log_header =
      "task_id\tepoch\tcell_id\t"
      "chrom\tburnin\tbarrier_occupancy\t"
      "num_active_lefs\tnum_stalls_rev\tnum_stalls_fwd\t"
      "num_stalls_both\tnum_lef_bar_collisions\tnum_primary_lef_lef_collisions\t"
      "num_secondary_lef_lef_collisions\tavg_loop_size\n";

  absl::Span<char*> args;
  std::string argv_json{};
};

}  // namespace modle
