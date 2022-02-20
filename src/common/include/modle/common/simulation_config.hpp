// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>  // for Span

#include <bitflags/bitflags.hpp>
#include <boost/filesystem/path.hpp>  // for path
#include <cmath>                      // for round
#include <limits>                     // for numeric_limits
#include <string>                     // for string
#include <thread>                     // for thread
#include <vector>                     // for vector

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

  // clang-format off
  // IO
  boost::filesystem::path path_to_chrom_sizes;
  boost::filesystem::path path_to_chrom_subranges;
  boost::filesystem::path path_to_output_prefix;
  boost::filesystem::path path_to_output_file_cool{};
  boost::filesystem::path path_to_config_file{};
  boost::filesystem::path path_to_log_file;
  boost::filesystem::path path_to_model_state_log_file;
  boost::filesystem::path path_to_extr_barriers;
  bool force{false};
  bool quiet{false};
  bool write_contacts_for_ko_chroms{false};
  bool write_tasks_to_disk{true};
  boost::filesystem::path path_to_reference_contacts{};
  std::vector<boost::filesystem::path> path_to_feature_bed_files{};
  boost::filesystem::path path_to_output_file_bedpe{};
  boost::filesystem::path path_to_deletion_bed{};
  boost::filesystem::path path_to_task_file;
  boost::filesystem::path path_to_task_filter_file{};
  bool write_header{true};

  // General settings
  bp_t bin_size{5'000};
  bp_t fwd_extrusion_speed{bin_size / 2};
  bp_t rev_extrusion_speed{fwd_extrusion_speed};
  double fwd_extrusion_speed_std{0.05};
  double rev_extrusion_speed_std{fwd_extrusion_speed_std};
  usize num_cells{512};
  usize nthreads{std::thread::hardware_concurrency()};
  bp_t diagonal_width{3'000'000 /* 3 Mbp */};
  usize target_simulation_epochs{2000};
  double target_contact_density{1.0};
  StoppingCriterion stopping_criterion{StoppingCriterion::contact_density};
  ContactSamplingStrategy contact_sampling_strategy{ContactSamplingStrategy::tad | ContactSamplingStrategy::loop | ContactSamplingStrategy::noisify};
  u64 number_of_loop_contacts_per_sampling_event{1};
  u64 number_of_tad_contacts_per_sampling_event{5};
  double number_of_lefs_per_mbp{20};
  double prob_of_lef_release{0.015};
  double hard_stall_lef_stability_multiplier{5.0};
  double soft_stall_lef_stability_multiplier{1.0};
  bp_t deletion_size{10'000};
  bool compute_reference_matrix{false};
  usize block_size{9};
  u64 seed{0};

  // Misc probabilities
  double probability_of_extrusion_unit_bypass{0.1};
  double extrusion_barrier_occupancy{0.825};
  double ctcf_occupied_self_prob{0.0};
  double ctcf_not_occupied_self_prob{0.70};
  double lef_bar_major_collision_pblock{1.0};
  double lef_bar_minor_collision_pblock{0.0};
  double lef_fraction_contact_sampling{0.025};
  double genextreme_mu{0};
  double genextreme_sigma{5'000};
  double genextreme_xi{0.001};

  // Burn-in
  bool skip_burnin{false};
  usize burnin_history_length{100};
  usize burnin_smoothing_window_size{5};
  usize max_burnin_epochs{(std::numeric_limits<i64>::max)()};
  usize burnin_target_epochs_for_lef_activation{320};
  double burnin_speed_coefficient{1.0};
  bp_t fwd_extrusion_speed_burnin{bp_t(std::round(double(fwd_extrusion_speed) * burnin_speed_coefficient))};
  bp_t rev_extrusion_speed_burnin{fwd_extrusion_speed_burnin};

  // Misc
  bool exclude_chrom_wo_extr_barriers{true};
  bool skip_output{false};
  bool log_model_internal_state{false};

  absl::Span<char*> args;

  std::string argv_json{};
  // clang-format on
};

}  // namespace modle