#pragma once

#include <boost/filesystem/path.hpp>  //  for operator<<, path
#include <cstdint>                    // for uint32_t, uint64_t
#include <cstdio>                     // for stderr
#include <limits>                     // for numeric_limits
#include <string>                     // for string, allocator, basic_string
#include <thread>                     // for thread::hardware_concurrency

#include "modle/common/common.hpp"  // for bp_t

namespace modle {
using namespace std::literals::string_view_literals;
class Cli;
class Simulation;

struct Config {  // NOLINT(altera-struct-pack-align)
  friend Cli;
  friend Simulation;
  // clang-format off
  // IO
  boost::filesystem::path path_to_chrom_sizes;
  boost::filesystem::path path_to_chrom_subranges;
  boost::filesystem::path path_to_output_prefix;
  boost::filesystem::path path_to_output_file_cool{};
  boost::filesystem::path path_to_config_file{};
  boost::filesystem::path path_to_log_file;
  boost::filesystem::path path_to_extr_barriers;
  bool force{false};
  bool quiet{false};
  bool write_contacts_for_ko_chroms{false};
  bool write_tasks_to_disk{true};
  std::vector<boost::filesystem::path> path_to_feature_bed_files{};
  boost::filesystem::path path_to_output_file_bedpe{};
  boost::filesystem::path path_to_deletion_bed{};
  boost::filesystem::path path_to_task_file;
  boost::filesystem::path path_to_task_filter_file{};
  bool write_header{true};

  // General settings
  bp_t bin_size{1'000};                        // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  bp_t fwd_extrusion_speed{(std::numeric_limits<bp_t>::max)()};
  bp_t rev_extrusion_speed{(std::numeric_limits<bp_t>::max)()};
  double fwd_extrusion_speed_std{0.05};        // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double rev_extrusion_speed_std{fwd_extrusion_speed_std};
  size_t num_cells{5'000};                     // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  size_t nthreads{std::thread::hardware_concurrency()};
  bp_t diagonal_width{3'000'000 /* 3 Mbp */};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  size_t simulation_iterations{200};           // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double target_contact_density{0.0};
  double number_of_lefs_per_mbp;
  bp_t average_lef_lifetime{600'000};          // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double hard_stall_multiplier{5.0};           // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double soft_stall_multiplier{0.6};           // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  bp_t deletion_size{10'000};                  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  bool compute_reference_matrix{false};
  size_t block_size{9};                        // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  uint64_t seed{0};

  // Misc probabilities
  double probability_of_extrusion_unit_bypass{0.1};      // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double extrusion_barrier_occupancy{0.825};             // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double ctcf_occupied_self_prob{0.0};                   // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double ctcf_not_occupied_self_prob{0.70};              // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double lef_hard_collision_pblock{1.0};                 // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double lef_soft_collision_pblock{0.0};                 // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double lef_fraction_contact_sampling{0.025};           // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  bool randomize_contacts{false};
  double genextreme_mu{0};                               // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double genextreme_sigma{12'500};                       // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double genextreme_xi{0.001};                           // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)

  // Burn-in
  bool skip_burnin{false};
  size_t burnin_lef_binding_epochs{0};
  size_t burnin_epochs{0};
  double burnin_speed_coefficient{5.0};
  bp_t fwd_extrusion_speed_burnin{(std::numeric_limits<bp_t>::max)()};
  bp_t rev_extrusion_speed_burnin{(std::numeric_limits<bp_t>::max)()};

  // Misc
  bool exclude_chrom_wo_extr_barriers{true};
  bool skip_output{false};

  int argc;
  char** argv;
  // clang-format on
};

}  // namespace modle
