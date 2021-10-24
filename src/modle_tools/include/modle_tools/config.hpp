// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/container/flat_hash_set.h>  // for BitMask, flat_hash_set, raw_hash_set
#include <absl/types/variant.h>            // for monostate, variant

#include <algorithm>                        // for max
#include <boost/filesystem/operations.hpp>  // for temp_directory_path
#include <boost/filesystem/path.hpp>        // for path
#include <string>                           // for string, basic_string, allocator
#include <thread>                           // for thread
#include <vector>                           // for vector

#include "modle/bed.hpp"            // for BED, BED::Dialect, BED::BED6
#include "modle/common/common.hpp"  // for bp_t, u64

namespace modle::tools {

struct eval_config {  // NOLINT
  // IO
  boost::filesystem::path path_to_input_matrix;
  boost::filesystem::path output_prefix;
  boost::filesystem::path path_to_reference_matrix;
  boost::filesystem::path path_to_chrom_subranges;
  boost::filesystem::path path_to_features_bed;
  boost::filesystem::path tmp_dir{boost::filesystem::temp_directory_path()};
  bool force{false};

  // Correlation methods
  bool compute_spearman{true};
  bool compute_pearson{true};
  bool compute_edist{true};

  // Reference contacts
  usize bin_size{0};
  usize diagonal_width;

  // Other
  usize nthreads{std::thread::hardware_concurrency()};
};

struct filter_barrier_config {  // NOLINT
  // IO
  boost::filesystem::path path_to_extrusion_barrier_motifs_bed;
  std::vector<boost::filesystem::path> path_to_bed_files_for_filtering;

  // Other
  std::string filtering_criterion{"intersection"};
  bed::BED::Dialect bed_dialect{bed::BED::Dialect::BED6};
  bool strict_bed_validation{true};
};

struct find_barrier_clusters_config {  // NOLINT
  // IO
  boost::filesystem::path path_to_input_barriers;
  boost::filesystem::path path_to_output;
  boost::filesystem::path path_to_breaking_points;
  bool force{false};
  bool quiet{false};

  // Cluster properties
  bp_t extension_window{5000};  // NOLINT
  bp_t min_cluster_span{0};
  bp_t max_cluster_span{0};
  bp_t min_cluster_size{0};
  bp_t max_cluster_size{0};
};

struct noisify_config {  // NOLINT
  // IO
  boost::filesystem::path path_to_input_matrix;
  boost::filesystem::path path_to_output_matrix;
  bool force{false};

  // Generic
  usize diagonal_width;
  usize bin_size{0};

  // Noise properties
  double genextreme_mu{0};
  double genextreme_sigma{7'500};  // NOLINT
  double genextreme_xi{0.001};     // NOLINT
  u64 seed{0};
};

struct stats_config {  // NOLINT
  // IO
  boost::filesystem::path path_to_input_matrix;
  boost::filesystem::path path_to_chrom_subranges;
  std::string output_path_for_histograms{};
  bool dump_depleted_matrices{false};
  bool force{false};

  // Contact matrix
  usize bin_size{0};
  usize diagonal_width;
  std::vector<std::string> chromosomes_excluded_vect{};
  double depletion_multiplier{1.0};
};

// clang-format off
using config = absl::variant<absl::monostate,
                             eval_config,
                             filter_barrier_config,
                             find_barrier_clusters_config,
                             noisify_config,
                             stats_config>;
// clang-format on

}  // namespace modle::tools
