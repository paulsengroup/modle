// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/variant.h>  // for monostate, variant

#include <boost/filesystem/operations.hpp>  // for temp_directory_path
#include <boost/filesystem/path.hpp>        // for path
#include <string>                           // for string, allocator, basic_string
#include <thread>                           // for thread
#include <vector>                           // for vector

#include "modle/bed/bed.hpp"           // for BED, BED::Dialect, BED::BED6
#include "modle/common/cli_utils.hpp"  // for CliEnumMappings
#include "modle/common/common.hpp"     // for usize, bp_t, u64
#include "modle/common/utils.hpp"      // for ConstMap

namespace modle::tools {

struct eval_config {
  // IO
  boost::filesystem::path path_to_input_matrix;
  boost::filesystem::path output_prefix;
  boost::filesystem::path path_to_reference_matrix;
  boost::filesystem::path path_to_chrom_subranges;
  boost::filesystem::path path_to_weights;
  bool force{false};
  bool normalize{false};

  // Metrics
  enum class Metric : u8f { custom, eucl_dist, pearson, rmse, spearman };
  Metric metric{Metric::custom};

  // NOLINTNEXTLINE(cert-err58-cpp)
  inline static const utils::CliEnumMappings<enum Metric> metric_map{
      // clang-format off
      {"custom", Metric::custom},
      {"eucl_dist", Metric::eucl_dist},
      {"pearson", Metric::pearson},
      {"rmse", Metric::rmse},
      {"spearman", Metric::spearman}
      // clang-format on
  };

  // Reference contacts
  usize bin_size{0};
  usize diagonal_width;

  // Other
  usize nthreads{std::thread::hardware_concurrency()};
  bool exclude_zero_pxls{false};
  std::string weight_column_name{"balanced.sum"};
  bool reciprocal_weights{false};
};

struct find_barrier_clusters_config {
  // IO
  boost::filesystem::path path_to_input_barriers;
  boost::filesystem::path path_to_output;
  boost::filesystem::path path_to_breaking_points;
  bool force{false};
  bool quiet{false};

  // Cluster properties
  bp_t extension_window{5000};
  bp_t min_cluster_span{0};
  bp_t max_cluster_span{0};
  bp_t min_cluster_size{0};
  bp_t max_cluster_size{0};
};

struct noisify_config {
  // IO
  boost::filesystem::path path_to_input_matrix;
  boost::filesystem::path path_to_output_matrix;
  bool force{false};

  // Generic
  usize diagonal_width;
  usize bin_size{0};

  // Noise properties
  double genextreme_mu{0};
  double genextreme_sigma{7'500};
  double genextreme_xi{0.001};
  u64 seed{0};
};

struct stats_config {
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

struct transform_config {
  // IO
  boost::filesystem::path path_to_input_matrix;
  boost::filesystem::path path_to_output_matrix;
  bool force{false};

  // Transformation methods
  enum class Transformation : u8f { normalize, gaussian_blur, difference_of_gaussians, discretize };
  Transformation method;

  // NOLINTNEXTLINE(cert-err58-cpp)
  inline static const utils::CliEnumMappings<enum Transformation> transformation_map{
      // clang-format off
      {"normalize", Transformation::normalize},
      {"gaussian_blur", Transformation::gaussian_blur},
      {"difference_of_gaussians", Transformation::difference_of_gaussians},
      {"discretize", Transformation::discretize}
      // clang-format on
  };

  // Transformation settings
  std::pair<double, double> normalization_range{0.0, 1.0};
  std::pair<double, double> saturation_range{-std::numeric_limits<double>::infinity(),
                                             std::numeric_limits<double>::infinity()};
  double gaussian_blur_sigma{1.0};
  double gaussian_blur_sigma_multiplier{1.6};  // NOLINT, approx. Laplacian of Gaussian

  double discretization_val{std::numeric_limits<double>::quiet_NaN()};
  boost::filesystem::path path_to_discretization_ranges_tsv{};
  bool floating_point{true};

  // Matrix options
  usize bin_size{0};
  usize diagonal_width;

  // Other
  usize nthreads{std::thread::hardware_concurrency()};
};

// clang-format off
using modle_tools_config = absl::variant<absl::monostate,
                                         eval_config,
                                         find_barrier_clusters_config,
                                         noisify_config,
                                         stats_config,
                                         transform_config>;
// clang-format on

}  // namespace modle::tools
