// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/variant.h>

#include <filesystem>
#include <hictk/cooler/uri.hpp>
#include <string>
#include <thread>
#include <vector>

#include "modle/bed/bed.hpp"
#include "modle/common/cli_utils.hpp"
#include "modle/common/common.hpp"
#include "modle/common/utils.hpp"

namespace modle::tools {

struct eval_config {
  // IO
  std::filesystem::path input_cooler_uri;
  std::filesystem::path output_prefix;
  std::filesystem::path reference_cooler_uri;
  std::filesystem::path path_to_chrom_sizes;
  std::filesystem::path path_to_regions_of_interest_bed;
  std::filesystem::path path_to_weights;
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

  // Matrix options
  usize diagonal_width;

  // Other
  usize nthreads{std::thread::hardware_concurrency()};
  bool exclude_zero_pxls{false};
  std::string weight_column_name{"balanced.sum"};
  bool reciprocal_weights{false};
  absl::Span<char*> args;
  std::string args_json{};
};

struct transform_config {
  // IO
  std::filesystem::path input_cooler_uri;
  std::filesystem::path output_cooler_uri;
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
  std::filesystem::path path_to_discretization_ranges_tsv{};
  bool floating_point{true};

  // Matrix options
  usize diagonal_width;

  // Other
  usize nthreads{std::thread::hardware_concurrency()};
  absl::Span<char*> args;
  std::string args_json{};
};

struct annotate_barriers_config {
  // IO
  std::filesystem::path path_to_bigwig{};
  std::filesystem::path path_to_bed{};

  // Annotation settings
  double occupancy_lb{0.6};
  double occupancy_ub{1.0};
  bool clamp_occupancy{false};
  // See https://www.cell.com/cell-reports/pdf/S2211-1247(16)30530-7.pdf
  double scaling_factor{20 - 3};

  // Other
  absl::Span<char*> args;
  std::string args_json{};
};

// clang-format off
using modle_tools_config = absl::variant<absl::monostate,
                                         annotate_barriers_config,
                                         eval_config,
                                         transform_config>;
// clang-format on

}  // namespace modle::tools
