// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/variant.h>  // for monostate, variant

#include <filesystem>  // for path
#include <string>      // for string, allocator, basic_string
#include <thread>      // for thread
#include <vector>      // for vector

#include "modle/bed/bed.hpp"           // for BED, BED::Dialect, BED::BED6
#include "modle/common/cli_utils.hpp"  // for CliEnumMappings
#include "modle/common/common.hpp"     // for usize, bp_t, u64
#include "modle/common/utils.hpp"      // for ConstMap

namespace modle::tools {

struct eval_config {
  // IO
  std::filesystem::path path_to_input_matrix;
  std::filesystem::path output_prefix;
  std::filesystem::path path_to_reference_matrix;
  std::filesystem::path path_to_chrom_subranges;
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

  // Reference contacts
  usize bin_size{0};
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
  std::filesystem::path path_to_input_matrix;
  std::filesystem::path path_to_output_matrix;
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
  usize bin_size{0};
  usize diagonal_width;

  // Other
  usize nthreads{std::thread::hardware_concurrency()};
  absl::Span<char*> args;
  std::string args_json{};
};

// clang-format off
using modle_tools_config = absl::variant<absl::monostate,
                                         eval_config,
                                         transform_config>;
// clang-format on

}  // namespace modle::tools
