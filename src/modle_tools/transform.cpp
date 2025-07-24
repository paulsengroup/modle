// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/time/clock.h>
#include <absl/time/time.h>
#include <absl/types/span.h>
#include <cpp-sort/comparators/natural_less.h>
#include <fmt/format.h>
#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>
#include <readerwriterqueue/readerwriterqueue.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <future>
#include <hictk/cooler.hpp>
#include <iosfwd>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "modle/bed/bed.hpp"
#include "modle/bigwig/bigwig.hpp"
#include "modle/common/common.hpp"
#include "modle/common/fmt_helpers.hpp"
#include "modle/common/string_utils.hpp"
#include "modle/common/utils.hpp"
#include "modle/compressed_io/compressed_io.hpp"
#include "modle/config/version.hpp"
#include "modle/contact_matrix_dense.hpp"
#include "modle/interval_tree.hpp"
#include "modle/io/contact_matrix_dense.hpp"
#include "modle/stats/correlation.hpp"
#include "modle_tools/modle_tools_config.hpp"
#include "modle_tools/tools.hpp"

namespace modle::tools {

static modle::IITree<double, double> import_discretization_ranges(const std::filesystem::path& p,
                                                                  const char sep = '\t') {
  IITree<double, double> ranges;
  if (p.empty()) {
    return ranges;
  }

  compressed_io::Reader r(p);

  std::string buff;
  std::vector<std::string_view> toks(3);
  usize i;  // NOLINT(cppcoreguidelines-init-variables)
  for (i = 1; r.getline(buff); ++i) {
    toks = str_split(buff, sep);
    try {
      if (toks.size() != 3) {
        throw std::runtime_error(
            fmt::format("expected 3 tokens, got {}. "
                        "Invalid record: \"{}\"",
                        i, p, toks.size(), buff));
      }

      const auto start = utils::parse_numeric_or_throw<double>(toks[0]);
      const auto end = utils::parse_numeric_or_throw<double>(toks[1]);
      const auto value = utils::parse_numeric_or_throw<double>(toks[2]);

      if (start >= end) {
        throw std::runtime_error(
            fmt::format("the begin of a range should be less or equal than its end, "
                        "found begin={}; end={}. Invalid record: \"{}\"",
                        start, end, buff));
      }
      ranges.insert(start, end, value);
    } catch (const std::exception& e) {
      throw std::runtime_error(
          fmt::format("found invalid record at line {} of file {}: {}", i, p, e.what()));
    }
  }
  ranges.make_BST();
  SPDLOG_INFO("imported {} ranges from file {}", i - 1, p);
  return ranges;
}

template <class N>
[[nodiscard]] static ContactMatrixDense<N> run_normalization(
    std::string_view chrom_name, ContactMatrixDense<N>& m,
    const std::pair<double, double> normalization_range,
    const std::pair<double, double> saturation_range) {
  const auto [lower_bound_norm, upper_bound_norm] = normalization_range;
  const auto [lower_bound_sat, upper_bound_sat] = saturation_range;

  SPDLOG_INFO("normalizing contacts for {} to the {:.4g}-{:.4g} range...", chrom_name,
              lower_bound_norm, upper_bound_norm);
  // Clamp contacts before normalizing
  if (!std::isinf(lower_bound_sat) || !std::isinf(upper_bound_sat)) {
    m.clamp_inplace(static_cast<N>(lower_bound_sat), static_cast<N>(upper_bound_sat));
  }
  m.normalize_inplace(lower_bound_norm, upper_bound_norm);
  return m;
}

template <class N>
[[nodiscard]] static ContactMatrixDense<double> run_gaussian_blur(
    BS::light_thread_pool& tpool, std::string_view chrom_name, ContactMatrixDense<N>& m,
    const double sigma, const std::pair<double, double> saturation_range) {
  const auto [lower_bound_sat, upper_bound_sat] = saturation_range;
  SPDLOG_INFO("applying Gaussian blur with sigma={:.4g} to contacts for {}...", sigma, chrom_name);
  // TODO: make truncate tunable
  auto m2 = m.blur(sigma, 3.5, &tpool);
  if (!std::isinf(lower_bound_sat) || !std::isinf(upper_bound_sat)) {
    m2.clamp_inplace(lower_bound_sat, upper_bound_sat);
  }
  return m2;
}

template <class N>
[[nodiscard]] static ContactMatrixDense<double> run_difference_of_gaussians(
    BS::light_thread_pool& tpool, std::string_view chrom_name, ContactMatrixDense<N>& m,
    const double sigma1, const double sigma2, const std::pair<double, double> saturation_range) {
  const auto [lower_bound_sat, upper_bound_sat] = saturation_range;
  SPDLOG_INFO(
      "computing the difference of Gaussians for {} (sigma1={:.4g}; "
      "sigma2={:.4g})...",
      chrom_name, sigma1, sigma2);
  // TODO: make truncate tunable
  return m.diff_of_gaussians(sigma1, sigma2, 3.5, lower_bound_sat, upper_bound_sat, &tpool);
}

[[nodiscard]] static ContactMatrixDense<double> process_chromosome(
    BS::light_thread_pool& tpool, const std::string_view chrom_name,
    const hictk::cooler::File& cooler, const modle::tools::transform_config& c,
    const modle::IITree<double, double>& discretization_ranges) {
  auto t0 = absl::Now();
  SPDLOG_INFO("processing contacts for {}...", chrom_name);
  auto matrix = [&]() {
    auto m = io::read_contact_matrix_from_cooler<double>(cooler, chrom_name, bp_t(0),
                                                         std::numeric_limits<bp_t>::max(),
                                                         static_cast<bp_t>(c.diagonal_width));
    using t = transform_config::Transformation;
    switch (c.method) {
      case t::normalize:
        return run_normalization(chrom_name, m, c.normalization_range, c.saturation_range);
      case t::gaussian_blur:
        return run_gaussian_blur(tpool, chrom_name, m, c.gaussian_blur_sigma, c.saturation_range);
      case t::difference_of_gaussians:
        return run_difference_of_gaussians(tpool, chrom_name, m, c.gaussian_blur_sigma,
                                           c.gaussian_blur_sigma * c.gaussian_blur_sigma_multiplier,
                                           c.saturation_range);
      default:
        utils::unreachable_code();
    }
  }();

  if (!discretization_ranges.empty()) {
    matrix.discretize_inplace(discretization_ranges);
  }

  SPDLOG_INFO("{} processing took {}", chrom_name, absl::FormatDuration(absl::Now() - t0));
  return matrix;
}

[[nodiscard]] static hictk::cooler::File init_output_cooler(
    const std::filesystem::path& output_cooler_uri, const hictk::cooler::File& input_cooler,
    bool floating_point, std::string_view args_json, bool force) {
  auto attrs = floating_point ? hictk::cooler::Attributes::init<double>(input_cooler.resolution())
                              : hictk::cooler::Attributes::init<i32>(input_cooler.resolution());
  attrs.metadata = args_json;
  attrs.generated_by = config::version::str_long("MoDLE-tools");
  if (const auto& assembly = input_cooler.attributes().assembly; assembly) {
    attrs.assembly = *assembly;
  }

  if (floating_point) {
    return io::init_cooler_file<double>(output_cooler_uri.string(), force,
                                        input_cooler.chromosomes(), std::move(attrs));
  }
  return io::init_cooler_file<i32>(output_cooler_uri.string(), force, input_cooler.chromosomes(),
                                   std::move(attrs));
}

void transform_subcmd(const modle::tools::transform_config& c) {
  const auto discretization_ranges = [&]() {
    if (!std::isnan(c.discretization_val)) {
      IITree<double, double> ranges;
      ranges.insert(std::numeric_limits<double>::lowest(), c.discretization_val, 0.0);
      ranges.insert(c.discretization_val, (std::numeric_limits<double>::max)(), 1.0);
      ranges.make_BST();
      return ranges;
    }
    return import_discretization_ranges(c.path_to_discretization_ranges_tsv);
  }();

  const hictk::cooler::File input_cooler(c.input_cooler_uri.string());
  if (const auto& output_dir = c.output_cooler_uri.parent_path(); !output_dir.empty()) {
    std::filesystem::create_directories(output_dir.string());
  }

  auto output_cooler =
      init_output_cooler(c.output_cooler_uri, input_cooler, c.floating_point, c.args_json, c.force);

  BS::light_thread_pool tpool(c.nthreads);
  const auto t0 = absl::Now();
  SPDLOG_INFO("transforming contacts from Cooler at URI \"{}\"...", input_cooler.uri());
  for (const auto& chrom : input_cooler.chromosomes()) {
    if (input_cooler.fetch(chrom.name()).empty()) {
      SPDLOG_WARN("read 0 contacts for {}. SKIPPING!", chrom.name());
      continue;
    }
    const auto transformed_matrix =
        process_chromosome(tpool, chrom.name(), input_cooler, c, discretization_ranges);
    if (c.floating_point) {
      io::append_contact_matrix_to_cooler(output_cooler, chrom.name(), transformed_matrix);
    } else {
      io::append_contact_matrix_to_cooler(output_cooler, chrom.name(),
                                          transformed_matrix.as<i32>());
    }
  }

  tpool.wait();
  SPDLOG_INFO("DONE! Processed {} chromosomes in {}!", input_cooler.chromosomes().size(),
              absl::FormatDuration(absl::Now() - t0));
  SPDLOG_INFO("Transformed contacts have been saved to file {}", c.output_cooler_uri);
}

}  // namespace modle::tools
