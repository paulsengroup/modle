// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/algorithm/container.h>
#include <absl/strings/ascii.h>  // AsciiStrToLower
#include <absl/strings/str_replace.h>
#include <absl/strings/strip.h>
#include <absl/time/clock.h>
#include <absl/time/time.h>
#include <absl/types/span.h>
#include <cpp-sort/comparators/natural_less.h>
#include <fmt/compile.h>
#include <fmt/format.h>
#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>
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
#include <vector>

#include "modle/bed/bed.hpp"
#include "modle/bigwig/bigwig.hpp"
#include "modle/chrom_sizes/chrom_sizes.hpp"
#include "modle/common/common.hpp"
#include "modle/common/fmt_helpers.hpp"
#include "modle/common/utils.hpp"
#include "modle/contact_matrix_dense.hpp"
#include "modle/interval_tree.hpp"
#include "modle/io/contact_matrix_dense.hpp"
#include "modle/stats/correlation.hpp"
#include "modle_tools/modle_tools_config.hpp"
#include "modle_tools/tools.hpp"

namespace modle::tools {

[[nodiscard]] static hictk::Reference import_chroms_from_coolers(const hictk::cooler::File &f1,
                                                                 const hictk::cooler::File &f2) {
  std::vector<hictk::Chromosome> chroms{};
  for (const auto &chrom : f1.chromosomes()) {
    if (f2.chromosomes().contains(chrom.name()) &&
        f2.chromosomes().at(chrom.name()).size() == chrom.size()) {
      chroms.push_back(chrom);
    }
  }
  return {chroms.begin(), chroms.end()};
}

[[nodiscard]] static std::vector<bed::BED> chromset_to_bed(const hictk::Reference &chroms) {
  std::vector<bed::BED> intervals{chroms.size()};
  std::transform(
      chroms.begin(), chroms.end(), intervals.begin(),
      [&](const hictk::Chromosome &chrom) { return bed::BED{chrom.name(), 0, chrom.size()}; });
  return intervals;
}

[[nodiscard]] static hictk::Reference import_chroms_from_chrom_sizes_file(
    const std::filesystem::path &path_to_chrom_sizes) {
  std::vector<std::string> chrom_names{};
  std::vector<u32> chrom_sizes{};
  for (auto &&bed : chrom_sizes::Parser(path_to_chrom_sizes).parse_all()) {
    chrom_names.emplace_back(std::move(bed.chrom));
    chrom_sizes.emplace_back(bed.chrom_end - bed.chrom_start);
  }

  return {chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()};
}

[[nodiscard]] static std::vector<bed::BED> import_regions_of_interest(
    const hictk::Reference &chroms, const std::filesystem::path &path_to_regions_of_interest) {
  auto regions_of_interest = bed::Parser(path_to_regions_of_interest).parse_all();
  std::sort(regions_of_interest.begin(), regions_of_interest.end(),
            [&](const bed::BED &r1, const bed::BED &r2) {
              const auto m1 = chroms.find(r1.chrom);
              const auto m2 = chroms.find(r2.chrom);

              // both regions overlap with the chrom annotation
              if (m1 != chroms.end() && m2 != chroms.end()) {
                return *m1 < *m2;
              }

              // both regions DO NOT overlap with the chrom annotation
              if (m1 == chroms.end() && m2 == chroms.end()) {
                return m1->size() < m2->size();
              }
              // only r2 overlaps with the chrom annotation
              // return false so that m1 is pushed towards the end of the vector
              if (m1 == chroms.end()) {
                return false;
              }

              // only r1 overlaps with the chrom annotation
              // return true so that m2 is pushed towards the end of the vector
              return true;
            });

  return regions_of_interest;
}

static void validate_regions_of_interest(const hictk::Reference &chroms,
                                         const std::vector<bed::BED> &regions_of_interest) {
  usize invalid_chrom = 0;
  usize invalid_range = 0;
  std::vector<bed::BED> invalid_intervals{};
  invalid_intervals.reserve(10);
  for (const auto &interval : regions_of_interest) {
    const auto match = chroms.find(interval.name);
    if (match == chroms.end()) {
      if (invalid_intervals.size() < invalid_intervals.capacity()) {
        invalid_intervals.push_back(interval);
      }
      ++invalid_chrom;
      continue;
    }
    if (interval.chrom_end > match->size()) {
      if (invalid_intervals.size() < invalid_intervals.capacity()) {
        invalid_intervals.push_back(interval);
      }
      ++invalid_range;
    }
  }

  if (!invalid_intervals.empty()) {
    const auto num_invalid_intervals = invalid_chrom + invalid_range;
    throw std::runtime_error(
        fmt::format("detected {} invalid intervals:\n"
                    " - {} intervals refer to missing chromosome(s)\n"
                    " - {} intervals span outside existing chromosome(s)\n"
                    "Showing {}{} invalid record(s):\n"
                    "- {}",
                    num_invalid_intervals, invalid_chrom, invalid_range,
                    invalid_intervals.size() == num_invalid_intervals ? "" : "the first ",
                    invalid_intervals.size(), fmt::join(invalid_intervals, "\n - ")));
  }
}

static void validate_chrom_sizes(const hictk::Reference &chrom_sizes_cooler,
                                 const hictk::Reference &chrom_sizes) {
  std::vector<std::string> errors{};
  for (const auto &chrom1 : chrom_sizes_cooler) {
    auto match = chrom_sizes.find(chrom1);
    if (match == chrom_sizes.end()) {
      errors.emplace_back(
          fmt::format("{} is present in Cooler files but not in .chrom.sizes", chrom1));
      continue;
    }

    if (match->size() != chrom1.size()) {
      errors.emplace_back(fmt::format(
          "{} has size {} according to Cooler files and {} according to .chrom.sizes file",
          chrom1.name(), chrom1.size(), match->size()));
    }
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format("cooler files and chrom.size()s file contain conflicting information:\n"
                    " - {}",
                    fmt::join(errors, "\n - ")));
  }
}

[[nodiscard]] static hictk::Reference generate_chrom_annotation(
    const hictk::cooler::File &f1, const hictk::cooler::File &f2,
    const std::filesystem::path &path_to_chrom_sizes) {
  // Import chromosomes from Cooler files and select chroms present in both files
  auto chroms = import_chroms_from_coolers(f1, f2);
  if (chroms.empty()) {
    throw std::runtime_error(
        fmt::format("coolers at URI \"{}\" and \"{}\" appear to have no chromosome in common. "
                    "Please make sure both Coolers are using the same genome assembly.\n"
                    "Chromosomes found in Cooler #1:\n - {}\n"
                    "Chromosomes found in Cooler #2:\n - {}",
                    f1.uri(), f2.uri(), fmt::join(f1.chromosomes(), "\n - "),
                    fmt::join(f2.chromosomes(), "\n - ")));
  }

  if (path_to_chrom_sizes.empty()) {
    return chroms;
  }
  if (!path_to_chrom_sizes.empty()) {
    validate_chrom_sizes(chroms, import_chroms_from_chrom_sizes_file(path_to_chrom_sizes));
  }
  return chroms;
}

[[nodiscard]] static std::vector<bed::BED> generate_regions_of_interest_for_eval(
    const hictk::Reference &chroms, const std::filesystem::path &path_to_regions_of_interest) {
  if (path_to_regions_of_interest.empty()) {
    return chromset_to_bed(chroms);
  }

  // Import regions of interest
  auto regions_of_interest = import_regions_of_interest(chroms, path_to_regions_of_interest);
  if (regions_of_interest.empty()) {
    return chromset_to_bed(chroms);
  }

  validate_regions_of_interest(chroms, regions_of_interest);

  return regions_of_interest;
}

[[nodiscard]] static isize find_col_idx(const std::filesystem::path &path_to_weights,
                                        std::string_view header, std::string_view col_name) {
  const auto toks = absl::StrSplit(header, '\t');
  const auto it = std::find(toks.begin(), toks.end(), col_name);
  if (it == toks.end()) {
    throw std::runtime_error(
        fmt::format("unable to find column \"{}\" in the header \"{}\" of file {}", col_name,
                    header, path_to_weights));
  }
  return std::distance(toks.begin(), it);
}

template <class Range>
[[nodiscard]] static isize find_col_idx(const std::filesystem::path &path_to_weights,
                                        std::string_view header, const Range &col_names) {
  assert(col_names.size() > 1);
  for (const auto &name : col_names) {
    try {
      return find_col_idx(path_to_weights, header, name);
    } catch (const std::exception &e) {
      if (!absl::StartsWith(e.what(), "Unable to find column")) {
        throw;
      }
    }
  }
  throw std::runtime_error(fmt::format(
      "Unable to find any columns named like \"{}\" or \"{}\" in header \"{}\" of file {}",
      fmt::join(col_names.begin(), col_names.end() - 1, "\", \""), col_names.back(), header,
      path_to_weights));
}

[[nodiscard]] static phmap::flat_hash_map<std::string, std::vector<double>> import_weights(
    const std::filesystem::path &path_to_weights, const std::string_view weight_column_name,
    const usize nbins, const bool reciprocal_weights) {
  assert(nbins != 0);
  phmap::flat_hash_map<std::string, std::vector<double>> weights;

  if (path_to_weights.empty()) {
    return weights;
  }

  assert(std::filesystem::exists(path_to_weights));
  compressed_io::Reader r(path_to_weights);

  std::string buff;
  r.getline(buff);
  if (buff.empty()) {
    throw std::runtime_error(fmt::format(
        "Failed to read header from file {}: file appears to be empty", path_to_weights));
  }
  using namespace std::string_view_literals;
  const auto col_chrom_idx =
      find_col_idx(path_to_weights, buff, std::array<std::string_view, 2>{"chrom"sv, "region"sv});
  const auto col_diag_idx = find_col_idx(path_to_weights, buff, "diag"sv);
  const auto col_weight_idx = find_col_idx(path_to_weights, buff, weight_column_name);

  for (usize i = 1; r.getline(buff); ++i) {
    const auto toks = absl::StrSplit(buff, '\t');
    if (const auto ntoks = std::distance(toks.begin(), toks.end()); ntoks < col_weight_idx) {
      throw std::runtime_error(fmt::format(
          "Failed to parse column \"{}\" at line #{}: expected at least {} columns, found {}",
          weight_column_name, i, col_weight_idx, ntoks));
    }
    const std::string_view chrom = *std::next(toks.begin(), col_chrom_idx);
    weights.try_emplace(chrom, nbins, 0);  // try_emplace a vector of doubles
    const auto diag = utils::parse_numeric_or_throw<usize>(*std::next(toks.begin(), col_diag_idx));
    if (diag >= nbins) {
      continue;
    }
    const auto weight =
        utils::parse_numeric_or_throw<double>(*std::next(toks.begin(), col_weight_idx));
    weights.at(chrom)[diag] = [&]() {
      if (std::isnan(weight) || weight == 0.0) {
        return 0.0;
      }
      if (reciprocal_weights) {
        return 1.0 / weight;
      }
      return weight;
    }();
  }

  return weights;
}

static void validate_weights(
    const std::vector<bed::BED> &regions_of_interest,
    const phmap::flat_hash_map<std::string, std::vector<double>> &weights) {
  std::vector<std::string> missing_chroms;
  for (const auto &bed : regions_of_interest) {
    if (!weights.contains(bed.chrom)) {
      missing_chroms.push_back(bed.chrom);
    }
  }
  if (!missing_chroms.empty()) {
    throw std::runtime_error(fmt::format("Unable to read weights for the following chromosomes: {}",
                                         fmt::join(missing_chroms, ", ")));
  }
}

[[nodiscard]] static io::bigwig::Writer create_bwig_file(
    const std::vector<std::string> &chrom_names, const std::vector<u32> &chrom_sizes,
    std::string_view base_name, std::string_view suffix) {
  auto bw = io::bigwig::Writer(fmt::format("{}_{}", base_name, absl::StripPrefix(suffix, "_")));
  bw.write_chromosomes(chrom_names, chrom_sizes);
  return bw;
}

enum class StripeDirection : u8f { vertical, horizontal };

template <class N, class = std::enable_if_t<std::is_arithmetic_v<N>>>
static size_t mask_zero_pixels(const std::vector<N> &v1, const std::vector<N> &v2,
                               std::vector<double> &weights) {
  assert(v1.size() == v2.size());
  assert(v1.size() == weights.size());

  usize masked_pixels = 0;
  for (usize i = 0; i < v1.size(); ++i) {
    if (v1[i] == N(0) || v2[i] == N(0)) {
      ++masked_pixels;
      weights[i] = 0;
    }
  }
  return masked_pixels;
}

template <class N>
[[nodiscard]] static auto compute_custom_metric(const std::vector<N> &ref_pixels,
                                                const std::vector<N> &tgt_pixels) {
  assert(ref_pixels.size() == tgt_pixels.size());
  // assert(std::all_of(ref_pixels.begin(), ref_pixels.end(),
  //                    [&](const auto n) { return n == 0 || n == 1; }));
  // assert(std::all_of(tgt_pixels.begin(), tgt_pixels.end(),
  //                    [&](const auto n) { return n == 0 || n == 1; }));

  // Do a backward search for the first non-zero pixel
  auto non_zero_backward_search = [](const auto &vect) -> usize {
    const auto it = std::find_if(vect.rbegin(), vect.rend(), [](const auto n) { return n != 0; });
    if (MODLE_UNLIKELY(it == vect.rend())) {
      return 0;
    }
    return static_cast<usize>(vect.rend() - 1 - it);
  };

  // the two backward searches find the index of the last non-zero pixel in the reference and
  // target vector of pixels respectively. These two indices are used to mark the end of a stripe
  // produced by DNA loop extrusion. NOTE: The curly braces are necessary as we want to call the
  // initializer-list overload to avoid dangling references. See Notes section:
  // https://en.cppreference.com/w/cpp/algorithm/minmax
  const auto [i0, i1] =
      std::minmax({non_zero_backward_search(ref_pixels), non_zero_backward_search(tgt_pixels)});

  usize score = 0;
  for (auto i = i0; i != i1; ++i) {  // Count mismatches of pixels between i0 and i1
    score += ref_pixels[i] != tgt_pixels[i];
  }

  struct CustomMetric {
    double correctly_classified_pixels;
    double incorrectly_classified_pixels;
  };

  return CustomMetric{static_cast<double>(i1 - i0 - score), static_cast<double>(score)};
}

struct MetricsBuff {
  std::vector<double> metric1;
  std::vector<double> metric2;
};

template <StripeDirection stripe_direction, class N, class WeightIt = utils::RepeatIterator<double>>
[[nodiscard]] static auto compute_metric(
    const enum eval_config::Metric metric, const ContactMatrixDense<N> &ref_contacts,
    const ContactMatrixDense<N> &tgt_contacts, const bool mask_zero_pixels_,
    [[maybe_unused]] WeightIt weight_first = utils::RepeatIterator<double>(1)) {
  assert(ref_contacts.nrows() == tgt_contacts.nrows());
  const auto nrows = ref_contacts.nrows();
  const auto ncols = ref_contacts.ncols();

  std::vector<double> ref_pixel_buff(nrows);
  std::vector<double> tgt_pixel_buff(nrows);
  std::vector<double> weight_buff;
  if (!std::is_same_v<WeightIt, utils::RepeatIterator<double>> || mask_zero_pixels_) {
    weight_buff.resize(nrows);
    std::copy_n(weight_first, nrows, weight_buff.begin());
  }

  MetricsBuff score_buff;
  score_buff.metric1.resize(ncols, 0.0);
  score_buff.metric2.resize(ncols, 0.0);

  // The scoring functions used down below either return a single double or a struct with two
  // fields with different names. This lambda uses structure binding to convert structs with two
  // fields into a pair. If val is of fundamental type (e.g. double), then this function returns a
  // pair of val and 0
  auto marshall_metric = [](const auto val) {
    using T = decltype(val);
    if constexpr (std::is_fundamental_v<T>) {
      return std::make_pair(val, static_cast<T>(0));
    } else {
      const auto [a, b] = val;
      return std::make_pair(a, b);
    }
  };

  auto compute_metric_once = [&]() {
    using m = eval_config::Metric;
    switch (metric) {
      case m::custom:
        return marshall_metric(compute_custom_metric(ref_pixel_buff, tgt_pixel_buff));

      case m::eucl_dist:
        return marshall_metric(
            weight_buff.empty()
                ? stats::sed(ref_pixel_buff.begin(), ref_pixel_buff.end(), tgt_pixel_buff.begin())
                : stats::weighted_sed(ref_pixel_buff.begin(), ref_pixel_buff.end(),
                                      tgt_pixel_buff.begin(), weight_buff.begin()));

      case m::pearson:
        return marshall_metric(
            weight_buff.empty() ? stats::Pearson<double>{}(ref_pixel_buff, tgt_pixel_buff) :

                                stats::Pearson<double>{}(ref_pixel_buff, tgt_pixel_buff,
                                                         weight_buff));

      case m::rmse:
        return marshall_metric(
            weight_buff.empty()
                ? stats::rmse(ref_pixel_buff.begin(), ref_pixel_buff.end(), tgt_pixel_buff.begin())
                : stats::weighted_rmse(ref_pixel_buff.begin(), ref_pixel_buff.end(),
                                       tgt_pixel_buff.begin(), weight_buff.begin()));

      case m::spearman:
        return marshall_metric(
            weight_buff.empty()
                ? stats::Spearman<double>{nrows}(ref_pixel_buff, tgt_pixel_buff)
                : stats::Spearman<double>{nrows}(ref_pixel_buff, tgt_pixel_buff, weight_buff));
    }
    MODLE_UNREACHABLE_CODE;
  };

  for (usize i = 0; i < ncols; ++i) {
    using d = StripeDirection;
    if constexpr (stripe_direction == d::vertical) {
      ref_contacts.unsafe_get_column(i, ref_pixel_buff);
      tgt_contacts.unsafe_get_column(i, tgt_pixel_buff);
    } else {
      static_assert(stripe_direction == d::horizontal);
      ref_contacts.unsafe_get_row(i, ref_pixel_buff);
      tgt_contacts.unsafe_get_row(i, tgt_pixel_buff);
    }
    ref_pixel_buff.resize(nrows);
    tgt_pixel_buff.resize(nrows);

    if (mask_zero_pixels_) {
      std::copy_n(weight_first, nrows, weight_buff.begin());
      mask_zero_pixels(ref_pixel_buff, tgt_pixel_buff, weight_buff);
    }

    const auto [s1, s2] = compute_metric_once();
    score_buff.metric1[i] = s1;
    score_buff.metric2[i] = s2;
  }
  return score_buff;
}

template <StripeDirection stripe_direction, class N>
[[nodiscard]] static auto compute_metric(const enum eval_config::Metric metric,
                                         const ContactMatrixDense<N> &ref_contacts,
                                         const ContactMatrixDense<N> &tgt_contacts,
                                         const bool mask_zero_pixels_,
                                         const std::vector<double> &weights) {
  assert(ref_contacts.nrows() == tgt_contacts.nrows());
  assert(weights.empty() || ref_contacts.nrows() == weights.size());
  return compute_metric<stripe_direction>(metric, ref_contacts, tgt_contacts, mask_zero_pixels_,
                                          weights.begin());
}

[[nodiscard]] [[maybe_unused]] static constexpr std::string_view corr_method_to_str(
    const enum eval_config::Metric metric, bool prettify = false) {
  using m = eval_config::Metric;
  if (metric == m::custom) {
    return prettify ? "Custom metric" : "custom_metric";
  }
  if (metric == m::pearson) {
    return prettify ? "Pearson correlation" : "pearson_corr";
  }
  if (metric == m::spearman) {
    return prettify ? "Spearman correlation" : "spearman_corr";
  }

  if (metric == m::eucl_dist) {
    return prettify ? "Euclidean distance" : "eucl_dist";
  }

  if (metric == m::rmse) {
    return prettify ? "RMSE" : "rmse";
  }
  std::abort();
}

[[nodiscard]] [[maybe_unused]] static std::string_view direction_to_str(const StripeDirection m) {
  if (m == StripeDirection::vertical) {
    return "Vertical";
  }
  if (m == StripeDirection::horizontal) {
    return "Horizontal";
  }
  MODLE_UNREACHABLE_CODE;
}

[[nodiscard]] static auto init_writers(const modle::tools::eval_config &c,
                                       const hictk::Reference &chroms, const bool weighted) {
  std::vector<std::string> chrom_names(utils::conditional_static_cast<usize>(chroms.size()));
  std::vector<u32> chrom_sizes(utils::conditional_static_cast<usize>(chroms.size()));

  std::transform(chroms.begin(), chroms.end(), chrom_names.begin(),
                 [](const auto &chrom) { return chrom.name(); });

  std::transform(chroms.begin(), chroms.end(), chrom_sizes.begin(),
                 [](const auto &chrom) { return static_cast<u32>(chrom.size()); });

  const auto name = corr_method_to_str(c.metric);
  const auto bname1 = fmt::format("{}{}_vertical", name, weighted ? "_weighted" : "");
  const auto bname2 = fmt::format("{}{}_horizontal", name, weighted ? "_weighted" : "");

  struct Writer {
    struct InternalWriterPair {
      std::unique_ptr<io::bigwig::Writer> bwig{nullptr};
      std::unique_ptr<compressed_io::Writer> tsv_gz{nullptr};
    };
    InternalWriterPair horizontal;
    InternalWriterPair vertical;
  };

  Writer writers;
  writers.vertical = Writer::InternalWriterPair{
      // clang-format off
          std::make_unique<io::bigwig::Writer>(create_bwig_file(chrom_names, chrom_sizes, c.output_prefix.string(), fmt::format("{}.bw", bname1))),
          std::make_unique<compressed_io::Writer>(fmt::format("{}_{}.tsv.gz", c.output_prefix.string(), bname1))
      // clang-format on
  };
  writers.horizontal = Writer::InternalWriterPair{
      // clang-format off
          std::make_unique<io::bigwig::Writer>(create_bwig_file(chrom_names, chrom_sizes, c.output_prefix.string(), fmt::format("{}.bw", bname2))),
          std::make_unique<compressed_io::Writer>(fmt::format("{}_{}.tsv.gz", c.output_prefix.string(), bname2))
      // clang-format on
  };

  switch (c.metric) {
    case eval_config::Metric::custom: {
      constexpr std::string_view header =
          "chrom\tchrom_start\tchrom_end\tcorrectly_classified_pixels\tincorrectly_classified_"
          "pixels\n";
      writers.horizontal.tsv_gz->write(header);
      writers.vertical.tsv_gz->write(header);
      break;
    }
    case eval_config::Metric::eucl_dist:
      [[fallthrough]];
    case eval_config::Metric::rmse: {
      constexpr std::string_view header = "chrom\tchrom_start\tchrom_end\tscore\n";
      writers.horizontal.tsv_gz->write(header);
      writers.vertical.tsv_gz->write(header);
      break;
    }
    case eval_config::Metric::spearman:
      [[fallthrough]];
    case eval_config::Metric::pearson: {
      constexpr std::string_view header =
          "chrom\tchrom_start\tchrom_end\tcorrelation\tsignificance\n";
      writers.horizontal.tsv_gz->write(header);
      writers.vertical.tsv_gz->write(header);
    }
  }

  return writers;
}

[[nodiscard]] static std::string format_tsv_record(const enum eval_config::Metric metric,
                                                   const std::string_view chrom_name_,
                                                   const bp_t start_pos_, const bp_t end_pos_,
                                                   const double score1, const double score2) {
  if (metric == eval_config::Metric::eucl_dist || metric == eval_config::Metric::rmse) {
    assert(score2 == 0.0);
    return fmt::format(FMT_COMPILE("{}\t{}\t{}\t{}\n"), chrom_name_, start_pos_, end_pos_, score1);
  }
  return fmt::format(FMT_COMPILE("{}\t{}\t{}\t{}\t{}\n"), chrom_name_, start_pos_, end_pos_, score1,
                     score2);
}

template <StripeDirection stripe_direction, class Writers, class N>
static void run_task(const enum eval_config::Metric metric, const bed::BED &interval,
                     Writers &writers, const ContactMatrixDense<N> &ref_contacts,
                     const ContactMatrixDense<N> &tgt_contacts, const bool exclude_zero_pixels,
                     const usize bin_size,
                     const phmap::flat_hash_map<std::string, std::vector<double>> &weights) {
  try {
    auto t0 = absl::Now();
    const auto metrics =
        weights.empty()
            ? compute_metric<stripe_direction>(metric, ref_contacts, tgt_contacts,
                                               exclude_zero_pixels)
            : compute_metric<stripe_direction>(metric, ref_contacts, tgt_contacts,
                                               exclude_zero_pixels, weights.at(interval.chrom));

    SPDLOG_INFO("{} for {} stripes from interval {}:{}-{} computed in {}.",
                corr_method_to_str(metric, true),
                absl::AsciiStrToLower(direction_to_str(stripe_direction)), interval.chrom,
                interval.chrom_start, interval.chrom_end, absl::FormatDuration(absl::Now() - t0));
    t0 = absl::Now();
    auto &writer =
        stripe_direction == StripeDirection::horizontal ? writers.horizontal : writers.vertical;
    writer.bwig->write_range(interval.chrom, absl::MakeSpan(metrics.metric1), bin_size, bin_size,
                             interval.chrom_start);

    std::string buff;
    for (usize i = 0; i < metrics.metric1.size(); ++i) {
      const auto start_pos = interval.chrom_start + (bin_size * i);
      const auto end_pos = std::min(start_pos + bin_size, interval.chrom_end);
      buff = format_tsv_record(metric, interval.chrom, start_pos, end_pos, metrics.metric1[i],
                               metrics.metric2[i]);
      writer.tsv_gz->write(buff);
    }
    SPDLOG_INFO("{} values have been written to files \"{}.{{tsv.gz,bw}}\" in {}.",
                metrics.metric1.size(), writer.bwig->path().stem().string(),
                absl::FormatDuration(absl::Now() - t0));
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        "The following error occurred while computing {} on {} stripes for {}:{}-{}: {}",
        corr_method_to_str(metric, true), absl::AsciiStrToLower(direction_to_str(stripe_direction)),
        interval.chrom, interval.chrom_start, interval.chrom_end, e.what()));
  }
}

static void log_regions_for_evaluation(const std::vector<bed::BED> &intervals) {
  if (intervals.size() == 1) {
    const auto &chrom = intervals.front();
    SPDLOG_INFO("Computing metric(s) for interval {}:{}-{}", chrom.chrom, chrom.chrom_start,
                chrom.chrom_end);
    return;
  }

  std::vector<std::string> printable_intervals(
      utils::conditional_static_cast<usize>(intervals.size()));

  std::transform(
      intervals.begin(), intervals.end(), printable_intervals.begin(), [](const auto &interval) {
        return fmt::format("{}:{}-{}", interval.chrom, interval.chrom_start, interval.chrom_end);
      });
  SPDLOG_INFO("Computing metric(s) for the following {} intervals:\n - {}", intervals.size(),
              fmt::join(printable_intervals, "\n - "));
}

void eval_subcmd(const modle::tools::eval_config &c) {
  const auto t0 = absl::Now();
  hictk::cooler::File ref_cooler(c.reference_cooler_uri.string());
  hictk::cooler::File tgt_cooler(c.input_cooler_uri.string());

  assert(ref_cooler.resolution() == tgt_cooler.resolution());

  const auto bin_size = ref_cooler.resolution();

  const auto chroms = generate_chrom_annotation(ref_cooler, tgt_cooler, c.path_to_chrom_sizes);

  const auto intervals =
      generate_regions_of_interest_for_eval(chroms, c.path_to_regions_of_interest_bed);

  // Import weights
  const auto weights =
      import_weights(c.path_to_weights, c.weight_column_name,
                     (c.diagonal_width + bin_size - 1) / bin_size, c.reciprocal_weights);

  if (!weights.empty()) {
    validate_weights(intervals, weights);
  }

  log_regions_for_evaluation(intervals);
  if (const auto &output_dir = c.output_prefix.parent_path(); !output_dir.empty()) {
    std::filesystem::create_directories(output_dir.string());
  }

  auto writers = init_writers(c, chroms, !weights.empty());

  BS::light_thread_pool tpool(c.nthreads);

  std::array<std::future<void>, 2> return_codes;

  for (const auto &interval : intervals) {
    const auto t1 = absl::Now();
    const auto coord_str = [&]() {
      if (interval.chrom_start == 0) {
        return fmt::format("{}:{}", interval.chrom, interval.chrom_end);
      }
      return fmt::format("{}:{}-{}", interval.chrom, interval.chrom_start, interval.chrom_end);
    }();

    if (io::query_returns_no_pixels(ref_cooler, interval.chrom, interval.chrom_start,
                                    interval.chrom_end) &&
        io::query_returns_no_pixels(tgt_cooler, interval.chrom, interval.chrom_start,
                                    interval.chrom_end)) {
      SPDLOG_WARN("Read 0 contacts for {}. SKIPPING!", coord_str);
      continue;
    }

    SPDLOG_INFO("Reading contacts for {}...", coord_str);
    auto ref_matrix = io::read_contact_matrix_from_cooler<double>(
        ref_cooler, interval.chrom, interval.chrom_start, interval.chrom_end,
        static_cast<bp_t>(c.diagonal_width));
    auto tgt_matrix = io::read_contact_matrix_from_cooler<double>(
        tgt_cooler, interval.chrom, interval.chrom_start, interval.chrom_end,
        static_cast<bp_t>(c.diagonal_width));
    SPDLOG_INFO("Read {} contacts for {} in {}",
                ref_matrix.get_tot_contacts() + tgt_matrix.get_tot_contacts(), coord_str,
                absl::FormatDuration(absl::Now() - t1));

    // Normalize contact matrix before computing the correlation/distance metrics
    if (c.normalize) {
      const auto t00 = absl::Now();
      SPDLOG_INFO("Normalizing contact matrices for {}...", coord_str);
      return_codes[0] = tpool.submit_task([&]() { ref_matrix.normalize_inplace(); });
      return_codes[1] = tpool.submit_task([&]() { tgt_matrix.normalize_inplace(); });
      tpool.wait();
      try {
        // Handle exceptions thrown inside worker threads
        std::ignore = return_codes[0];
        std::ignore = return_codes[1];
      } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format(
            "The following error occurred while normalizing contact matrices for {}: {}", coord_str,
            e.what()));
      }
      SPDLOG_INFO("DONE! Normalization took {}.", absl::FormatDuration(absl::Now() - t00));
    }

    using d = StripeDirection;
    return_codes[0] = tpool.submit_task([&, interval = interval]() {
      run_task<d::horizontal>(c.metric, interval, writers, ref_matrix, tgt_matrix,
                              c.exclude_zero_pxls, bin_size, weights);
    });
    return_codes[1] = tpool.submit_task([&, interval = interval]() {
      run_task<d::vertical>(c.metric, interval, writers, ref_matrix, tgt_matrix,
                            c.exclude_zero_pxls, bin_size, weights);
    });
    tpool.wait();
    // Raise exceptions thrown inside worker threads (if any)
    std::ignore = return_codes[0];
    std::ignore = return_codes[1];
  }
  SPDLOG_INFO("DONE in {}!", absl::FormatDuration(absl::Now() - t0));
}

}  // namespace modle::tools
