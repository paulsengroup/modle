// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/algorithm/container.h>           // for c_set_intersection
#include <absl/container/btree_map.h>           // for btree_map, btree_iterator, map_params<>::...
#include <absl/container/flat_hash_map.h>       // for flat_hash_map
#include <absl/strings/ascii.h>                 // AsciiStrToLower
#include <absl/strings/str_replace.h>           // for ReplaceAll
#include <absl/strings/strip.h>                 // for StripPrefix
#include <absl/time/clock.h>                    // for Now
#include <absl/time/time.h>                     // for FormatDuration, operator-, Time
#include <absl/types/span.h>                    // for MakeSpan
#include <cpp-sort/comparators/natural_less.h>  // for natural_less_t
#include <fmt/compile.h>                        // for FMT_COMPILE
#include <fmt/format.h>                         // for format, make_format_args, vformat_to, FMT...
#include <spdlog/spdlog.h>                      // for info

#include <BS_thread_pool.hpp>  // for BS::thread_pool
#include <algorithm>           // for all_of, find_if, transform, minmax, max
#include <cassert>             // for assert
#include <cstdio>              // for stderr
#include <exception>           // for exception
#include <filesystem>          // for operator<<, path
#include <future>              // for future
#include <iosfwd>              // for streamsize
#include <iterator>            // for insert_iterator, inserter
#include <memory>              // for unique_ptr, shared_ptr, __shared_ptr_access
#include <stdexcept>           // for runtime_error, overflow_error
#include <string>              // for string, basic_string
#include <string_view>         // for string_view
#include <utility>             // for tuple_element<>::type, pair, make_pair
#include <vector>              // for vector

#include "modle/bed/bed.hpp"        // for BED_tree, BED_tree<>::value_type, Parser
#include "modle/bigwig/bigwig.hpp"  // for Writer
#include "modle/common/common.hpp"  // for u32, usize, bp_t, u8, i64
#include "modle/common/fmt_helpers.hpp"
#include "modle/common/utils.hpp"              // for identity::operator()
#include "modle/contact_matrix_dense.hpp"      // for ContactMatrixDense
#include "modle/cooler/cooler.hpp"             // for Cooler, Cooler::READ_ONLY
#include "modle/interval_tree.hpp"             // for IITree, IITree::IITree<I, T>, IITree::empty
#include "modle/stats/correlation.hpp"         // for Pearson, Spearman
#include "modle_tools/modle_tools_config.hpp"  // for eval_config
#include "modle_tools/tools.hpp"               // for eval_subcmd

namespace modle::tools {

using ChromSet = absl::btree_map<std::string, std::pair<bp_t, bp_t>, cppsort::natural_less_t>;

[[nodiscard]] static ChromSet import_chrom_subranges(
    const std::filesystem::path &path_to_chrom_subranges) {
  if (path_to_chrom_subranges.empty()) {
    return ChromSet{};
  }

  ChromSet chroms;
  const auto records = bed::Parser(path_to_chrom_subranges).parse_all();
  std::transform(records.begin(), records.end(), std::inserter(chroms, chroms.end()),
                 [](const auto &r) -> std::pair<std::string, std::pair<bp_t, bp_t>> {
                   return {r.chrom, {r.chrom_start, r.chrom_end}};
                 });
  return chroms;
}

template <class N>
[[nodiscard]] static ChromSet import_chroms_from_cool(cooler::Cooler<N> &cooler_file) {
  try {
    ChromSet chrom_set;
    auto chroms = cooler_file.get_chroms();

    std::transform(chroms.begin(), chroms.end(), std::inserter(chrom_set, chrom_set.begin()),
                   [](const auto &chrom) {
                     return std::make_pair(chrom.first,
                                           std::make_pair(i64(0), static_cast<i64>(chrom.second)));
                   });
    return chrom_set;

  } catch (const std::runtime_error &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("An error occurred while reading file {}: {}"),
                                         cooler_file.get_path(), e.what()));
  }
}

template <class N>
[[nodiscard]] static ChromSet select_chromosomes_for_eval(
    cooler::Cooler<N> &c1, cooler::Cooler<N> &c2,
    const std::filesystem::path &path_to_chrom_subranges) {
  // Import chromosomes from Cooler files and select shared chromosomes
  const auto chrom_set1 = import_chroms_from_cool(c1);
  const auto chrom_set2 = import_chroms_from_cool(c2);

  ChromSet chrom_intersection;
  absl::c_set_intersection(chrom_set1, chrom_set2,
                           std::inserter(chrom_intersection, chrom_intersection.begin()),
                           [&](const auto &chrom1, const auto &chrom2) {
                             return cppsort::natural_less(chrom1.first, chrom2.first);
                           });

  if (chrom_intersection.empty()) {
    std::vector<std::string> chrom_printable1(static_cast<usize>(chrom_set1.size()));
    std::vector<std::string> chrom_printable2(static_cast<usize>(chrom_set2.size()));

    std::transform(
        chrom_set1.begin(), chrom_set1.end(), chrom_printable1.begin(), [&](const auto &c) {
          return fmt::format(FMT_STRING("{}:{}-{}"), c.first, c.second.first, c.second.second);
        });
    std::transform(
        chrom_set2.begin(), chrom_set2.end(), chrom_printable2.begin(), [&](const auto &c) {
          return fmt::format(FMT_STRING("{}:{}-{}"), c.first, c.second.first, c.second.second);
        });
    throw std::runtime_error(
        fmt::format(FMT_STRING("Files {} and {} have 0 chromosomes in common. "
                               "Please make sure both files are using the same genome assembly as "
                               "reference.\nChromosomes found in file #1:\n - {}\nChromosomes "
                               "found in file #2:\n - {}"),
                    c1.get_path(), c2.get_path(), fmt::join(chrom_printable1, "\n - "),
                    fmt::join(chrom_printable2, "\n - ")));
  }

  // Look for conflicting information in Cooler files
  std::string err;
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_UNUSED_VARIABLE
  for (const auto &[name, _] : chrom_intersection) {
    if (chrom_set1.at(name).second != chrom_set2.at(name).second) {
      err += fmt::format(FMT_STRING("\n - \"{}\" has size {} in file #1 and {} in file #2"), name,
                         chrom_set1.at(name).second, chrom_set2.at(name).second);
    }
  }
  DISABLE_WARNING_POP

  if (!err.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("File #1 {} and file #2 {} contain conflicting information:{}"),
                    c1.get_path(), c2.get_path(), err));
  }

  // Import chromosome subranges
  auto chrom_subranges = import_chrom_subranges(path_to_chrom_subranges);
  if (chrom_subranges.empty()) {
    return chrom_intersection;
  }

  // Look for conflicting information between Cooler files and BED file
  for (const auto &[name, range] : chrom_subranges) {
    const auto it = chrom_intersection.find(name);
    if (it == chrom_intersection.end()) {
      err += fmt::format(
          FMT_STRING("\n - \"{}\" is present in BED file but missing in one or both Cooler files"),
          name);
      continue;
    }
    const auto &chrom_size = it->second.second;
    const auto &[subrange_start, subrange_end] = range;
    if (subrange_end > chrom_size) {
      err += fmt::format(FMT_STRING("\n - Subrange for \"{}-{}:{}\" lies outside of the entire "
                                    "chromosome (\"{}-{}:{}\")"),
                         name, subrange_start, subrange_end, name, 0, chrom_size);
    }
  }
  if (!err.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Files {}, {} and {} contain conflicting information:{}"),
                    c1.get_path(), c2.get_path(), path_to_chrom_subranges, err));
  }

  return chrom_subranges;
}

[[nodiscard]] static isize find_col_idx(const std::filesystem::path &path_to_weights,
                                        std::string_view header, std::string_view col_name) {
  const auto toks = absl::StrSplit(header, '\t');
  const auto it = std::find(toks.begin(), toks.end(), col_name);
  if (it == toks.end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find column \"{}\" in the header \"{}\" of file {}"),
                    col_name, header, path_to_weights));
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
      FMT_STRING(
          "Unable to find any columns named like \"{}\" or \"{}\" in header \"{}\" of file {}"),
      fmt::join(col_names.begin(), col_names.end() - 1, "\", \""), col_names.back(), header,
      path_to_weights));
}

[[nodiscard]] static absl::flat_hash_map<std::string, std::vector<double>> import_weights(
    const std::filesystem::path &path_to_weights, const std::string_view weight_column_name,
    const usize nbins, const bool reciprocal_weights) {
  assert(nbins != 0);
  absl::flat_hash_map<std::string, std::vector<double>> weights;

  if (path_to_weights.empty()) {
    return weights;
  }

  assert(std::filesystem::exists(path_to_weights));
  compressed_io::Reader r(path_to_weights);

  std::string buff;
  r.getline(buff);
  if (buff.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Falied to read header from file {}: file appears to be empty"),
                    path_to_weights));
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
          FMT_STRING(
              "Failed to parse column \"{}\" at line #{}: expected at least {} columns, found {}"),
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

[[nodiscard]] static io::bigwig::Writer create_bwig_file(
    const std::vector<std::string> &chrom_names, const std::vector<u32> &chrom_sizes,
    std::string_view base_name, std::string_view suffix) {
  auto bw = io::bigwig::Writer(
      fmt::format(FMT_STRING("{}_{}"), base_name, absl::StripPrefix(suffix, "_")));
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

  // the two backward searches find the index of the last non-zero pixel in the reference and target
  // vector of pixels respectively. These two indices are used to mark the end of a stripe produced
  // by DNA loop extrusion.
  // NOTE: The curly braces are necessary as we want to call the initializer-list overload to avoid
  // dangling references. See Notes section: https://en.cppreference.com/w/cpp/algorithm/minmax
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
    const enum eval_config::Metric metric, std::shared_ptr<ContactMatrixDense<N>> ref_contacts,
    std::shared_ptr<ContactMatrixDense<N>> tgt_contacts, const bool mask_zero_pixels_,
    [[maybe_unused]] WeightIt weight_first = utils::RepeatIterator<double>(1)) {
  assert(ref_contacts->nrows() == tgt_contacts->nrows());
  const auto nrows = ref_contacts->nrows();
  const auto ncols = ref_contacts->ncols();

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

  // The scoring functions used down below either return a single double or a struct with two fields
  // with different names. This lambda uses structure binding to convert structs with two fields
  // into a pair.
  // If val is of fundamental type (e.g. double), then this function returns a pair of val and 0
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
      ref_contacts->unsafe_get_column(i, ref_pixel_buff);
      tgt_contacts->unsafe_get_column(i, tgt_pixel_buff);
    } else {
      static_assert(stripe_direction == d::horizontal);
      ref_contacts->unsafe_get_row(i, ref_pixel_buff);
      tgt_contacts->unsafe_get_row(i, tgt_pixel_buff);
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
                                         std::shared_ptr<ContactMatrixDense<N>> ref_contacts,
                                         std::shared_ptr<ContactMatrixDense<N>> tgt_contacts,
                                         const bool mask_zero_pixels_,
                                         const std::vector<double> &weights) {
  assert(ref_contacts->nrows() == tgt_contacts->nrows());
  assert(weights.empty() || ref_contacts->nrows() == weights.size());
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

[[nodiscard]] [[maybe_unused]] static constexpr std::string_view direction_to_str(
    const StripeDirection m) {
  if (m == StripeDirection::vertical) {
    return "Vertical";
  }
  if (m == StripeDirection::horizontal) {
    return "Horizontal";
  }
  std::abort();
}

[[nodiscard]] static auto init_writers(const modle::tools::eval_config &c, const ChromSet &chroms,
                                       const bool weighted) {
  std::vector<std::string> chrom_names(static_cast<usize>(chroms.size()));
  std::vector<u32> chrom_sizes(static_cast<usize>(chroms.size()));

  std::transform(chroms.begin(), chroms.end(), chrom_names.begin(),
                 [](const auto &chrom) { return chrom.first; });

  std::transform(chroms.begin(), chroms.end(), chrom_sizes.begin(),
                 [](const auto &chrom) { return static_cast<u32>(chrom.second.second); });

  const auto name = corr_method_to_str(c.metric);
  const auto bname1 = fmt::format(FMT_STRING("{}{}_vertical"), name, weighted ? "_weighted" : "");
  const auto bname2 = fmt::format(FMT_STRING("{}{}_horizontal"), name, weighted ? "_weighted" : "");

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
          std::make_unique<io::bigwig::Writer>(create_bwig_file(chrom_names, chrom_sizes, c.output_prefix.string(), fmt::format(FMT_STRING("{}.bw"), bname1))),
          std::make_unique<compressed_io::Writer>(fmt::format(FMT_STRING("{}_{}.tsv.gz"), c.output_prefix.string(), bname1))
      // clang-format on
  };
  writers.horizontal = Writer::InternalWriterPair{
      // clang-format off
          std::make_unique<io::bigwig::Writer>(create_bwig_file(chrom_names, chrom_sizes, c.output_prefix.string(), fmt::format(FMT_STRING("{}.bw"), bname2))),
          std::make_unique<compressed_io::Writer>(fmt::format(FMT_STRING("{}_{}.tsv.gz"), c.output_prefix.string(), bname2))
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

template <StripeDirection stripe_direction, class Writers, class N>
static void run_task(const enum eval_config::Metric metric, const std::string &chrom_name,
                     const std::pair<usize, usize> chrom_range, Writers &writers,
                     std::shared_ptr<ContactMatrixDense<N>> &ref_contacts,
                     std::shared_ptr<ContactMatrixDense<N>> &tgt_contacts,
                     const bool exclude_zero_pixels, const usize bin_size,
                     const absl::flat_hash_map<std::string, std::vector<double>> &weights) {
  try {
    auto t0 = absl::Now();
    const auto metrics =
        weights.empty()
            ? compute_metric<stripe_direction>(metric, ref_contacts, tgt_contacts,
                                               exclude_zero_pixels)
            : compute_metric<stripe_direction>(metric, ref_contacts, tgt_contacts,
                                               exclude_zero_pixels, weights.at(chrom_name));

    spdlog::info(FMT_STRING("{} for {} stripes from chrom {} computed in {}."),
                 corr_method_to_str(metric, true),
                 absl::AsciiStrToLower(direction_to_str(stripe_direction)), chrom_name,
                 absl::FormatDuration(absl::Now() - t0));
    t0 = absl::Now();
    const auto [chrom_start, chrom_end] = chrom_range;
    auto &writer =
        stripe_direction == StripeDirection::horizontal ? writers.horizontal : writers.vertical;
    writer.bwig->write_range(chrom_name, absl::MakeSpan(metrics.metric1), bin_size, bin_size,
                             chrom_start);

    std::string buff;
    auto format_tsv_record = [&](const std::string_view chrom_name_, const auto start_pos_,
                                 const auto end_pos_, const auto score1, const auto score2) {
      switch (metric) {
        case eval_config::Metric::eucl_dist:
          [[fallthrough]];
        case eval_config::Metric::rmse:
          assert(score2 == 0.0);
          return fmt::format(FMT_COMPILE("{}\t{}\t{}\t{}\n"), chrom_name_, start_pos_, end_pos_,
                             score1);
        case eval_config::Metric::spearman:
          [[fallthrough]];
        case eval_config::Metric::pearson:
          [[fallthrough]];
        case eval_config::Metric::custom: {
          return fmt::format(FMT_COMPILE("{}\t{}\t{}\t{}\t{}\n"), chrom_name_, start_pos_, end_pos_,
                             score1, score2);
        }
      }
      MODLE_UNREACHABLE_CODE;
    };

    for (usize i = 0; i < metrics.metric1.size(); ++i) {
      const auto start_pos = chrom_start + (bin_size * i);
      const auto end_pos = std::min(start_pos + bin_size, chrom_end);
      buff =
          format_tsv_record(chrom_name, start_pos, end_pos, metrics.metric1[i], metrics.metric2[i]);
      writer.tsv_gz->write(buff);
    }
    spdlog::info(FMT_STRING("{} values have been written to files \"{}.{{tsv.gz,bw}}\" in {}."),
                 metrics.metric1.size(), writer.bwig->path().stem().string(),
                 absl::FormatDuration(absl::Now() - t0));
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("The following error occurred while computing {} on {} stripes for {}: {}"),
        corr_method_to_str(metric, true), absl::AsciiStrToLower(direction_to_str(stripe_direction)),
        chrom_name, e.what()));
  }
}

void eval_subcmd(const modle::tools::eval_config &c) {
  auto ref_cooler = cooler::Cooler<double>(c.path_to_reference_matrix,
                                           cooler::Cooler<double>::IO_MODE::READ_ONLY, c.bin_size);
  auto tgt_cooler = cooler::Cooler<double>(c.path_to_input_matrix,
                                           cooler::Cooler<double>::IO_MODE::READ_ONLY, c.bin_size);

  if (c.bin_size == 0) {  // Bin size was not specified
    const auto &ref_bin_size = ref_cooler.get_bin_size();
    const auto &tgt_bin_size = tgt_cooler.get_bin_size();
    if (ref_bin_size != tgt_bin_size) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Reference and input matrices appear to have different bin sizes:\n"
                     " - {}: {} bp\n"
                     " - {}: {} bp"),
          c.path_to_reference_matrix, ref_bin_size, c.path_to_input_matrix, tgt_bin_size));
    }
  }

  const auto bin_size = ref_cooler.get_bin_size();

  // This cannot be made const
  auto chromosomes = select_chromosomes_for_eval(ref_cooler, tgt_cooler, c.path_to_chrom_subranges);

  const auto weights =
      import_weights(c.path_to_weights, c.weight_column_name,
                     (c.diagonal_width + bin_size - 1) / bin_size, c.reciprocal_weights);

  // Import weights
  if (!weights.empty()) {
    std::vector<std::string> missing_chroms;
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_UNUSED_VARIABLE
    for (const auto &[chrom, _] : chromosomes) {
      if (!weights.contains(chrom)) {
        missing_chroms.push_back(chrom);
      }
    }
    DISABLE_WARNING_POP
    if (!missing_chroms.empty()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Failed to import weights from file {} for the following chromosomes: {}"),
          c.path_to_weights, fmt::join(missing_chroms, ", ")));
    }
  }

  if (chromosomes.size() == 1) {
    spdlog::info(FMT_STRING("Computing metric(s) for chromosome: \"{}\""),
                 chromosomes.begin()->first);
  } else {
    std::vector<std::string> chrom_names(static_cast<usize>(chromosomes.size()));
    std::transform(chromosomes.begin(), chromosomes.end(), chrom_names.begin(),
                   [](const auto &chrom) { return chrom.first; });
    spdlog::info(FMT_STRING("Computing metric(s) for the following {} chromosomes: \"{}\""),
                 chromosomes.size(), fmt::join(chrom_names, "\", \""));
  }

  if (const auto &output_dir = c.output_prefix.parent_path(); !output_dir.empty()) {
    std::filesystem::create_directories(output_dir.string());
  }

  auto writers = init_writers(c, chromosomes, !weights.empty());

  const auto nrows = (c.diagonal_width + bin_size - 1) / bin_size;
  assert(nrows != 0);

  const auto t0 = absl::Now();
  BS::thread_pool tpool(static_cast<u32>(c.nthreads));

  std::array<std::future<void>, 2> return_codes;

  for (const auto &[chrom_name_sv, chrom_range] : chromosomes) {
    const std::string chrom_name{chrom_name_sv};

    const auto t1 = absl::Now();
    spdlog::info(FMT_STRING("Reading contacts for {}..."), chrom_name_sv);

    auto ref_contacts = std::make_shared<ContactMatrixDense<double>>(
        ref_cooler.cooler_to_cmatrix(chrom_name, nrows, chrom_range));
    auto tgt_contacts = std::make_shared<ContactMatrixDense<double>>(
        tgt_cooler.cooler_to_cmatrix(chrom_name, nrows, chrom_range));
    spdlog::info(FMT_STRING("Read {} contacts for {} in {}"),
                 ref_contacts->get_tot_contacts() + tgt_contacts->get_tot_contacts(), chrom_name_sv,
                 absl::FormatDuration(absl::Now() - t1));

    // Normalize contact matrix before computing the correlation/distance metrics
    if (c.normalize) {
      const auto t00 = absl::Now();
      spdlog::info(FMT_STRING("Normalizing contact matrices for {}..."), chrom_name_sv);
      return_codes[0] = tpool.submit([&]() { ref_contacts->normalize_inplace(); });
      return_codes[1] = tpool.submit([&]() { tgt_contacts->normalize_inplace(); });
      tpool.wait_for_tasks();
      try {
        // Handle exceptions thrown inside worker threads
        std::ignore = return_codes[0];
        std::ignore = return_codes[1];
      } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format(
            FMT_STRING(
                "The following error occurred while normalizing contact matrices for {}: {}"),
            chrom_name_sv, e.what()));
      }
      spdlog::info(FMT_STRING("DONE! Normalization took {}."),
                   absl::FormatDuration(absl::Now() - t00));
    }

    using d = StripeDirection;
    return_codes[0] = tpool.submit([&, chrom_range = chrom_range]() {
      run_task<d::horizontal>(c.metric, chrom_name, chrom_range, writers, ref_contacts,
                              tgt_contacts, c.exclude_zero_pxls, bin_size, weights);
    });
    return_codes[1] = tpool.submit([&, chrom_range = chrom_range]() {
      run_task<d::vertical>(c.metric, chrom_name, chrom_range, writers, ref_contacts, tgt_contacts,
                            c.exclude_zero_pxls, bin_size, weights);
    });
    tpool.wait_for_tasks();
    // Raise exceptions thrown inside worker threads (if any)
    std::ignore = return_codes[0];
    std::ignore = return_codes[1];
  }
  spdlog::info(FMT_STRING("DONE in {}!"), absl::FormatDuration(absl::Now() - t0));
}

}  // namespace modle::tools
