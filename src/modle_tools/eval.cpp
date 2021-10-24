// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/algorithm/container.h>           // for c_set_intersection
#include <absl/container/btree_map.h>           // for btree_map, btree_iterator, map_params<>::...
#include <absl/strings/strip.h>                 // for StripPrefix
#include <absl/time/clock.h>                    // for Now
#include <absl/time/time.h>                     // for FormatDuration, operator-, Time
#include <absl/types/span.h>                    // for MakeSpan
#include <cpp-sort/comparators/natural_less.h>  // for natural_less_t
#include <fmt/format.h>                         // for format, make_format_args, vformat_to, FMT...
#include <fmt/ostream.h>                        // for formatbuf<>::int_type
#include <spdlog/spdlog.h>                      // for info

#include <algorithm>                        // for transform, max
#include <boost/filesystem/operations.hpp>  // for create_directories
#include <boost/filesystem/path.hpp>        // for operator<<, path
#include <cassert>                          // for assert
#include <cstdint>                          // for uint_fast8_t
#include <cstdio>                           // for stderr
#include <exception>                        // for exception
#include <iosfwd>                           // for streamsize
#include <iterator>                         // for insert_iterator, inserter
#include <memory>                           // for unique_ptr, shared_ptr, __shared_ptr_access
#include <stdexcept>                        // for runtime_error, overflow_error
#include <string>                           // for string, basic_string
#include <string_view>                      // for string_view
#include <thread_pool/thread_pool.hpp>      // for thread_pool
#include <utility>                          // for tuple_element<>::type, pair, make_pair
#include <vector>                           // for vector

#include "modle/bed.hpp"                // for BED_tree, BED_tree<>::value_type, Parser
#include "modle/bigwig.hpp"             // for Writer
#include "modle/common/common.hpp"      // for u32, usize, bp_t, u8, i64
#include "modle/common/utils.hpp"       // for identity::operator()
#include "modle/contacts.hpp"           // for ContactMatrix
#include "modle/cooler.hpp"             // for Cooler, Cooler::READ_ONLY
#include "modle/interval_tree.hpp"      // for IITree, IITree::IITree<I, T>, IITree::empty
#include "modle/stats/correlation.hpp"  // for Pearson, Spearman
#include "modle_tools/config.hpp"       // for eval_config
#include "modle_tools/tools.hpp"        // for eval_subcmd

namespace modle::tools {

using ChromSet = absl::btree_map<std::string, std::pair<bp_t, bp_t>, cppsort::natural_less_t>;

[[nodiscard]] static ChromSet import_chrom_subranges(
    const boost::filesystem::path &path_to_chrom_subranges) {
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

[[nodiscard]] static ChromSet import_chroms_from_cool(cooler::Cooler &cooler_file) {
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

[[nodiscard]] static ChromSet select_chromosomes_for_eval(
    cooler::Cooler &c1, cooler::Cooler &c2,
    const boost::filesystem::path &path_to_chrom_subranges) {
  // Import chromosomes from Cooler files and select shared chromosomes
  const auto chrom_set1 = import_chroms_from_cool(c1);
  const auto chrom_set2 = import_chroms_from_cool(c2);

  ChromSet chrom_intersection;
  absl::c_set_intersection(chrom_set1, chrom_set2,
                           std::inserter(chrom_intersection, chrom_intersection.begin()));

  if (chrom_intersection.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Files {} and {} have 0 chromosomes in common. "
                   "Make sure both files are using the same genome assembly as reference"),
        c1.get_path(), c2.get_path()));
  }

  // Look for conflicting information in Cooler files
  std::string err;
  for (const auto &[name, range] : chrom_intersection) {
    if (chrom_set1.at(name).second != chrom_set2.at(name).second) {
      err += fmt::format(FMT_STRING("\n - \"{}\" has size {} in file #1 and {} in file #2"), name,
                         chrom_set1.at(name).second, chrom_set2.at(name).second);
    }
  }

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

  // Look for conflicting information in Cooler files and BED file
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

[[nodiscard]] static std::unique_ptr<io::bigwig::Writer> create_bwig_file(
    std::vector<std::pair<std::string, usize>> &chrom_list, std::string_view base_name,
    std::string_view suffix, bool skip) {
  if (skip) {
    return nullptr;
  }
  auto bw = std::make_unique<io::bigwig::Writer>(
      fmt::format(FMT_STRING("{}_{}"), base_name, absl::StripPrefix(suffix, "_")));
  bw->write_chromosomes(chrom_list);
  return bw;
}

enum CorrMethod : std::uint_fast8_t { pearson, spearman };
enum StripeDirection : std::uint_fast8_t { vertical, horizontal };

template <CorrMethod correlation_method, StripeDirection stripe_direction>
static std::vector<double> compute_correlation(const bp_t bin_size,
                                               std::shared_ptr<ContactMatrix<>> ref_contacts,
                                               std::shared_ptr<ContactMatrix<>> tgt_contacts,
                                               const bed::BED_tree<>::value_type &features,
                                               const std::vector<double> &weights = {}) {
  assert(ref_contacts->nrows() == tgt_contacts->nrows());              // NOLINT
  assert(weights.empty() || weights.size() == ref_contacts->nrows());  // NOLINT
  const auto nrows = ref_contacts->nrows();
  const auto ncols = ref_contacts->ncols();

  std::vector<u32> ref_pixel_buff(nrows);
  std::vector<u32> tgt_pixel_buff(nrows);
  std::vector<double> correlation_buff(ncols, 0.0);

  auto corr_fx = [&]() {
    if constexpr (correlation_method == pearson) {
      return stats::Pearson<double>{};
    } else {
      static_assert(correlation_method == spearman);
      return stats::Spearman<double>{nrows};
    }
  }();

  bp_t start_pos = 0;
  bp_t end_pos = bin_size;
  for (usize i = 0; i < ncols; ++i) {
    if (features.empty() || features.overlaps_with(start_pos, end_pos)) {
      if constexpr (stripe_direction == vertical) {
        ref_contacts->unsafe_get_column(i, ref_pixel_buff);
        tgt_contacts->unsafe_get_column(i, tgt_pixel_buff);
      } else {
        static_assert(stripe_direction == horizontal);
        ref_contacts->unsafe_get_row(i, ref_pixel_buff);
        tgt_contacts->unsafe_get_row(i, tgt_pixel_buff);
      }
      ref_pixel_buff.resize(nrows);
      tgt_pixel_buff.resize(nrows);

      [[maybe_unused]] const auto &[corr, pval] =
          weights.empty() ? corr_fx(ref_pixel_buff, tgt_pixel_buff)
                          : corr_fx(ref_pixel_buff, tgt_pixel_buff, weights);
      correlation_buff[i] = corr;
    }
    start_pos += bin_size;
    end_pos += bin_size;
  }
  return correlation_buff;
}

void eval_subcmd(const modle::tools::eval_config &c) {
  assert(c.compute_spearman || c.compute_pearson || c.compute_edist);  // NOLINT
  auto ref_cooler =
      cooler::Cooler(c.path_to_reference_matrix, cooler::Cooler::READ_ONLY, c.bin_size);
  auto tgt_cooler = cooler::Cooler(c.path_to_input_matrix, cooler::Cooler::READ_ONLY, c.bin_size);

  const auto bin_size = ref_cooler.get_bin_size();

  // This cannot be made const
  auto chromosomes = select_chromosomes_for_eval(ref_cooler, tgt_cooler, c.path_to_chrom_subranges);

  const auto features = [&]() {
    if (c.path_to_features_bed.empty()) {
      return bed::BED_tree<>{};
    }

    // TODO drop chromosomes not found in chromosomes
    return bed::Parser(c.path_to_features_bed).parse_all_in_interval_tree();
  }();

  if (chromosomes.size() == 1) {
    spdlog::info(FMT_STRING("Computing correlation for chromosome: \"{}\""),
                 chromosomes.begin()->first);
  } else {
    std::vector<std::string> chrom_names(static_cast<usize>(chromosomes.size()));
    std::transform(chromosomes.begin(), chromosomes.end(), chrom_names.begin(),
                   [](const auto &chrom) { return chrom.first; });
    spdlog::info(FMT_STRING("Computing correlation for the following {} chromosomes: \"{}\""),
                 chromosomes.size(), fmt::join(chrom_names, "\", \""));
  }

  if (const auto &output_dir = c.output_prefix.parent_path(); !output_dir.empty()) {
    boost::filesystem::create_directories(output_dir.string());
  }

  // Init files and write bw header
  std::vector<std::pair<std::string, usize>> chrom_vect(static_cast<usize>(chromosomes.size()));
  std::transform(chromosomes.begin(), chromosomes.end(), chrom_vect.begin(), [](const auto &chrom) {
    return std::make_pair(chrom.first, chrom.second.second);
  });
  auto init_bw = [&, bn = c.output_prefix.string()](const auto &suffix, bool skip) {
    struct bw {
      std::vector<double> buff{};
      std::unique_ptr<io::bigwig::Writer> writer{nullptr};
    };
    bw bw;
    bw.writer = create_bwig_file(chrom_vect, bn, suffix, skip);
    return bw;
  };

  auto bw_pearson_vertical = init_bw("pearson_vertical.bw", !c.compute_pearson);
  auto bw_pearson_horizontal = init_bw("pearson_horizontal.bw", !c.compute_pearson);
  auto bw_spearman_vertical = init_bw("spearman_vertical.bw", !c.compute_pearson);
  auto bw_spearman_horizontal = init_bw("spearman_horizontal.bw", !c.compute_pearson);
  auto bw_sed_vertical = init_bw("eucl_dist_vertical.bw", !c.compute_edist);
  auto bw_sed_horizontal = init_bw("eucl_dist_horizontal.bw", !c.compute_edist);

  const auto nrows = (c.diagonal_width + bin_size - 1) / bin_size;
  assert(nrows != 0);  // NOLINT

  const auto t0 = absl::Now();
  thread_pool tpool(6);  // NOLINT
  for (const auto &[chrom_name_sv, chrom_range] : chromosomes) {
    const auto t00 = absl::Now();
    const std::string chrom_name{chrom_name_sv};
    auto ref_contacts = std::make_shared<ContactMatrix<>>(
        ref_cooler.cooler_to_cmatrix(chrom_name, nrows, chrom_range));
    auto tgt_contacts = std::make_shared<ContactMatrix<>>(
        tgt_cooler.cooler_to_cmatrix(chrom_name, nrows, chrom_range));
    fmt::print(stderr, FMT_STRING("Read contacts for {} in {}\n"), chrom_name,
               absl::FormatDuration(absl::Now() - t00));

    const auto &chrom_features = [&]() {
      if (!features.contains(chrom_name)) {
        return bed::BED_tree<>::value_type{};
      }
      return features.at(chrom_name);
    }();

    // TODO handle exceptions in threads
    if (const auto &bw_fp = bw_pearson_vertical.writer; bw_fp) {
      tpool.push_task([&]() {
        auto t1 = absl::Now();
        auto corr = compute_correlation<pearson, vertical>(bin_size, ref_contacts, tgt_contacts,
                                                           chrom_features);
        spdlog::info(
            FMT_STRING("Pearson correlation for vertical stripes from chrom {} computed in {}."),
            chrom_name, absl::FormatDuration(absl::Now() - t1));
        t1 = absl::Now();
        bw_fp->write_range(chrom_name, absl::MakeSpan(corr), bin_size, bin_size);
        spdlog::info(FMT_STRING("Correlation values written to file {} in {}."), bw_fp->path(),
                     absl::FormatDuration(absl::Now() - t1));
      });
    }
    if (const auto &bw_fp = bw_pearson_horizontal.writer; bw_fp) {
      tpool.push_task([&]() {
        auto t1 = absl::Now();
        auto corr = compute_correlation<pearson, horizontal>(bin_size, ref_contacts, tgt_contacts,
                                                             chrom_features);
        spdlog::info(
            FMT_STRING("Pearson correlation for horizontal stripes from chrom {} computed in {}."),
            chrom_name, absl::FormatDuration(absl::Now() - t1));
        t1 = absl::Now();
        bw_fp->write_range(chrom_name, absl::MakeSpan(corr), bin_size, bin_size);
        spdlog::info(FMT_STRING("Correlation values written to file {} in {}."), bw_fp->path(),
                     absl::FormatDuration(absl::Now() - t1));
      });
    }
    if (const auto &bw_fp = bw_spearman_vertical.writer; bw_fp) {
      tpool.push_task([&]() {
        auto t1 = absl::Now();
        auto corr = compute_correlation<spearman, vertical>(bin_size, ref_contacts, tgt_contacts,
                                                            chrom_features);
        spdlog::info(
            FMT_STRING("Spearman correlation for vertical stripes from chrom {} computed in {}."),
            chrom_name, absl::FormatDuration(absl::Now() - t1));
        t1 = absl::Now();
        bw_fp->write_range(chrom_name, absl::MakeSpan(corr), bin_size, bin_size);
        spdlog::info(FMT_STRING("Correlation values written to file {} in {}."), bw_fp->path(),
                     absl::FormatDuration(absl::Now() - t1));
      });
    }
    if (const auto &bw_fp = bw_spearman_horizontal.writer; bw_fp) {
      tpool.push_task([&]() {
        auto t1 = absl::Now();
        auto corr = compute_correlation<spearman, horizontal>(bin_size, ref_contacts, tgt_contacts,
                                                              chrom_features);
        spdlog::info(
            FMT_STRING("Spearman correlation for horizontal stripes from chrom {} computed in {}."),
            chrom_name, absl::FormatDuration(absl::Now() - t1));
        t1 = absl::Now();
        bw_fp->write_range(chrom_name, absl::MakeSpan(corr), bin_size, bin_size);
        spdlog::info(FMT_STRING("Correlation values written to file {} in {}."), bw_fp->path(),
                     absl::FormatDuration(absl::Now() - t1));
      });
    }
    tpool.wait_for_tasks();
  }
  spdlog::info(FMT_STRING("DONE in {}!"), absl::FormatDuration(absl::Now() - t0));
}

}  // namespace modle::tools
