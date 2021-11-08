// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/algorithm/container.h>           // for c_set_intersection
#include <absl/container/btree_map.h>           // for btree_map, btree_iterator, map_params<>::...
#include <absl/container/flat_hash_map.h>       // for flat_hash_map
#include <absl/strings/ascii.h>                 // AsciiStrToLower
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
#include <future>                           // for future
#include <iosfwd>                           // for streamsize
#include <iterator>                         // for insert_iterator, inserter
#include <memory>                           // for unique_ptr, shared_ptr, __shared_ptr_access
#include <stdexcept>                        // for runtime_error, overflow_error
#include <string>                           // for string, basic_string
#include <string_view>                      // for string_view
#include <thread_pool/thread_pool.hpp>      // for thread_pool
#include <utility>                          // for tuple_element<>::type, pair, make_pair
#include <vector>                           // for vector

#include "modle/bed/bed.hpp"            // for BED_tree, BED_tree<>::value_type, Parser
#include "modle/bigwig/bigwig.hpp"      // for Writer
#include "modle/common/common.hpp"      // for u32, usize, bp_t, u8, i64
#include "modle/common/utils.hpp"       // for identity::operator()
#include "modle/contacts.hpp"           // for ContactMatrix
#include "modle/cooler/cooler.hpp"      // for Cooler, Cooler::READ_ONLY
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

[[nodiscard]] static ChromSet import_chroms_from_cool(cooler::Cooler<> &cooler_file) {
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
    cooler::Cooler<> &c1, cooler::Cooler<> &c2,
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

[[nodiscard]] static isize find_col_idx(const boost::filesystem::path &path_to_weights,
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
[[nodiscard]] static isize find_col_idx(const boost::filesystem::path &path_to_weights,
                                        std::string_view header, const Range &col_names) {
  assert(col_names.size() > 1);  // NOLINT
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
    const boost::filesystem::path &path_to_weights, const std::string_view weight_column_name,
    const usize nbins, const bool reciprocal_weights) {
  assert(nbins != 0);  // NOLINT
  absl::flat_hash_map<std::string, std::vector<double>> weights;

  if (path_to_weights.empty()) {
    return weights;
  }

  assert(boost::filesystem::exists(path_to_weights));  // NOLINT
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
    std::vector<std::pair<std::string, usize>> &chrom_list, std::string_view base_name,
    std::string_view suffix) {
  auto bw = io::bigwig::Writer(
      fmt::format(FMT_STRING("{}_{}"), base_name, absl::StripPrefix(suffix, "_")));
  bw.write_chromosomes(chrom_list);
  return bw;
}

enum CorrMethod : std::uint_fast8_t { pearson, spearman, eucl_dist };
enum StripeDirection : std::uint_fast8_t { vertical, horizontal };

template <class N, class = std::enable_if_t<std::is_arithmetic_v<N>>>
static size_t mask_zero_pixels(const std::vector<N> &v1, const std::vector<N> &v2,
                               std::vector<double> &weights) {
  assert(v1.size() == v2.size());       // NOLINT
  assert(v1.size() == weights.size());  // NOLINT

  usize masked_pixels = 0;
  for (usize i = 0; i < v1.size(); ++i) {
    if (v1[i] == N(0) || v2[i] == N(0)) {
      ++masked_pixels;
      weights[i] = 0;
    }
  }
  return masked_pixels;
}

template <CorrMethod correlation_method, StripeDirection stripe_direction, class N,
          class WeightIt = utils::RepeatIterator<double>>
static std::vector<double> compute_correlation(
    std::shared_ptr<ContactMatrix<N>> ref_contacts, std::shared_ptr<ContactMatrix<N>> tgt_contacts,
    const bool mask_zero_pixels_, WeightIt weight_first = utils::RepeatIterator<double>(1)) {
  assert(ref_contacts->nrows() == tgt_contacts->nrows());  // NOLINT
  const auto nrows = ref_contacts->nrows();
  const auto ncols = ref_contacts->ncols();

  std::vector<double> ref_pixel_buff(nrows);
  std::vector<double> tgt_pixel_buff(nrows);
  std::vector<double> weight_buff;
  if (!std::is_same_v<WeightIt, utils::RepeatIterator<double>> || mask_zero_pixels_) {
    weight_buff.resize(nrows);
    std::copy_n(weight_first, nrows, weight_buff.begin());
  }

  std::vector<double> correlation_buff(ncols, 0.0);

  auto corr_fx = [&]() {
    if constexpr (correlation_method == pearson) {
      return stats::Pearson<double>{};
    } else {
      static_assert(correlation_method == spearman);
      return stats::Spearman<double>{nrows};
    }
  }();

  for (usize i = 0; i < ncols; ++i) {
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

    if (mask_zero_pixels_) {
      std::copy_n(weight_first, nrows, weight_buff.begin());
      mask_zero_pixels(ref_pixel_buff, tgt_pixel_buff, weight_buff);
    }

    [[maybe_unused]] const auto &[corr, pval] =
        !weight_buff.empty() ? corr_fx(ref_pixel_buff.begin(), ref_pixel_buff.end(),
                                       tgt_pixel_buff.begin(), weight_buff.begin())
                             : corr_fx(ref_pixel_buff.begin(), ref_pixel_buff.end(),
                                       tgt_pixel_buff.begin(), weight_first);
    correlation_buff[i] = corr;
  }
  return correlation_buff;
}

template <CorrMethod correlation_method, StripeDirection stripe_direction, class N>
static std::vector<double> compute_correlation(std::shared_ptr<ContactMatrix<N>> ref_contacts,
                                               std::shared_ptr<ContactMatrix<N>> tgt_contacts,
                                               const bool mask_zero_pixels_,
                                               const std::vector<double> &weights) {
  assert(ref_contacts->nrows() == tgt_contacts->nrows());              // NOLINT
  assert(weights.empty() || ref_contacts->nrows() == weights.size());  // NOLINT
  return compute_correlation<correlation_method, stripe_direction>(
      ref_contacts, tgt_contacts, mask_zero_pixels_, weights.begin());
}

absl::flat_hash_map<std::pair<CorrMethod, StripeDirection>, io::bigwig::Writer> init_bwigs(
    const modle::tools::eval_config &c, const ChromSet &chroms, const bool weighted) {
  std::vector<std::pair<std::string, usize>> chrom_vect(static_cast<usize>(chroms.size()));
  std::transform(chroms.begin(), chroms.end(), chrom_vect.begin(), [](const auto &chrom) {
    return std::make_pair(chrom.first, chrom.second.second);
  });

  absl::flat_hash_map<std::pair<CorrMethod, StripeDirection>, io::bigwig::Writer> bwigs;
  if (c.compute_pearson) {
    bwigs.emplace(std::make_pair(pearson, vertical),
                  create_bwig_file(chrom_vect, c.output_prefix.string(),
                                   fmt::format(FMT_STRING("pearson{}_vertical.bw"),
                                               weighted ? "_weighted" : "")));
    bwigs.emplace(std::make_pair(pearson, horizontal),
                  create_bwig_file(chrom_vect, c.output_prefix.string(),
                                   fmt::format(FMT_STRING("pearson{}_horizontal.bw"),
                                               weighted ? "_weighted" : "")));
  }

  if (c.compute_spearman) {
    bwigs.emplace(std::make_pair(spearman, vertical),
                  create_bwig_file(chrom_vect, c.output_prefix.string(),
                                   fmt::format(FMT_STRING("spearman{}_vertical.bw"),
                                               weighted ? "_weighted" : "")));
    bwigs.emplace(std::make_pair(spearman, horizontal),
                  create_bwig_file(chrom_vect, c.output_prefix.string(),
                                   fmt::format(FMT_STRING("spearman{}_horizontal.bw"),
                                               weighted ? "_weighted" : "")));
  }
  /*
  bwigs.emplace(std::make_pair(eucl_dist, vertical),
                create_bwig_file(chrom_vect, c.output_prefix.string(), "eucl_dist_vertical.bw",
                                 !c.compute_eucl_dist));
  bwigs.emplace(std::make_pair(eucl_dist, horizontal),
                create_bwig_file(chrom_vect, c.output_prefix.string(), "eucl_dist_horizontal.bw",
                                 !c.compute_eucl_dist));
 */
  return bwigs;
}

[[nodiscard]] [[maybe_unused]] static constexpr std::string_view corr_method_to_str(
    const CorrMethod m) {
  if (m == pearson) {
    return "Pearson correlation";
  }
  if (m == spearman) {
    return "Spearman correlation";
  }

  if (m == eucl_dist) {
    return "Euclidean distance";
  }
  std::abort();
}

[[nodiscard]] [[maybe_unused]] static constexpr std::string_view direction_to_str(
    const StripeDirection m) {
  if (m == vertical) {
    return "Vertical";
  }
  if (m == horizontal) {
    return "Horizontal";
  }
  std::abort();
}

template <CorrMethod correlation_method, StripeDirection stripe_direction, class N>
[[nodiscard]] static std::future<bool> submit_task(
    const std::string &chrom_name,
    absl::flat_hash_map<std::pair<CorrMethod, StripeDirection>, io::bigwig::Writer> &bwigs,
    thread_pool &tpool, std::shared_ptr<ContactMatrix<N>> &ref_contacts,
    std::shared_ptr<ContactMatrix<N>> &tgt_contacts, const bool exclude_zero_pixels,
    const usize bin_size, const absl::flat_hash_map<std::string, std::vector<double>> &weights) {
  return tpool.submit([&]() {
    try {
      auto t0 = absl::Now();
      auto corr = [&]() {
        if (weights.empty()) {
          return compute_correlation<correlation_method, stripe_direction>(
              ref_contacts, tgt_contacts, exclude_zero_pixels);
        }
        return compute_correlation<correlation_method, stripe_direction>(
            ref_contacts, tgt_contacts, exclude_zero_pixels, weights.at(chrom_name));
      }();

      spdlog::info(FMT_STRING("{} for {} stripes from chrom {} computed in {}."),
                   corr_method_to_str(correlation_method),
                   absl::AsciiStrToLower(direction_to_str(stripe_direction)), chrom_name,
                   absl::FormatDuration(absl::Now() - t0));
      t0 = absl::Now();
      auto &bw = bwigs.at({correlation_method, stripe_direction});
      bw.write_range(chrom_name, absl::MakeSpan(corr), bin_size, bin_size);
      spdlog::info(FMT_STRING("{} values have been written to file {} in {}."), corr.size(),
                   bw.path(), absl::FormatDuration(absl::Now() - t0));
    } catch (const std::exception &e) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("The following error occurred while computing {} on {} stripes for {}: {}"),
          corr_method_to_str(correlation_method),
          absl::AsciiStrToLower(direction_to_str(stripe_direction)), chrom_name, e.what()));
    }
  });
}

void eval_subcmd(const modle::tools::eval_config &c) {
  assert(c.compute_spearman || c.compute_pearson /*|| c.compute_eucl_dist*/);  // NOLINT
  auto ref_cooler =
      cooler::Cooler(c.path_to_reference_matrix, cooler::Cooler<>::IO_MODE::READ_ONLY, c.bin_size);
  auto tgt_cooler =
      cooler::Cooler(c.path_to_input_matrix, cooler::Cooler<>::IO_MODE::READ_ONLY, c.bin_size);

  const auto bin_size = ref_cooler.get_bin_size();

  // This cannot be made const
  auto chromosomes = select_chromosomes_for_eval(ref_cooler, tgt_cooler, c.path_to_chrom_subranges);

  const auto weights =
      import_weights(c.path_to_weights, c.weight_column_name,
                     (c.diagonal_width + bin_size - 1) / bin_size, c.reciprocal_weights);
  if (!weights.empty()) {
    std::vector<std::string> missing_chroms;
    for (const auto &[chrom, _] : chromosomes) {
      if (!weights.contains(chrom)) {
        missing_chroms.push_back(chrom);
      }
    }
    if (!missing_chroms.empty()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Unable to import weights from file {} for the following chromosomes: {}"),
          c.path_to_weights, fmt::join(missing_chroms, ", ")));
    }
  }

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

  auto bwigs = init_bwigs(c, chromosomes, !weights.empty());

  const auto nrows = (c.diagonal_width + bin_size - 1) / bin_size;
  assert(nrows != 0);  // NOLINT

  const auto t0 = absl::Now();
  thread_pool tpool(c.nthreads);
  std::vector<std::future<bool>> return_codes(6);
  for (const auto &[chrom_name_sv, chrom_range] : chromosomes) {
    const std::string chrom_name{chrom_name_sv};

    const auto t1 = absl::Now();
    spdlog::info(FMT_STRING("Reading contacts for {}..."), chrom_name_sv);

    auto ref_contacts = std::make_shared<ContactMatrix<double>>(
        ref_cooler.cooler_to_cmatrix(chrom_name, nrows, chrom_range)
            .unsafe_gaussian_diff(1.0, 1.6));
    auto tgt_contacts = std::make_shared<ContactMatrix<double>>(
        tgt_cooler.cooler_to_cmatrix(chrom_name, nrows, chrom_range)
            .unsafe_gaussian_diff(1.0, 1.6));
    spdlog::info(FMT_STRING("Read {} contacts for {} in {}"),
                 ref_contacts->get_tot_contacts() + tgt_contacts->get_tot_contacts(), chrom_name_sv,
                 absl::FormatDuration(absl::Now() - t1));

    return_codes.clear();

    if (bwigs.contains({spearman, horizontal})) {
      return_codes.emplace_back(
          submit_task<spearman, horizontal>(chrom_name, bwigs, tpool, ref_contacts, tgt_contacts,
                                            c.exclude_zero_pxls, bin_size, weights));
    }

    if (bwigs.contains({spearman, vertical})) {
      return_codes.emplace_back(
          submit_task<spearman, vertical>(chrom_name, bwigs, tpool, ref_contacts, tgt_contacts,
                                          c.exclude_zero_pxls, bin_size, weights));
    }

    if (bwigs.contains({pearson, horizontal})) {
      return_codes.emplace_back(
          submit_task<pearson, horizontal>(chrom_name, bwigs, tpool, ref_contacts, tgt_contacts,
                                           c.exclude_zero_pxls, bin_size, weights));
    }

    if (bwigs.contains({pearson, vertical})) {
      return_codes.emplace_back(
          submit_task<pearson, vertical>(chrom_name, bwigs, tpool, ref_contacts, tgt_contacts,
                                         c.exclude_zero_pxls, bin_size, weights));
    }

    tpool.wait_for_tasks();
    // Catch exception raised in worker threads
    for (auto &code : return_codes) {
      code.get();
    }
  }
  spdlog::info(FMT_STRING("DONE in {}!"), absl::FormatDuration(absl::Now() - t0));
}

}  // namespace modle::tools
