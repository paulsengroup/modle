// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/algorithm/container.h>           // for c_set_intersection
#include <absl/container/btree_map.h>           // for btree_map, btree_iterator, map_params<>::...
#include <absl/container/flat_hash_map.h>       // for flat_hash_map
#include <absl/strings/ascii.h>                 // AsciiStrToLower
#include <absl/strings/str_split.h>             // for StrSplit
#include <absl/strings/strip.h>                 // for StripPrefix
#include <absl/time/clock.h>                    // for Now
#include <absl/time/time.h>                     // for FormatDuration, operator-, Time
#include <absl/types/span.h>                    // for MakeSpan
#include <absl/types/variant.h>                 // for get, variant
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

#include "modle/bed/bed.hpp"                      // for BED_tree, BED_tree<>::value_type, Parser
#include "modle/bigwig/bigwig.hpp"                // for Writer
#include "modle/common/common.hpp"                // for u32, usize, bp_t, u8, i64
#include "modle/common/utils.hpp"                 // for identity::operator()
#include "modle/compressed_io/compressed_io.hpp"  // for Reader
#include "modle/contacts.hpp"                     // for ContactMatrix
#include "modle/cooler/cooler.hpp"                // for Cooler, Cooler::READ_ONLY
#include "modle/interval_tree.hpp"                // for IITree, IITree::IITree<I, T>, IITree::empty
#include "modle/stats/correlation.hpp"            // for Pearson, Spearman
#include "modle_tools/config.hpp"                 // for eval_config
#include "modle_tools/tools.hpp"                  // for eval_subcmd

namespace modle::tools {

static modle::IITree<double, double> import_discretization_ranges(const boost::filesystem::path& p,
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
    toks = absl::StrSplit(buff, sep);
    try {
      if (toks.size() != 3) {
        throw std::runtime_error(fmt::format(FMT_STRING("expected 3 tokens, got {}. "
                                                        "Invalid record: \"{}\""),
                                             i, p, toks.size(), buff));
      }

      const auto start = utils::parse_numeric_or_throw<double>(toks[0]);
      const auto end = utils::parse_numeric_or_throw<double>(toks[1]);
      const auto value = utils::parse_numeric_or_throw<double>(toks[2]);

      if (start >= end) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("the begin of a range should be less or equal than its end, "
                                   "found begin={}; end={}. Invalid record: \"{}\""),
                        start, end, buff));
      }
      ranges.insert(start, end, value);
    } catch (const std::exception& e) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Found an invalid record at line {} of file {}: {}"), i, p, e.what()));
    }
  }
  ranges.make_BST();
  spdlog::info(FMT_STRING("Imported {} ranges from file {}"), i - 1, p);
  return ranges;
}

template <class N>
[[nodiscard]] static ContactMatrix<N> run_normalization(
    std::string_view chrom_name, ContactMatrix<N>& m,
    const std::pair<double, double> normalization_range,
    const std::pair<double, double> saturation_range) {
  const auto [lower_bound_norm, upper_bound_norm] = normalization_range;
  const auto [lower_bound_sat, upper_bound_sat] = saturation_range;

  spdlog::info(FMT_STRING("Normalizing contacts for {} to the {:.4g}-{:.4g} range..."), chrom_name,
               lower_bound_norm, upper_bound_norm);
  // Clamp contacts before normalizing
  if (!std::isinf(lower_bound_sat) || !std::isinf(upper_bound_sat)) {
    m.clamp_inplace(static_cast<N>(lower_bound_sat), static_cast<N>(upper_bound_sat));
  }
  if constexpr (std::is_same_v<N, double>) {
    m.normalize_inplace(lower_bound_norm, upper_bound_norm);
    return m;
  } else {
    return m.normalize(lower_bound_norm, upper_bound_norm);
  }
}

template <class N>
[[nodiscard]] static ContactMatrix<double> run_gaussian_blur(
    thread_pool& tpool, std::string_view chrom_name, ContactMatrix<N>& m, const double sigma,
    const std::pair<double, double> saturation_range) {
  const auto [lower_bound_sat, upper_bound_sat] = saturation_range;
  spdlog::info(FMT_STRING("Applying Gaussian blur with sigma={:.4g} to contacts for {}..."), sigma,
               chrom_name);
  auto m2 = m.blur(sigma, 0.005, &tpool);
  if (!std::isinf(lower_bound_sat) || !std::isinf(upper_bound_sat)) {
    m2.clamp_inplace(lower_bound_sat, upper_bound_sat);
  }
  return m2;
}

template <class N>
[[nodiscard]] static ContactMatrix<double> run_difference_of_gaussians(
    thread_pool& tpool, std::string_view chrom_name, ContactMatrix<N>& m, const double sigma1,
    const double sigma2, const std::pair<double, double> saturation_range) {
  const auto [lower_bound_sat, upper_bound_sat] = saturation_range;
  spdlog::info(FMT_STRING("Computing the difference of Gaussians for {} (sigma1={:.4g}; "
                          "sigma2={:.4g})..."),
               chrom_name, sigma1, sigma2);
  return m.gaussian_diff(sigma1, sigma2, lower_bound_sat, upper_bound_sat, &tpool);
}

using ContactMatrixVariant = absl::variant<ContactMatrix<double>, ContactMatrix<i32>>;
template <class N>
[[nodiscard]] static ContactMatrixVariant process_chromosome(
    thread_pool& tpool, const std::string_view chrom_name, const bp_t bin_size,
    cooler::Cooler<N>& cooler, std::mutex& cooler_mtx, const modle::tools::transform_config& c,
    const modle::IITree<double, double>& discretization_ranges) {
  auto t1 = absl::Now();
  spdlog::info(FMT_STRING("Processing contacts for {}..."), chrom_name);
  auto matrix = [&, chrom_name = chrom_name]() {
    std::unique_lock<std::mutex> lck(cooler_mtx);
    auto m = cooler.cooler_to_cmatrix(chrom_name, c.diagonal_width, bin_size);
    lck.unlock();
    using t = transform_config::transformation;
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
        throw std::logic_error("Unreachable code");
    }
  }();

  if (!discretization_ranges.empty()) {
    matrix.discretize_inplace(discretization_ranges);
  }

  spdlog::info(FMT_STRING("{} processing took {}"), chrom_name,
               absl::FormatDuration(absl::Now() - t1));
  if (c.floating_point) {
    return matrix;
  }
  return matrix.template as<i32>();
}

static void write_contacts_to_file(const modle::tools::transform_config& c, std::mutex& cooler_mtx,
                                   const bp_t bin_size,
                                   std::deque<std::future<ContactMatrixVariant>>& matrix_queue,
                                   std::mutex& matrix_queue_mtx,
                                   const std::vector<std::pair<std::string, bp_t>>& chroms) {
  assert(bin_size != 0);  // NOLINT
  usize i = 0;
  try {
    const auto max_chrom_name_size = [&]() {
      const auto it = std::max_element(
          chroms.begin(), chroms.end(),
          [](const auto& c1, const auto& c2) { return c1.first.size() < c2.first.size(); });
      return it->first.size();
    }();

    assert(c.force || !boost::filesystem::exists(c.path_to_output_matrix));  // NOLINT
    boost::filesystem::create_directories(c.path_to_output_matrix.parent_path());
    auto cooler_variant = [&]() {
      std::scoped_lock<std::mutex> lck(cooler_mtx);
      using CoolerV = absl::variant<cooler::Cooler<double>, cooler::Cooler<i32>>;
      if (c.floating_point) {
        return CoolerV(absl::in_place_type<cooler::Cooler<double>>, c.path_to_output_matrix,
                       cooler::Cooler<double>::IO_MODE::WRITE_ONLY, bin_size, max_chrom_name_size);
      }
      return CoolerV(absl::in_place_type<cooler::Cooler<i32>>, c.path_to_output_matrix,
                     cooler::Cooler<i32>::IO_MODE::WRITE_ONLY, bin_size, max_chrom_name_size);
    }();

    while (true) {
      const auto queue_is_empty = [&]() {
        std::scoped_lock<std::mutex> lck(matrix_queue_mtx);
        return matrix_queue.empty();
      }();

      if (queue_is_empty) {
        // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        continue;
      }

      const auto matrix_variant = [&]() {
        // Wait for the next matrix to become available. Also catch exceptions from other worker
        // threads
        auto m = matrix_queue.front().get();
        std::scoped_lock<std::mutex> lck(matrix_queue_mtx);
        matrix_queue.pop_front();
        return m;
      }();

      const auto& [chrom_name, chrom_size] = chroms[i];
      if (c.floating_point) {
        using M = double;
        const auto& matrix = absl::get<ContactMatrix<M>>(matrix_variant);
        if (matrix.empty()) {  // EOQ signal
          return;
        }
        auto& cooler = absl::get<cooler::Cooler<M>>(cooler_variant);
        std::scoped_lock<std::mutex> lck(cooler_mtx);
        cooler.write_or_append_cmatrix_to_file(matrix, chrom_name, usize(0), chrom_size,
                                               chrom_size);
      } else {
        using M = i32;
        const auto& matrix = absl::get<ContactMatrix<M>>(matrix_variant);
        if (matrix.empty()) {  // EOQ signal
          return;
        }
        auto& cooler = absl::get<cooler::Cooler<M>>(cooler_variant);
        std::scoped_lock<std::mutex> lck(cooler_mtx);
        cooler.write_or_append_cmatrix_to_file(matrix, chrom_name, usize(0), chrom_size,
                                               chrom_size);
      }
      ++i;
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("The following error occurred while writing contacts for {} to file {}: {}"),
        chroms[i].first, c.path_to_output_matrix, e.what()));
  }
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

  // The underlying HDF5 API is not yet thread safe, so we need to use this global mutex
  std::mutex input_cooler_mtx;  // Protect access to input and output coolers
  // std::mutex output_cooler_mtx;  // Protect access to input and output coolers
  auto input_cooler = cooler::Cooler<double>(
      c.path_to_input_matrix, cooler::Cooler<double>::IO_MODE::READ_ONLY, c.bin_size);
  const auto chroms = input_cooler.get_chroms();
  const auto bin_size = input_cooler.get_bin_size();

  thread_pool tpool(std::max(c.nthreads, usize(2)));

  std::mutex output_matrices_mtx;
  std::deque<std::future<ContactMatrixVariant>> output_matrices;

  // TODO rewrite this part using a conditional variable.
  auto writer_return_code = tpool.submit([&]() {
    write_contacts_to_file(c, input_cooler_mtx, bin_size, output_matrices, output_matrices_mtx,
                           chroms);
  });

  const auto t0 = absl::Now();
  spdlog::info(FMT_STRING("Transforming contacts from file {}..."), input_cooler.get_path());
  for (usize i = 0; i < chroms.size(); ++i) {
    /*
    auto result_fut = tpool.submit([&, chrom_name = chroms[i].first]() -> ContactMatrixVariant {
      try {
        return process_chromosome(tpool, chrom_name, bin_size, input_cooler, input_cooler_mtx, c,
                                  discretization_ranges);
      } catch (const std::exception& e) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("The following error occurred while processing {}: {}"),
                        chrom_name, e.what()));
      }
    });
     */
    auto result_fut =
        utils::make_ready_future(process_chromosome(tpool, chroms[i].first, bin_size, input_cooler,
                                                    input_cooler_mtx, c, discretization_ranges));
    std::scoped_lock<std::mutex> lck(output_matrices_mtx);
    output_matrices.emplace_back(std::move(result_fut));
  }

  {  // Signal EOQ
    ContactMatrixVariant eoq;
    if (c.floating_point) {
      eoq = ContactMatrix<double>{};
    } else {
      eoq = ContactMatrix<i32>{};
    }
    auto fut = utils::make_ready_future(std::move(eoq));
    std::scoped_lock<std::mutex> lck(output_matrices_mtx);
    output_matrices.emplace_back(std::move(fut));
  }

  tpool.wait_for_tasks();
  spdlog::info(FMT_STRING("DONE! Processed {} chromosomes in {}!"), chroms.size(),
               absl::FormatDuration(absl::Now() - t0));
  spdlog::info(FMT_STRING("Transformed contacts have been saved to file {}"),
               c.path_to_output_matrix);
}

}  // namespace modle::tools
