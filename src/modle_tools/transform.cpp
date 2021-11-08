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

void transform_subcmd(const modle::tools::transform_config& c) {
  auto input_cooler =
      cooler::Cooler(c.path_to_input_matrix, cooler::Cooler<>::IO_MODE::READ_ONLY, c.bin_size);
  const auto max_chrom_name_size = [&]() {
    const auto& chrom_names = input_cooler.get_chrom_names();
    const auto it =
        std::max_element(chrom_names.begin(), chrom_names.end(),
                         [](const auto& c1, const auto& c2) { return c1.size() < c2.size(); });
    return it->size();
  }();

  assert(c.force || !boost::filesystem::exists(c.path_to_output_matrix));  // NOLINT
  auto output_cooler =
      cooler::Cooler<double>(c.path_to_output_matrix, cooler::Cooler<double>::IO_MODE::WRITE_ONLY,
                             c.bin_size, max_chrom_name_size);

  const auto bin_size = input_cooler.get_bin_size();

  const auto t0 = absl::Now();
  spdlog::info(FMT_STRING("Transforming contacts from file {}..."), input_cooler.get_path());
  for (const auto& [chrom_name, chrom_size] : input_cooler.get_chroms()) {
    const auto m = [&, chrom_name = chrom_name]() {
      using t = transform_config::transformation;
      switch (c.method) {
        case t::normalize: {
          const auto lb = c.normalization_range.first;
          const auto ub = c.normalization_range.second;
          spdlog::info(FMT_STRING("Normalizing contacts for {} to the {:.4g}-{:.4g} range..."),
                       chrom_name, lb, ub);
          return input_cooler.cooler_to_cmatrix(chrom_name, c.diagonal_width, bin_size)
              .normalize(lb, ub);
        }
        case t::gaussian_blur:
          spdlog::info(FMT_STRING("Applying Gaussian blur with sigma={:.4g} to contacts for {}..."),
                       c.gaussian_blur_sigma, chrom_name);
          return input_cooler.cooler_to_cmatrix(chrom_name, c.diagonal_width, bin_size)
              .blur(c.gaussian_blur_sigma);
        case t::difference_of_gaussians: {
          const auto sigma1 = c.gaussian_blur_sigma;
          const auto sigma2 = sigma1 * c.gaussian_blur_sigma_multiplier;
          spdlog::info(
              FMT_STRING(
                  "Computing the difference of Gaussians for {} (sigma1={:.4g}; sigma2={:.4g})..."),
              chrom_name, sigma1, sigma2);
          return input_cooler.cooler_to_cmatrix(chrom_name, c.diagonal_width, bin_size)
              .unsafe_gaussian_diff(sigma1, sigma2);
        }
        default:
          throw std::logic_error("Unreachable code");
      }
    }();

    output_cooler.write_or_append_cmatrix_to_file(m, chrom_name, usize(0), chrom_size, chrom_size);
  }
  spdlog::info(FMT_STRING("DONE in {}!"), absl::FormatDuration(absl::Now() - t0));
}

}  // namespace modle::tools
