#include <absl/container/flat_hash_map.h>  // for flat_hash_map, BitMask
#include <absl/container/flat_hash_set.h>  // for flat_hash_set
#include <absl/strings/str_cat.h>          // for StrCat
#include <absl/strings/str_join.h>         // for StrJoin
#include <absl/strings/strip.h>            // for StripSuffix
#include <absl/types/span.h>               // for Span
#include <fmt/format.h>                    // for print, FMT_STRING
#include <fmt/ostream.h>                   // for print

#include <algorithm>                                // for transform
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/filesystem/path.hpp>                // for path, create_directories, remove_all
#include <cassert>                                  // for assert
#include <cmath>                                    // for round
#include <cstddef>                                  // for size_t
#include <cstdint>                                  // for uint64_t
#include <cstdio>                                   // for stdout
#include <fstream>                                  // for ofstream, basic_ofstream
#include <iterator>                                 // for insert_iterator, inserter
#include <memory>                                   // for unique_ptr, make_unique
#include <numeric>                                  // for iota
#include <string>                                   // for string, basic_string
#include <utility>                                  // for make_pair
#include <vector>                                   // for vector

#include "modle/bed.hpp"  // for Parser
#include "modle/common/suppress_compiler_warnings.hpp"
#include "modle/contacts.hpp"      // for ContactMatrix
#include "modle/cooler.hpp"        // for Cooler, ContactMatrix, Cooler::READ_ONLY
#include "modle_tools/config.hpp"  // for config
#include "modle_tools/tools.hpp"   // for stats_subcmd

namespace modle::tools {

template <typename I>
size_t compute_number_of_contacts_after_depletion(const ContactMatrix<I>& cmatrix,
                                                  absl::Span<const uint64_t> hist,
                                                  size_t effective_nbins,
                                                  double depletion_multiplier) {
  assert(hist.size() == cmatrix.nrows());
  // This histogram contains the average contact number (instead of the total)
  std::vector<uint64_t> row_wise_avg_contacts(hist.size());
  std::transform(hist.begin(), hist.end(), row_wise_avg_contacts.begin(), [&](const auto n) {
    return static_cast<uint64_t>(std::round((depletion_multiplier * static_cast<double>(n)) /
                                            static_cast<double>(effective_nbins)));
  });

  size_t depl_contacts{0};
  for (auto i = 0UL; i < cmatrix.ncols(); ++i) {
    for (auto j = i; j < i + cmatrix.nrows() && j < cmatrix.ncols(); ++j) {
      // j - i corresponds to the distance from the diagonal
      // Only process bins that more contacts than average
      if (const auto n = cmatrix.get(j, i); n > row_wise_avg_contacts[j - i]) {
        depl_contacts += n - row_wise_avg_contacts[j - i];
      }
    }
  }

  return depl_contacts;
}

void stats_subcmd(const modle::tools::config& c) {
  const auto& path_to_output_hist = c.output_path_for_histograms;

  cooler::Cooler m1(c.path_to_input_matrix, cooler::Cooler::READ_ONLY, c.bin_size);

  std::unique_ptr<cooler::Cooler> m2{nullptr};
  std::unique_ptr<std::ofstream> hist_file{nullptr};

  const auto nchroms = m1.get_nchroms();
  const auto bin_size = m1.get_bin_size();
  const auto chrom_names = m1.get_chrom_names();
  const auto chrom_sizes = m1.get_chrom_sizes();

  absl::flat_hash_map<std::string, size_t> chrom_start_offsets;
  if (!c.path_to_chrom_subranges.empty()) {
    modle::bed::Parser p(c.path_to_chrom_subranges);
    const auto buff = p.parse_all();
    std::transform(
        buff.begin(), buff.end(), std::inserter(chrom_start_offsets, chrom_start_offsets.end()),
        [](const auto& record) { return std::make_pair(record.chrom, record.chrom_start); });
  }

  if (c.dump_depleted_matrices) {  // Create cooler file to write depl. contacts
    const auto ext = c.path_to_input_matrix.extension().string();
    const auto path =
        absl::StrCat(absl::StripSuffix(c.path_to_input_matrix.string(), ext), "_depl.cool");
    boost::filesystem::remove_all(path);
    m2 = std::make_unique<cooler::Cooler>(path, cooler::Cooler::WRITE_ONLY, bin_size);
  }

  if (!path_to_output_hist.empty()) {  // Create hist. file
    boost::filesystem::create_directories(
        boost::filesystem::path(path_to_output_hist).parent_path());
    hist_file = std::make_unique<std::ofstream>(path_to_output_hist);
    std::vector<size_t> buff(c.diagonal_width / bin_size);
    std::iota(buff.begin(), buff.end(), 0);
    // TODO: Consider whether it make sense to write this header
    fmt::print(*hist_file, "#{}\n", absl::StrJoin(buff, "\t"));
  }

  // Write header
  fmt::print(
      stdout,
      "chrom_name\ttot_number_of_contacts_1\ttot_number_of_contacts_2\tavg_number_of_contacts_"
      "1\tavg_number_of_contacts_2\tfraction_of_graylisted_bins\n");

  size_t tot_contacts = 0;
  size_t tot_contacts_after_depl = 0;
  size_t tot_number_of_pixels = 0;

  for (auto i = 0UL; i < nchroms; ++i) {
    const auto& chrom_name = chrom_names[i];
    const auto& chrom_size = chrom_sizes[i];

    // Skip chromosome found in the exclusion list
    if (c.chromosomes_excluded.contains(chrom_name)) {
      continue;
    }

    // Read contacts for chrom_name into memory
    auto cmatrix = m1.cooler_to_cmatrix(chrom_name, c.diagonal_width, bin_size);

    const auto hist = cmatrix.compute_row_wise_contact_histogram();
    const auto mask = cmatrix.generate_mask_for_bins_without_contacts();

    const auto chrom_contacts = cmatrix.get_tot_contacts();
    const auto chrom_contacts_after_depl = compute_number_of_contacts_after_depletion(
        cmatrix, hist, mask.count(), c.depletion_multiplier);
    assert(chrom_contacts_after_depl <= chrom_contacts);  // NOLINT

    const auto chrom_avg_contacts =
        static_cast<double>(chrom_contacts) / static_cast<double>(cmatrix.npixels_after_masking());
    const auto chrom_avg_contacts_after_depl = static_cast<double>(chrom_contacts_after_depl) /
                                               static_cast<double>(cmatrix.npixels_after_masking());
    const auto fraction_of_graylisted_bins =
        1.0 - (static_cast<double>(cmatrix.npixels_after_masking()) /  // NOLINT
               static_cast<double>(cmatrix.npixels()));

    tot_contacts += chrom_contacts;
    tot_contacts_after_depl += chrom_contacts_after_depl;
    tot_number_of_pixels += cmatrix.npixels_after_masking();

    // clang-format off
    fmt::print(stdout, FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\n"),
               chrom_name,
               chrom_contacts,
               chrom_contacts_after_depl,
               chrom_avg_contacts,
               chrom_avg_contacts_after_depl,
               fraction_of_graylisted_bins);
    // clang-format on

    if (hist_file) {
      fmt::print(*hist_file, "{}\n", absl::StrJoin(hist, "\t"));
    }

    if (m2) {
      cmatrix.deplete_contacts(c.depletion_multiplier);
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_USELESS_CAST
      m2->write_or_append_cmatrix_to_file(cmatrix, chrom_name, int64_t(0), chrom_size, chrom_size,
                                          true);
      DISABLE_WARNING_POP
    }
  }

  fmt::print(
      stdout, FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\n"), "grand_average", tot_contacts,
      tot_contacts_after_depl,
      static_cast<double>(tot_contacts) / static_cast<double>(tot_number_of_pixels),
      static_cast<double>(tot_contacts_after_depl) / static_cast<double>(tot_number_of_pixels), 0);
}

}  // namespace modle::tools
