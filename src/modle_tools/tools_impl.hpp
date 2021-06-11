#pragma once

// IWYU pragma: private, include "modle_tools/tools.hpp"

#include <absl/container/flat_hash_map.h>  // for flat_hash_map, BitMask, operator!=
#include <absl/container/flat_hash_set.h>  // for flat_hash_set
#include <absl/meta/type_traits.h>         // for remove_reference_t
#include <absl/strings/match.h>            // for StartsWith
#include <absl/strings/str_cat.h>          // for StrCat
#include <absl/strings/str_join.h>         // for StrJoin
#include <absl/strings/string_view.h>      // for string_view
#include <absl/strings/strip.h>            // for StripPrefix, StripSuffix
#include <absl/time/clock.h>               // for Now
#include <absl/time/time.h>                // for FormatDuration, operator-, Time
#include <absl/types/span.h>               // for Span
#include <fmt/format.h>                    // for print, FMT_STRING, format
#include <fmt/ostream.h>                   // for print, formatbuf<>::int_type

#include <algorithm>                                // for fill, transform, max
#include <array>                                    // for array, array<>::value_type
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cassert>                                  // for assert
#include <cgranges/IITree.hpp>                      // for IITree
#include <cmath>                                    // for sqrt
#include <cstdint>                                  // for uint32_t
#include <cstdio>                                   // for stderr, stdout
#include <filesystem>                               // for path, operator<<, create_directories
#include <fstream>                                  // for size_t, ofstream, basic_ofstream, str...
#include <functional>                               // for ref
#include <initializer_list>                         // for initializer_list
#include <iterator>                                 // for insert_iterator, inserter
#include <memory>                                   // for unique_ptr, make_unique
#include <numeric>                                  // for iota
#include <sstream>                                  // for basic_stringbuf<>::int_type, basic_st...
#include <stdexcept>                                // for runtime_error
#include <string>                                   // for string, basic_string
#include <thread>                                   // for thread
#include <utility>                                  // for pair, make_pair, move
#include <vector>                                   // for vector

#include "modle/bed.hpp"            // for Parser
#include "modle/bigwig.hpp"         // for write_range, init_bigwig_file, close_...
#include "modle/contacts.hpp"       // for ContactMatrix
#include "modle/cooler.hpp"         // for Cooler, Cooler::READ_ONLY, Cooler::WR...
#include "modle/hdf5.hpp"           // for read_attribute_int
#include "modle_tools/config.hpp"   // for config
#include "modle_tools/eval.hpp"     // for Transformation, Cross, Linear, comput...
#include "modle_tools/noisify.hpp"  // for noisify_contacts
#include "modle_tools/stats.hpp"    // for compute_number_of_contacts_after_depl...

namespace modle::tools {

void eval_subcmd(const modle::tools::config& c) {
  assert(c.compute_spearman || c.compute_pearson);  // NOLINT
  const auto bin_size =
      static_cast<size_t>(hdf5::read_attribute_int(c.path_to_input_matrix.string(), "bin-size"));

  auto chrom_list =  // This cannot be made const
      select_chromosomes_for_eval(c.path_to_input_matrix.string(),
                                  c.path_to_reference_matrix.string(), bin_size);
  if (chrom_list.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Files {} and {} have 0 chromosomes in common. Make sure you are not trying "
                   "to compare different genome assemblies (chromosomes needs to have the same "
                   "name and size in order to qualify for comparison)"),
        c.path_to_input_matrix, c.path_to_reference_matrix));
  }

  absl::flat_hash_map<std::string, std::pair<size_t, size_t>> chrom_subranges;
  if (!c.path_to_chrom_subranges.empty()) {
    const auto records = bed::Parser(c.path_to_chrom_subranges).parse_all();
    chrom_subranges.reserve(records.size());
    std::transform(records.begin(), records.end(),
                   std::inserter(chrom_subranges, chrom_subranges.end()),
                   [](const auto& r) -> std::pair<std::string, std::pair<size_t, size_t>> {
                     return {r.chrom, {r.chrom_start, r.chrom_end}};
                   });
  }

  {
    std::vector<std::string_view> chrom_names(chrom_list.size());
    std::transform(chrom_list.begin(), chrom_list.end(), chrom_names.begin(), [](const auto& p) {
      return std::string_view{p.first.data(), p.first.size()};
    });
    if (chrom_list.size() == 1) {
      fmt::print(stderr, FMT_STRING("Computing correlation for chromosome: '{}'\n"),
                 chrom_names.front());
    } else {
      fmt::print(stderr,
                 FMT_STRING("Computing correlation for the following {} chromosomes: '{}'\n"),
                 chrom_list.size(), absl::StrJoin(chrom_names, "', '"));
    }
  }

  const auto bn = c.output_base_name.string();

  auto create_bwig_file = [&](std::string_view fname_suffix, bool skip) -> bigwig::file {
    if (skip) {
      return {nullptr, bigwig::close_bigwig_file};
    }
    return bigwig::init_bigwig_file(absl::StrCat(bn, "_", absl::StripPrefix(fname_suffix, "_")),
                                    chrom_list);
  };

  // Init files and write bw header
  auto bw_corr_linear_pearson = create_bwig_file("pearson_r_linear.bw", !c.compute_pearson);
  auto bw_pv_linear_pearson = create_bwig_file("pearson_pv_linear.bw", !c.compute_pearson);
  auto bw_corr_cross_pearson = create_bwig_file("pearson_r_cross.bw", !c.compute_pearson);
  auto bw_pv_cross_pearson = create_bwig_file("pearson_pv_cross.bw", !c.compute_pearson);
  auto bw_corr_linear_spearman = create_bwig_file("spearman_r_linear.bw", !c.compute_spearman);
  auto bw_pv_linear_spearman = create_bwig_file("spearman_pv_linear.bw", !c.compute_spearman);
  auto bw_corr_cross_spearman = create_bwig_file("spearman_r_cross.bw", !c.compute_spearman);
  auto bw_pv_cross_spearman = create_bwig_file("spearman_pv_cross.bw", !c.compute_spearman);
  auto bw_linear_sed = create_bwig_file("eucl_dist_linear.bw", !c.compute_edist);
  auto bw_cross_sed = create_bwig_file("eucl_dist_cross.bw", !c.compute_edist);

  auto ref_cooler = cooler::Cooler(c.path_to_reference_matrix, cooler::Cooler::READ_ONLY, bin_size);
  auto input_cooler = cooler::Cooler(c.path_to_input_matrix, cooler::Cooler::READ_ONLY, bin_size);
  const auto nrows = (c.diagonal_width / bin_size) + (c.diagonal_width % bin_size != 0);  // NOLINT
  assert(nrows != 0);                                                                     // NOLINT

  std::vector<double> pc_linear_corr_buff;
  std::vector<double> pc_linear_pval_buff;

  std::vector<double> pc_cross_corr_buff;
  std::vector<double> pc_cross_pval_buff;

  std::vector<double> sc_linear_corr_buff;
  std::vector<double> sc_linear_pval_buff;

  std::vector<double> ed_linear_buff;
  std::vector<double> ed_cross_buff;

  std::vector<double> sc_cross_corr_buff;
  std::vector<double> sc_cross_pval_buff;

  auto pcc = [&](std::string_view chrom_name, absl::Span<const uint32_t> v1,
                 absl::Span<const uint32_t> v2, size_t ncols, size_t offset,
                 std::vector<double>& corr_buff, std::vector<double>& pval_buff, Transformation t) {
    if (!c.compute_pearson) {
      return;
    }
    const auto t0 = absl::Now();
    compute_pearson_over_range(v1, v2, corr_buff, pval_buff, nrows, ncols, t);
    switch (t) {
      case Transformation::Linear:
        bigwig::write_range(std::string{chrom_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_linear_pearson);
        bigwig::write_range(std::string{chrom_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_linear_pearson);
        fmt::print(stderr,
                   FMT_STRING("Pearson \"linear\" correlation calculation "
                              "completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
      case Transformation::Cross:
        bigwig::write_range(std::string{chrom_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_cross_pearson);
        bigwig::write_range(std::string{chrom_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_cross_pearson);
        fmt::print(stderr,
                   FMT_STRING("Pearson \"cross\" correlation calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
    }
  };

  auto src = [&](std::string_view chrom_name, absl::Span<const uint32_t> v1,
                 absl::Span<const uint32_t> v2, size_t ncols, size_t offset,
                 std::vector<double>& corr_buff, std::vector<double>& pval_buff, Transformation t) {
    if (!c.compute_spearman) {
      return;
    }
    const auto t0 = absl::Now();
    compute_spearman_over_range(v1, v2, corr_buff, pval_buff, nrows, ncols, t);
    switch (t) {
      case Transformation::Linear:
        bigwig::write_range(std::string{chrom_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_linear_spearman);
        bigwig::write_range(std::string{chrom_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_linear_spearman);
        fmt::print(stderr,
                   FMT_STRING("Spearman \"linear\" correlation calculation "
                              "completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
      case Transformation::Cross:
        bigwig::write_range(std::string{chrom_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_cross_spearman);
        bigwig::write_range(std::string{chrom_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_cross_spearman);
        fmt::print(stderr,
                   FMT_STRING("Spearman \"cross\" correlation calculation "
                              "completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
    }
  };

  auto edist = [&](std::string_view chrom_name, absl::Span<const uint32_t> v1,
                   absl::Span<const uint32_t> v2, size_t ncols, size_t offset,
                   std::vector<double>& sed_buff, Transformation t) {
    if (!c.compute_edist) {
      return;
    }
    const auto t0 = absl::Now();
    compute_euc_dist_over_range(v1, v2, sed_buff, nrows, ncols, t);
    switch (t) {
      case Transformation::Linear:
        bigwig::write_range(std::string{chrom_name}, sed_buff, offset, bin_size, bin_size,
                            bw_linear_sed);
        fmt::print(stderr, FMT_STRING("Euclidean dist. \"linear\" calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
      case Transformation::Cross:
        bigwig::write_range(std::string{chrom_name}, sed_buff, offset, bin_size, bin_size,
                            bw_cross_sed);
        fmt::print(stderr, FMT_STRING("Euclidean dist. \"cross\" calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
    }
  };

  absl::flat_hash_map<std::string, size_t> ref_chrom_idxes;
  absl::flat_hash_map<std::string, size_t> inp_chrom_idxes;

  size_t i = 0;
  for (auto& name : ref_cooler.get_chrom_names()) {
    ref_chrom_idxes.emplace(std::move(name), i++);
  }
  i = 0;
  for (auto& name : input_cooler.get_chrom_names()) {
    inp_chrom_idxes.emplace(std::move(name), i++);
  }

  for (const auto& chrom : chrom_list) {
    std::string q;
    for (const auto& prefix : {"chrom", "chrom", "chrom"}) {
      if (absl::StartsWith(chrom.first, prefix)) {
        q = absl::StripPrefix(chrom.first, prefix);
      } else {
        q = absl::StrCat(prefix, chrom.first);
      }
      if (auto it = ref_chrom_idxes.find(q); it != ref_chrom_idxes.end()) {
        const auto idx = it->second;
        ref_chrom_idxes.emplace(chrom.first, idx);
        ref_chrom_idxes.erase(q);
      }
      if (auto it = inp_chrom_idxes.find(q); it != inp_chrom_idxes.end()) {
        const auto idx = it->second;
        inp_chrom_idxes.emplace(chrom.first, idx);
        inp_chrom_idxes.erase(q);
      }
    }
  }

  std::array<std::thread, 6> threads;  // NOLINT
  for (const auto& chrom : chrom_list) {
    const auto& chrom_name = chrom.first;
    auto chrom_subrange = std::make_pair(0UL, static_cast<size_t>(chrom.second));
    if (!chrom_subranges.empty()) {
      auto it = chrom_subranges.find(chrom_name);

      if (it != chrom_subranges.end()) {
        chrom_subrange = it->second;
      } else {  // Try common prefixes
        for (const auto& prefix : {"chrom", "chrom", "chrom"}) {
          if (absl::StartsWith(chrom_name, prefix)) {
            it = chrom_subranges.find(absl::StripPrefix(chrom_name, prefix));
          } else {
            it = chrom_subranges.find(absl::StrCat(prefix, chrom_name));
          }
          if (it != chrom_subranges.end()) {
            chrom_subrange = it->second;
            break;
          }
        }
      }
    }

    auto t0 = absl::Now();
    fmt::print(stderr, FMT_STRING("Reading contacts for '{}' into memory...\n"), chrom_name);
    if (!ref_cooler.has_contacts_for_chrom(ref_chrom_idxes.at(chrom_name))) {
      fmt::print(stderr,
                 FMT_STRING("WARNING: reference contact matrix doesn't have any contacts for "
                            "'{}'. SKIPPING!\n"),
                 chrom_name);
      continue;
    }

    if (!input_cooler.has_contacts_for_chrom(inp_chrom_idxes.at(chrom_name))) {
      fmt::print(stderr,
                 FMT_STRING("WARNING: contact matrix doesn't have any contacts "
                            "for '{}'. SKIPPING!\n"),
                 chrom_name);
      continue;
    }

    auto cmatrix1 = ref_cooler.cooler_to_cmatrix(chrom_name, nrows, chrom_subrange);
    if (c.deplete_contacts_from_reference) {
      cmatrix1.deplete_contacts(c.depletion_multiplier);
    }
    fmt::print(stderr,
               FMT_STRING("Read {:.2f}M contacts for a {}x{} reference matrix in {} using "
                          "{:.2f} MB of RAM.\n"),
               static_cast<double>(cmatrix1.get_tot_contacts()) / 1.0e6, cmatrix1.nrows(),
               cmatrix1.ncols(), absl::FormatDuration(absl::Now() - t0),
               cmatrix1.get_matrix_size_in_mb());
    t0 = absl::Now();

    const auto cmatrix2 = input_cooler.cooler_to_cmatrix(chrom_name, nrows, chrom_subrange);
    fmt::print(stderr,
               FMT_STRING("Read {:.2f}M contacts for a {}x{} input matrix in {} using "
                          "{:.2f} MB of RAM.\n"),
               static_cast<double>(cmatrix2.get_tot_contacts()) / 1.0e6, cmatrix2.nrows(),
               cmatrix2.ncols(), absl::FormatDuration(absl::Now() - t0),
               cmatrix2.get_matrix_size_in_mb());

    if (cmatrix1.ncols() != cmatrix2.ncols() || cmatrix1.nrows() != cmatrix2.nrows()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("An error occurred while computing the correlation for "
                                 "'{}' between files "
                                 "{} and {}: Contact matrices should have the same shape "
                                 "m1=[{}][{}], m2=[{}][{}]"),
                      chrom_name, c.path_to_reference_matrix, c.path_to_input_matrix,
                      cmatrix1.nrows(), cmatrix1.ncols(), cmatrix2.nrows(), cmatrix2.ncols()));
    }

    const auto ncols = cmatrix1.ncols();

    const auto& v1 = cmatrix1.get_raw_count_vector();
    const auto& v2 = cmatrix2.get_raw_count_vector();
    assert(v1.size() == v2.size());  // NOLINT

    fmt::print(stderr, FMT_STRING("Computing correlation(s) for '{}'...\n"), chrom_name);

    threads[0] = std::thread(pcc, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(pc_linear_corr_buff), std::ref(pc_linear_pval_buff),
                             Transformation::Linear);
    threads[1] = std::thread(pcc, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(pc_cross_corr_buff), std::ref(pc_cross_pval_buff),
                             Transformation::Cross);
    threads[2] = std::thread(src, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(sc_linear_corr_buff), std::ref(sc_linear_pval_buff),
                             Transformation::Linear);
    threads[3] = std::thread(src, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(sc_cross_corr_buff), std::ref(sc_cross_pval_buff),
                             Transformation::Cross);

    threads[4] = std::thread(edist, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(ed_linear_buff), Transformation::Linear);

    threads[5] = std::thread(edist, chrom_name, v1, v2, ncols, chrom_subrange.first,
                             std::ref(ed_cross_buff), Transformation::Cross);
    for (auto& t : threads) {
      t.join();
    }
  }
}

void filter_barriers(const modle::tools::config& c) {
  using IITree_t = IITree<uint32_t, uint32_t>;
  assert(!c.path_to_bed_files_for_filtering.empty());  // NOLINT
  absl::flat_hash_map<std::string, std::vector<IITree_t>> intervals(
      c.path_to_bed_files_for_filtering.size());

  for (auto i = 0UL; i < c.path_to_bed_files_for_filtering.size(); ++i) {
    for (const auto& record : bed::Parser(c.path_to_bed_files_for_filtering[i].string(),
                                          c.bed_dialect, c.strict_bed_validation)
                                  .parse_all(false)) {
      auto [node, new_insertion] = intervals.try_emplace(record.chrom, std::vector<IITree_t>{});
      if (new_insertion) {
        node->second.resize(c.path_to_bed_files_for_filtering.size());
      }
      node->second[i].add(record.chrom_start, record.chrom_end, uint8_t(0));
    }

    for (auto& [_, trees] : intervals) {
      trees[i].index();
    }
  }

  std::vector<size_t> buff;
  for (const auto& record :
       bed::Parser(c.path_to_extrusion_barrier_motifs_bed, c.bed_dialect, c.strict_bed_validation)
           .parse_all(false)) {
    if (auto node = intervals.find(record.chrom); node != intervals.end()) {
      auto print = true;
      for (const auto& tree : node->second) {
        if (tree.size() == 0) {
          print = false;
          break;
        }
        tree.overlap(record.chrom_start, record.chrom_end, buff);
        if (buff.empty()) {
          print = false;
          break;
        }
      }
      if (print) {
        fmt::print(stdout, "{}\n", record.to_string());
      }
    }
  }
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
    std::filesystem::remove_all(path);
    m2 = std::make_unique<cooler::Cooler>(path, cooler::Cooler::WRITE_ONLY, bin_size);
  }

  if (!path_to_output_hist.empty()) {  // Create hist. file
    std::filesystem::create_directories(std::filesystem::path(path_to_output_hist).parent_path());
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
      m2->write_or_append_cmatrix_to_file(cmatrix, chrom_name, 0L, chrom_size, chrom_size, true);
    }
  }

  fmt::print(
      stdout, FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\n"), "grand_average", tot_contacts,
      tot_contacts_after_depl,
      static_cast<double>(tot_contacts) / static_cast<double>(tot_number_of_pixels),
      static_cast<double>(tot_contacts_after_depl) / static_cast<double>(tot_number_of_pixels), 0);
}

void noisify_subcmd(const modle::tools::config& c) { modle::tools::noisify_contacts(c); }

}  // namespace modle::tools
