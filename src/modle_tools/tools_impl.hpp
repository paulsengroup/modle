#pragma once

#include <absl/strings/match.h>  // for EndsWith
#include <absl/time/clock.h>     // for Now
#include <absl/time/time.h>      // for FormatDuration, operator-, Time
#include <fmt/format.h>          // for format, print

#include <array>
#include <boost/asio/thread_pool.hpp>
#include <cassert>
#include <cstdint>     // for uint32_t
#include <cstdio>      // for stderr
#include <filesystem>  // for create_directories
#include <mutex>

#include "modle/bed.hpp"
#include "modle/bigwig.hpp"    // for write_range
#include "modle/contacts.hpp"  // for ContactMatrix<>::Header, ContactMatrix
#include "modle/cooler.hpp"
#include "modle/hdf5.hpp"
#include "modle_tools/config.hpp"   // for config
#include "modle_tools/eval.hpp"     // for Transformation, Cross, Linear
#include "modle_tools/stats.hpp"

namespace modle::tools {

void eval_subcmd(const modle::tools::config& c) {
  assert(c.compute_spearman || c.compute_pearson);  // NOLINT
  assert(c.path_to_input_matrices.size() == 1);     // NOLINT
  const auto& path_to_input_cmatrix = c.path_to_input_matrices.front();

  const auto bin_size =
      static_cast<std::size_t>(hdf5::read_attribute_int(path_to_input_cmatrix, "bin-size"));

  auto chr_list =  // This cannot be made const
      select_chromosomes_for_eval(path_to_input_cmatrix, c.path_to_reference_matrix, bin_size);
  if (chr_list.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Files '{}' and '{}' have 0 chromosomes in common. Make sure you are not trying "
                   "to compare different genome assemblies (chromosomes needs to have the same "
                   "name and size in order to qualify for comparison)"),
        path_to_input_cmatrix, c.path_to_reference_matrix));
  }

  absl::flat_hash_map<std::string, std::pair<std::size_t, std::size_t>> chr_subranges;
  if (!c.path_to_chr_subranges.empty()) {
    const auto records = bed::Parser(c.path_to_chr_subranges).parse_all();
    chr_subranges.reserve(records.size());
    std::transform(
        records.begin(), records.end(), std::inserter(chr_subranges, chr_subranges.end()),
        [](const auto& r) -> std::pair<std::string, std::pair<std::size_t, std::size_t>> {
          return {r.chrom, {r.chrom_start, r.chrom_end}};
        });
  }

  {
    std::vector<std::string_view> chr_names(chr_list.size());
    std::transform(chr_list.begin(), chr_list.end(), chr_names.begin(), [](const auto& p) {
      return std::string_view{p.first.data(), p.first.size()};
    });
    if (chr_list.size() == 1) {
      fmt::print(stderr, FMT_STRING("Computing correlation for chromosome: '{}'\n"),
                 chr_names.front());
    } else {
      fmt::print(stderr,
                 FMT_STRING("Computing correlation for the following {} chromosomes: '{}'\n"),
                 chr_list.size(), absl::StrJoin(chr_names, "', '"));
    }
  }

  c.compute_pearson ? fmt::format(FMT_STRING("{}_spearman_linear_r.bw"), c.output_base_name) : "";
  const auto out_path_pv_linear_spearman =
      c.compute_pearson ? fmt::format(FMT_STRING("{}_spearman_linear_pv.bw"), c.output_base_name)
                        : "";
  const auto out_path_corr_cross_spearman =
      c.compute_pearson ? fmt::format(FMT_STRING("{}_spearman_cross_r.bw"), c.output_base_name)
                        : "";
  const auto out_path_pv_cross_spearman =
      c.compute_pearson ? fmt::format(FMT_STRING("{}_spearman_cross_pv.bw"), c.output_base_name)
                        : "";
  const auto& bn = c.output_base_name;
  // Init files and write bw header
  auto* bw_corr_linear_pearson =
      c.compute_pearson
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_pearson_linear_r.bw"), chr_list)
          : nullptr;
  auto* bw_pv_linear_pearson =
      c.compute_pearson
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_pearson_linear_pv.bw"), chr_list)
          : nullptr;
  auto* bw_corr_cross_pearson =
      c.compute_pearson
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_pearson_cross_r.bw"), chr_list)
          : nullptr;
  auto* bw_pv_cross_pearson =
      c.compute_pearson
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_pearson_cross_pv.bw"), chr_list)
          : nullptr;
  auto* bw_corr_linear_spearman =
      c.compute_spearman
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_spearman_linear_r.bw"), chr_list)
          : nullptr;
  auto* bw_pv_linear_spearman =
      c.compute_spearman
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_spearman_linear_pv.bw"), chr_list)
          : nullptr;
  auto* bw_corr_cross_spearman =
      c.compute_spearman
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_spearman_cross_r.bw"), chr_list)
          : nullptr;
  auto* bw_pv_cross_spearman =
      c.compute_spearman
          ? bigwig::init_bigwig_file(absl::StrCat(bn, "_spearman_cross_pv.bw"), chr_list)
          : nullptr;

  auto ref_cooler = cooler::Cooler(c.path_to_reference_matrix, cooler::Cooler::READ_ONLY, bin_size);
  auto input_cooler = cooler::Cooler(path_to_input_cmatrix, cooler::Cooler::READ_ONLY, bin_size);
  const auto nrows = (c.diagonal_width / bin_size) + (c.diagonal_width % bin_size != 0);  // NOLINT
  assert(nrows != 0);                                                                     // NOLINT

  std::vector<double> pc_linear_corr_buff;
  std::vector<double> pc_linear_pval_buff;

  std::vector<double> pc_cross_corr_buff;
  std::vector<double> pc_cross_pval_buff;

  std::vector<double> sc_linear_corr_buff;
  std::vector<double> sc_linear_pval_buff;

  std::vector<double> sc_cross_corr_buff;
  std::vector<double> sc_cross_pval_buff;

  auto pcc = [&](std::string_view chr_name, absl::Span<const uint32_t> v1,
                 absl::Span<const uint32_t> v2, std::size_t ncols, std::size_t offset,
                 std::vector<double>& corr_buff, std::vector<double>& pval_buff, Transformation t) {
    if (!c.compute_pearson) {
      return;
    }
    const auto t0 = absl::Now();
    compute_pearson_over_range(v1, v2, corr_buff, pval_buff, nrows, ncols, t);
    switch (t) {
      case Transformation::Linear:
        bigwig::write_range(std::string{chr_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_linear_pearson);
        bigwig::write_range(std::string{chr_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_linear_pearson);
        fmt::print(stderr,
                   FMT_STRING("Pearson \"linear\" correlation calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
      case Transformation::Cross:
        bigwig::write_range(std::string{chr_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_cross_pearson);
        bigwig::write_range(std::string{chr_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_cross_pearson);
        fmt::print(stderr,
                   FMT_STRING("Pearson \"cross\" correlation calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
    }
  };

  auto src = [&](std::string_view chr_name, absl::Span<const uint32_t> v1,
                 absl::Span<const uint32_t> v2, std::size_t ncols, std::size_t offset,
                 std::vector<double>& corr_buff, std::vector<double>& pval_buff, Transformation t) {
    if (!c.compute_spearman) {
      return;
    }
    const auto t0 = absl::Now();
    compute_spearman_over_range(v1, v2, corr_buff, pval_buff, nrows, ncols, t);
    switch (t) {
      case Transformation::Linear:
        bigwig::write_range(std::string{chr_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_linear_spearman);
        bigwig::write_range(std::string{chr_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_linear_spearman);
        fmt::print(stderr,
                   FMT_STRING("Spearman \"linear\" correlation calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
      case Transformation::Cross:
        bigwig::write_range(std::string{chr_name}, corr_buff, offset, bin_size, bin_size,
                            bw_corr_cross_spearman);
        bigwig::write_range(std::string{chr_name}, pval_buff, offset, bin_size, bin_size,
                            bw_pv_cross_spearman);
        fmt::print(stderr,
                   FMT_STRING("Spearman \"cross\" correlation calculation completed in {}.\n"),
                   absl::FormatDuration(absl::Now() - t0));
        break;
    }
  };

  absl::flat_hash_map<std::string, std::size_t> ref_chr_idxes;
  absl::flat_hash_map<std::string, std::size_t> inp_chr_idxes;

  std::size_t i = 0;
  for (auto& name : ref_cooler.get_chr_names()) {
    ref_chr_idxes.emplace(std::move(name), i++);
  }
  i = 0;
  for (auto& name : input_cooler.get_chr_names()) {
    inp_chr_idxes.emplace(std::move(name), i++);
  }

  for (const auto& chr : chr_list) {
    std::string_view q;
    for (const auto& prefix : {"chr", "CHR", "Chr"}) {
      if (absl::StartsWith(chr.first, prefix)) {
        q = absl::StripPrefix(chr.first, prefix);
      } else {
        q = absl::StrCat(prefix, chr.first);
      }
      if (auto it = ref_chr_idxes.find(q); it != ref_chr_idxes.end()) {
        ref_chr_idxes.emplace(chr.first, it->second);
        ref_chr_idxes.erase(q);
      }
      if (auto it = inp_chr_idxes.find(q); it != inp_chr_idxes.end()) {
        inp_chr_idxes.emplace(chr.first, it->second);
        inp_chr_idxes.erase(q);
      }
    }
  }

  std::array<std::thread, 4> threads;
  for (const auto& chr : chr_list) {
    const auto& chr_name = chr.first;
    auto chr_subrange = std::make_pair(0UL, static_cast<std::size_t>(chr.second));
    if (!chr_subranges.empty()) {
      auto it = chr_subranges.find(chr_name);

      if (it != chr_subranges.end()) {
        chr_subrange = it->second;
      } else {  // Try common prefixes
        for (const auto& prefix : {"chr", "CHR", "Chr"}) {
          if (absl::StartsWith(chr_name, prefix)) {
            it = chr_subranges.find(absl::StripPrefix(chr_name, prefix));
          } else {
            it = chr_subranges.find(absl::StrCat(prefix, chr_name));
          }
          if (it != chr_subranges.end()) {
            chr_subrange = it->second;
            break;
          }
        }
      }
    }

    auto t0 = absl::Now();
    fmt::print(stderr, FMT_STRING("Reading contacts for '{}' into memory...\n"), chr_name);
    if (!ref_cooler.has_contacts_for_chr(ref_chr_idxes.at(chr_name))) {
      fmt::print(stderr,
                 FMT_STRING("WARNING: reference contact matrix doesn't have any contacts for "
                            "'{}'. SKIPPING!\n"),
                 chr_name);
      continue;
    }

    if (!input_cooler.has_contacts_for_chr(inp_chr_idxes.at(chr_name))) {
      fmt::print(
          stderr,
          FMT_STRING("WARNING: contact matrix doesn't have any contacts for '{}'. SKIPPING!\n"),
          chr_name);
      continue;
    }

    auto cmatrix1 = ref_cooler.cooler_to_cmatrix(chr_name, nrows, chr_subrange);
    if (c.deplete_contacts_from_reference) {
      cmatrix1.deplete_contacts(c.depletion_multiplier);
    }
    fmt::print(
        stderr,
        FMT_STRING(
            "Read {:.2f}M contacts for a {}x{} reference matrix in {} using {:.2f} MB of RAM.\n"),
        static_cast<double>(cmatrix1.get_tot_contacts()) / 1.0e6, cmatrix1.nrows(),
        cmatrix1.ncols(), absl::FormatDuration(absl::Now() - t0), cmatrix1.get_matrix_size_in_mb());
    t0 = absl::Now();

    const auto cmatrix2 = input_cooler.cooler_to_cmatrix(chr_name, nrows, chr_subrange);
    fmt::print(
        stderr,
        FMT_STRING(
            "Read {:.2f}M contacts for a {}x{} input matrix in {} using {:.2f} MB of RAM.\n"),
        static_cast<double>(cmatrix2.get_tot_contacts()) / 1.0e6, cmatrix2.nrows(),
        cmatrix2.ncols(), absl::FormatDuration(absl::Now() - t0), cmatrix2.get_matrix_size_in_mb());

    if (cmatrix1.ncols() != cmatrix2.ncols() || cmatrix1.nrows() != cmatrix2.nrows()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("An error occurred while computing the correlation for '{}' between files "
                     "'{}' and '{}': Contact matrices should have the same shape "
                     "m1=[{}][{}], m2=[{}][{}]"),
          chr_name, c.path_to_reference_matrix, path_to_input_cmatrix, cmatrix1.nrows(),
          cmatrix1.ncols(), cmatrix2.nrows(), cmatrix2.ncols()));
    }

    const auto ncols = cmatrix1.ncols();

    const auto& v1 = cmatrix1.get_raw_count_vector();
    const auto& v2 = cmatrix2.get_raw_count_vector();
    assert(v1.size() == v2.size());  // NOLINT

    fmt::print(stderr, FMT_STRING("Computing correlation(s) for '{}'...\n"), chr_name);

    threads[0] =
        std::thread(pcc, chr_name, v1, v2, ncols, chr_subrange.first, std::ref(pc_linear_corr_buff),
                    std::ref(pc_linear_pval_buff), Transformation::Linear);
    threads[1] =
        std::thread(pcc, chr_name, v1, v2, ncols, chr_subrange.first, std::ref(pc_cross_corr_buff),
                    std::ref(pc_cross_pval_buff), Transformation::Cross);
    threads[2] =
        std::thread(src, chr_name, v1, v2, ncols, chr_subrange.first, std::ref(sc_linear_corr_buff),
                    std::ref(sc_linear_pval_buff), Transformation::Linear);
    threads[3] =
        std::thread(src, chr_name, v1, v2, ncols, chr_subrange.first, std::ref(sc_cross_corr_buff),
                    std::ref(sc_cross_pval_buff), Transformation::Cross);
    for (auto& t : threads) {
      t.join();
    }
  }

  bigwig::close_bigwig_file(bw_corr_linear_pearson);
  bigwig::close_bigwig_file(bw_pv_linear_pearson);
  bigwig::close_bigwig_file(bw_corr_cross_pearson);
  bigwig::close_bigwig_file(bw_pv_cross_pearson);
  bigwig::close_bigwig_file(bw_corr_linear_spearman);
  bigwig::close_bigwig_file(bw_pv_linear_spearman);
  bigwig::close_bigwig_file(bw_corr_cross_spearman);
  bigwig::close_bigwig_file(bw_pv_cross_spearman);
}

void stats_subcmd(const modle::tools::config& c) {
  assert(c.path_to_input_matrices.size() == 1);  // NOLINT
  const auto& path_to_input_matrix = c.path_to_input_matrices.front();
  const auto& path_to_output_hist = c.output_path_for_histograms;

  cooler::Cooler m1(path_to_input_matrix, cooler::Cooler::READ_ONLY, c.bin_size);

  std::unique_ptr<cooler::Cooler> m2{nullptr};
  std::unique_ptr<std::ofstream> hist_file{nullptr};

  const auto nchroms = m1.get_nchroms();
  const auto bin_size = m1.get_bin_size();
  const auto chr_names = m1.get_chr_names();
  const auto chr_sizes = m1.get_chr_sizes();

  absl::flat_hash_map<std::string, std::size_t> chr_start_offsets;
  if (!c.path_to_chr_subranges.empty()) {
    modle::bed::Parser p(c.path_to_chr_subranges);
    const auto buff = p.parse_all();
    std::transform(
        buff.begin(), buff.end(), std::inserter(chr_start_offsets, chr_start_offsets.end()),
        [](const auto& record) { return std::make_pair(record.chrom, record.chrom_start); });
  }

  if (c.dump_depleted_matrices) {  // Create cooler file to write depl. contacts
    const auto ext = std::filesystem::path(path_to_input_matrix).extension().string();
    const auto path = absl::StrCat(absl::StripSuffix(path_to_input_matrix, ext), "_depl.cool");
    std::filesystem::remove_all(path);
    m2 = std::make_unique<cooler::Cooler>(path, cooler::Cooler::WRITE_ONLY, bin_size);
  }

  if (!path_to_output_hist.empty()) {  // Create hist. file
    std::filesystem::create_directories(std::filesystem::path(path_to_output_hist).parent_path());
    hist_file = std::make_unique<std::ofstream>(path_to_output_hist);
    std::vector<std::size_t> buff(c.diagonal_width / bin_size);
    std::iota(buff.begin(), buff.end(), 0);
    // TODO: Consider whether it make sense to write this header
    fmt::print(*hist_file, "#{}\n", absl::StrJoin(buff, "\t"));
  }

  // Write header
  fmt::print(stdout,
             "chr_name\ttot_number_of_contacts_1\ttot_number_of_contacts_2\tavg_number_of_contacts_"
             "1\tavg_number_of_contacts_2\tfraction_of_graylisted_bins\n");

  std::size_t tot_contacts = 0;
  std::size_t tot_contacts_after_depl = 0;
  std::size_t tot_number_of_pixels = 0;

  for (auto i = 0UL; i < nchroms; ++i) {
    const auto& chr_name = chr_names[i];
    const auto& chr_size = chr_sizes[i];

    // Skip chromosome found in the exclusion list
    if (c.chromosomes_excluded.contains(chr_name)) {
      continue;
    }

    // Read contacts for chr_name into memory
    auto cmatrix = m1.cooler_to_cmatrix(chr_name, c.diagonal_width, bin_size);

    const auto hist = cmatrix.compute_row_wise_contact_histogram();
    const auto mask = cmatrix.generate_mask_for_bins_without_contacts();

    const auto chr_contacts = cmatrix.get_tot_contacts();
    const auto chr_contacts_after_depl = compute_number_of_contacts_after_depletion(
        cmatrix, hist, mask.count(), c.depletion_multiplier);
    assert(chr_contacts_after_depl <= chr_contacts);  // NOLINT

    const auto chr_avg_contacts =
        static_cast<double>(chr_contacts) / static_cast<double>(cmatrix.npixels_after_masking());
    const auto chr_avg_contacts_after_depl = static_cast<double>(chr_contacts_after_depl) /
                                             static_cast<double>(cmatrix.npixels_after_masking());
    const auto fraction_of_graylisted_bins =
        1.0 - (static_cast<double>(cmatrix.npixels_after_masking()) /  // NOLINT
               static_cast<double>(cmatrix.npixels()));

    tot_contacts += chr_contacts;
    tot_contacts_after_depl += chr_contacts_after_depl;
    tot_number_of_pixels += cmatrix.npixels_after_masking();

    // clang-format off
    fmt::print(stdout, FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\n"),
               chr_name,
               chr_contacts,
               chr_contacts_after_depl,
               chr_avg_contacts,
               chr_avg_contacts_after_depl,
               fraction_of_graylisted_bins);
    // clang-format on

    if (hist_file) {
      fmt::print(*hist_file, "{}\n", absl::StrJoin(hist, "\t"));
    }

    if (m2) {
      cmatrix.deplete_contacts(c.depletion_multiplier);
      m2->write_or_append_cmatrix_to_file(cmatrix, chr_name, 0L, chr_size, chr_size, true);
    }
  }

  fmt::print(
      stdout, FMT_STRING("{}\t{}\t{}\t{}\t{}\t{}\n"), "grand_average", tot_contacts,
      tot_contacts_after_depl,
      static_cast<double>(tot_contacts) / static_cast<double>(tot_number_of_pixels),
      static_cast<double>(tot_contacts_after_depl) / static_cast<double>(tot_number_of_pixels), 0);
}

}  // namespace modle::tools
