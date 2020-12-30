#pragma once
#include <absl/strings/match.h>
#include <fmt/format.h>
#include <libBigWig/bigWig.h>

#include <cstdint>
#include <filesystem>
#include <limits>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace modle::bigwig {

template <typename I>
bigWigFile_t* init_bigwig_file(std::string_view output_path,
                               std::vector<std::pair<std::string, I>>& chromosomes,
                               int32_t zoom_levels, std::size_t buff_size) {
  static_assert(std::is_integral<I>::value, "I should be an integral numeric type.");
  std::vector<char*> chr_names(chromosomes.size());
  std::vector<uint32_t> chr_sizes(chromosomes.size());

  for (auto i = 0UL; i < chromosomes.size(); ++i) {
    chr_names[i] = const_cast<char*>(
        chromosomes[i].first.data());  // NOLINT this should be fine as long as libBigWig doesn't
                                       // try to modify the data stored in the char*
    if constexpr (const auto max_val = std::numeric_limits<uint32_t>::max();
                  std::numeric_limits<I>::max() > max_val) {
      if (chromosomes[i].second > max_val) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("We currently don't support writing chromosomes longer than ~4.29 Gbp (2^32 "
                       "bp) to bigWig files. '{}' has a length of {:.4f} Gbp"),
            chr_names[i], static_cast<double>(chromosomes[i].second) / 1e9));  // NOLINT
      }
    }
    if constexpr (std::is_signed<I>()) {
      if (chromosomes[i].second < 0) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("'{}' appears to have a negative length of {} bp. This is not allowed"),
            chr_names[i], chromosomes[i].second));
      }
    }
    chr_sizes[i] = static_cast<uint32_t>(chromosomes[i].second);
  }

  return (init_bigwig_file(output_path, chr_names, chr_sizes, zoom_levels, buff_size));
}

template <typename N1, typename N2>
void write_range(absl::flat_hash_map<std::pair<std::string, N1>, std::vector<N2>>& data,
                 uint64_t offset, uint64_t span, uint64_t step, bigWigFile_t* bigwig_fp) {
  static_assert(std::is_arithmetic<N1>::value, "N1 should be a numeric type.");
  static_assert(std::is_arithmetic<N2>::value, "N2 should be a numeric type.");
  try {
    for (const auto& [chr_data, values] : data) {
      std::vector<float> fvalues(values.begin(), values.end());
      if (bwAddIntervalSpanSteps(
              bigwig_fp,
              const_cast<char*>(
                  chr_data.first.c_str()),  // NOLINT this should be fine as long as libBigWig
                                            // doesn't try to modify the data stored in the char*
              static_cast<uint32_t>(offset), static_cast<uint32_t>(span),
              static_cast<uint32_t>(step), fvalues.data(), static_cast<uint32_t>(fvalues.size()))) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Failed to write data for chr '{}'"), chr_data.first));
      }
    }
    // Closing the file causes the zoom levels to be created
    bwClose(bigwig_fp);
    bwCleanup();
  } catch (...) {
    bwClose(bigwig_fp);
    bwCleanup();
    throw;
  }
}

template <typename N1, typename N2>
uint64_t write_range(absl::flat_hash_map<std::pair<std::string, N1>, std::vector<N2>>& data,
                     uint64_t offset, uint64_t span, uint64_t step, std::string_view output_file) {
  std::vector<std::pair<std::string, N1>> chromosomes(data.size());

  std::transform(data.begin(), data.end(), chromosomes.begin(),
                 [](const auto& p) { return p.first; });
  bigWigFile_t* f{nullptr};
  try {
    f = init_bigwig_file(output_file, chromosomes);
    write_range(data, offset, span, step, f);
  } catch (const std::runtime_error& e) {
    bwClose(f);
    bwCleanup();
    if (absl::StartsWith(e.what(), "Failed to write data for chr '")) {
      throw std::runtime_error(fmt::format(FMT_STRING("{} to file '{}'"), e.what(), output_file));
    }
    throw;
  }
  return std::filesystem::file_size(output_file);
}

}  // namespace modle::bigwig