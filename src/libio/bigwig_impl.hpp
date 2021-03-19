#pragma once

// IWYU pragma: private, include "modle/bigwig.hpp"

#include <absl/container/flat_hash_map.h>  // for flat_hash_map
#include <absl/strings/match.h>            // for StartsWith
#include <fmt/format.h>                    // for format, FMT_STRING
#include <libBigWig/bigWig.h>              // for bwAddIntervalSpanSteps, bwCleanup, bwClose

#include <cassert>      // for assert
#include <cstddef>      // IWYU pragm: keep for size_t
#include <cstdint>      // for uint32_t, uint64_t, int32_t, int64_t
#include <filesystem>   // for file_size
#include <iterator>     // for pair
#include <limits>       // for numeric_limits
#include <memory>       // for unique_ptr
#include <stdexcept>    // for runtime_error
#include <string>       // for string
#include <string_view>  // for string_view
#include <type_traits>  // for is_arithmetic, is_signed
#include <vector>       // for vector, allocator

#include "bigwig_impl.hpp"  // for FMT_COMPILE_STRING

namespace modle::bigwig {

void close_bigwig_file(bigWigFile_t* fp) {
  if (fp) {
    bwClose(fp);
    bwCleanup();
  }
}

template <typename I>
bigwig::file init_bigwig_file(std::string_view output_path,
                              std::vector<std::pair<std::string, I>>& chromosomes,
                              int32_t zoom_levels, std::size_t buff_size) {
  static_assert(std::is_integral_v<I>, "I should be an integral numeric type.");
  std::vector<char*> chr_names(chromosomes.size());
  std::vector<uint32_t> chr_sizes(chromosomes.size());

  for (auto i = 0UL; i < chromosomes.size(); ++i) {
    chr_names[i] = const_cast<char*>(  // NOLINT this should be fine as long as libBigWig doesn't
        chromosomes[i].first.data());  // try to modify the data stored in the char*

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
                 uint64_t offset, uint64_t span, uint64_t step, bigwig::file& bigwig_fp) {
  static_assert(std::is_arithmetic<N1>::value, "N1 should be a numeric type.");
  static_assert(std::is_arithmetic<N2>::value, "N2 should be a numeric type.");
  for (const auto& [chr_data, values] : data) {
    std::vector<float> fvalues(values.begin(), values.end());
    // NOLINTNEXTLINE(readability-implicit-bool-conversion)
    if (bwAddIntervalSpanSteps(bigwig_fp.get(),
                               // Casting away constness should be fine as long as libBigWig
                               // doesn't try to modify the data stored in the char*
                               const_cast<char*>(chr_data.first.c_str()),  // NOLINT
                               static_cast<uint32_t>(offset), static_cast<uint32_t>(span),
                               static_cast<uint32_t>(step), fvalues.data(),
                               static_cast<uint32_t>(fvalues.size()))) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Failed to write data for chr '{}'"), chr_data.first));
    }
  }
}

template <typename N1, typename N2>
uint64_t write_range(absl::flat_hash_map<std::pair<std::string, N1>, std::vector<N2>>& data,
                     uint64_t offset, uint64_t span, uint64_t step, std::string_view output_file) {
  std::vector<std::pair<std::string, N1>> chromosomes(data.size());
  std::transform(data.begin(), data.end(), chromosomes.begin(),
                 [](const auto& p) { return p.first; });

  bigwig::file f{nullptr, &close_bigwig_file};
  try {
    f = init_bigwig_file(output_file, chromosomes);
    write_range(data, offset, span, step, f);
  } catch (const std::runtime_error& e) {
    if (absl::StartsWith(e.what(), "Failed to write data for chr '")) {
      throw std::runtime_error(fmt::format(FMT_STRING("{} to file '{}'"), e.what(), output_file));
    }
    throw;
  }
  return std::filesystem::file_size(output_file);
}

void write_range(const std::string& chr_name, const std::vector<double>& vals, uint64_t offset,
                 uint64_t span, uint64_t step, bigwig::file& bigwig_fp) {
  assert(bigwig_fp);
  std::vector<float> fvalues(vals.begin(), vals.end());
  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  if (bwAddIntervalSpanSteps(bigwig_fp.get(),
                             // this should be fine as long as libBigWig doesn't try to
                             // modify the data stored in the char*
                             const_cast<char*>(chr_name.c_str()),  // NOLINT
                             static_cast<uint32_t>(offset), static_cast<uint32_t>(span),
                             static_cast<uint32_t>(step), fvalues.data(),
                             static_cast<uint32_t>(fvalues.size()))) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to write data for chr '{}'"), chr_name));
  }
}

bigwig::file init_bigwig_file(std::string_view output_path, std::vector<char*>& chr_names,
                              std::vector<uint32_t>& chr_sizes, int32_t zoom_levels,
                              std::size_t buff_size) {
  bigwig::file bw_fp{nullptr, &close_bigwig_file};
  if (bwInit(buff_size)) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(
        fmt::format("Failed to initialize a buffer for file: '{}'", output_path));
  }

  std::string output_path_tmp = output_path.data();  // This is just a workaround, because bwOpen
  // takes a char* ... not a const char*
  bw_fp.reset(bwOpen(output_path_tmp.data(), nullptr, "w"));
  if (!bw_fp) {
    throw std::runtime_error(
        fmt::format("An error occurred while opening file '{}' for writing", output_path));
  }

  if (bwCreateHdr(bw_fp.get(), zoom_levels)) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(
        fmt::format("Failed to initialize the file header for file '{}'", output_path));
  }

  bw_fp->cl =
      bwCreateChromList(chr_names.data(), chr_sizes.data(), static_cast<int64_t>(chr_sizes.size()));
  if (!bw_fp->cl) {
    throw std::runtime_error(
        fmt::format("Failed to create the chromosome list for file '{}'", output_path));
  }

  if (bwWriteHdr(bw_fp.get())) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(fmt::format("Failed to write file header to file '{}'", output_path));
  }
  return bw_fp;
}

bigwig::file init_bigwig_file(std::string_view output_path, std::string& chr_name,
                              uint64_t chr_size, int32_t zoom_levels, std::size_t buff_size) {
  // Create the chromosome lists
  std::vector<char*> chr_names{chr_name.data()};
  std::vector<uint32_t> chr_sizes{static_cast<uint32_t>(chr_size)};
  return init_bigwig_file(output_path, chr_names, chr_sizes, zoom_levels, buff_size);
}

void close_bigwig_file(bigwig::file fp) {
  if (fp) {
    bwClose(fp.get());
    bwCleanup();
  }
}

}  // namespace modle::bigwig
