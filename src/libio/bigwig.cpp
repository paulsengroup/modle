#include "modle/bigwig.hpp"  // for write_range

#include <absl/container/flat_hash_map.h>
#include <fmt/format.h>
#include <libBigWig/bigWig.h>  // for bwCleanup, bwClose, bigWigFile_t, bwAddIntervalSpanSteps, bwCreateChromList, bwCreateHdr, bwInit, bwOpen, bwWriteHdr

#include <cstdint>     // for uint*_t, int*_t
#include <filesystem>  // for file_size
#include <iosfwd>      // for size_t
#include <stdexcept>   // for runtime_error
#include <string_view>
#include <type_traits>  // for is_arithmetic
#include <utility>      // for pair, forward, move
#include <vector>

namespace modle::bigwig {

[[nodiscard]] bigWigFile_t* init_bigwig(std::string_view output_path, std::vector<char*>& chr_names,
                                        std::vector<uint32_t>& chr_sizes, int32_t zoom_levels = 10,
                                        std::size_t buff_size = 1U << 17U /* 128KiB */) {
  bigWigFile_t* bw_fp = nullptr;
  try {
    if (bwInit(buff_size)) {  // NOLINT(readability-implicit-bool-conversion)
      throw std::runtime_error(
          fmt::format("Failed to initialize a buffer for file: '{}'", output_path));
    }

    std::string output_path_tmp = output_path.data();  // This is just a workaround, because bwOpen
                                                       // takes a char* ... not a const char*
    bw_fp = bwOpen(output_path_tmp.data(), nullptr, "w");
    if (!bw_fp) {
      throw std::runtime_error(
          fmt::format("An error occurred while opening file '{}' for writing", output_path));
    }

    if (bwCreateHdr(bw_fp, zoom_levels)) {  // NOLINT(readability-implicit-bool-conversion)
      throw std::runtime_error(
          fmt::format("Failed to initialize the file header for file '{}'", output_path));
    }

    bw_fp->cl = bwCreateChromList(chr_names.data(), chr_sizes.data(),
                                  static_cast<int64_t>(chr_sizes.size()));
    if (!bw_fp->cl) {
      throw std::runtime_error(
          fmt::format("Failed to create the chromosome list for file '{}'", output_path));
    }

    if (bwWriteHdr(bw_fp)) {  // NOLINT(readability-implicit-bool-conversion)
      throw std::runtime_error(
          fmt::format("Failed to write file header to file '{}'", output_path));
    }
  } catch (...) {
    bwClose(bw_fp);
    bwCleanup();
    throw;
  }
  return bw_fp;
}

[[nodiscard]] bigWigFile_t* init_bigwig(std::string_view output_path, std::string& chr_name,
                                        uint64_t chr_size, int32_t zoom_levels = 10,
                                        std::size_t buff_size = 1U << 17U /* 128KiB */) {
  // Create the chromosome lists
  std::vector<char*> chr_names{chr_name.data()};
  std::vector<uint32_t> chr_sizes{static_cast<uint32_t>(chr_size)};
  return init_bigwig(output_path, chr_names, chr_sizes, zoom_levels, buff_size);
}

template <typename N>
uint64_t write_range(
    absl::flat_hash_map<std::pair<std::string, N>, const std::vector<double>*>& data,
    uint64_t offset, uint64_t span, uint64_t step, std::string_view output_file) {
  static_assert(std::is_arithmetic<N>::value, "N should be a numeric type.");
  std::vector<char*> chr_names(data.size());
  std::vector<uint32_t> chr_lengths(data.size());
  bigWigFile_t* bw_fp = nullptr;
  std::size_t i = 0;
  for (auto& [chr_data, _] : data) {
    chr_names[i] =
        const_cast<char*>(chr_data.first.c_str());  // Figure out how to avoid using const_cast
    chr_lengths[i++] = static_cast<uint32_t>(chr_data.second);
  }

  try {
    bw_fp = init_bigwig(output_file, chr_names, chr_lengths);
    for (const auto& [chr_data, values] : data) {
      std::vector<float> fvalues(values->begin(), values->end());
      if (bwAddIntervalSpanSteps(bw_fp, const_cast<char*>(chr_data.first.c_str()),
                                 static_cast<uint32_t>(offset), static_cast<uint32_t>(span),
                                 static_cast<uint32_t>(step), fvalues.data(),
                                 static_cast<uint32_t>(fvalues.size()))) {
        throw std::runtime_error(fmt::format("Failed to write data for chr '{}' to file '{}'",
                                             chr_data.first, output_file));
      }
    }
    // Closing the file causes the zoom levels to be created
    bwClose(bw_fp);
    bwCleanup();
  } catch (...) {
    bwClose(bw_fp);
    bwCleanup();
    throw;
  }
  return std::filesystem::file_size(output_file);
}

uint64_t write_range(const std::string& chr_name, uint64_t chr_length,
                     const std::vector<double>& vals, uint64_t offset, uint64_t span, uint64_t step,
                     std::string_view output_file) {
  absl::flat_hash_map<std::pair<std::string /* chr_name */, uint64_t /* chr_leng */>,
                      const std::vector<double>*>
      data{{{chr_name, chr_length}, &vals}};
  return write_range(data, offset, span, step, output_file);
}

}  // namespace modle::bigwig