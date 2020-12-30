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

bigWigFile_t* init_bigwig_file(std::string_view output_path, std::vector<char*>& chr_names,
                               std::vector<uint32_t>& chr_sizes, int32_t zoom_levels,
                               std::size_t buff_size) {
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

bigWigFile_t* init_bigwig_file(std::string_view output_path, std::string& chr_name,
                               uint64_t chr_size, int32_t zoom_levels, std::size_t buff_size) {
  // Create the chromosome lists
  std::vector<char*> chr_names{chr_name.data()};
  std::vector<uint32_t> chr_sizes{static_cast<uint32_t>(chr_size)};
  return init_bigwig_file(output_path, chr_names, chr_sizes, zoom_levels, buff_size);
}

uint64_t write_range(const std::string& chr_name, uint64_t chr_length,
                     const std::vector<double>& vals, uint64_t offset, uint64_t span, uint64_t step,
                     std::string_view output_file) {
  absl::flat_hash_map<std::pair<std::string /* chr_name */, uint64_t /* chr_leng */>,
                      std::vector<double>>
      data{{{chr_name, chr_length}, vals}};
  return write_range(data, offset, span, step, output_file);
}

}  // namespace modle::bigwig