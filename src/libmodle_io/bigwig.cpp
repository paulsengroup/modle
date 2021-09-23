// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/bigwig.hpp"

#include <fmt/format.h>  // for format, FMT_STRING

#include <cassert>      // for assert
#include <cstddef>      // for size_t
#include <cstdint>      // for uint32_t, uint64_t, int32_t, int64_t
#include <memory>       // for allocator, unique_ptr
#include <stdexcept>    // for runtime_error
#include <string>       // for string
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "libBigWig/bigWig.h"  // for bwCleanup, bwClose, bwAddIntervalSpanSteps, bwCreateChromList

namespace modle::bigwig {

void close_bigwig_file(bigWigFile_t* fp) {
  if (fp) {
    bwClose(fp);
    bwCleanup();
  }
}

void write_range(const std::string& chrom_name, const std::vector<double>& vals, uint64_t offset,
                 uint64_t span, uint64_t step, bigwig::file& bigwig_fp) {
  assert(bigwig_fp);  // NOLINT
  std::vector<float> fvalues(vals.begin(), vals.end());
  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  if (bwAddIntervalSpanSteps(bigwig_fp.get(),
                             // this should be fine as long as libBigWig doesn't try to
                             // modify the data stored in the char*
                             const_cast<char*>(chrom_name.c_str()),  // NOLINT
                             static_cast<uint32_t>(offset), static_cast<uint32_t>(span),
                             static_cast<uint32_t>(step), fvalues.data(),
                             static_cast<uint32_t>(fvalues.size()))) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to write data for chrom '{}'"), chrom_name));
  }
}

bigwig::file init_bigwig_file(std::string_view output_path, std::vector<char*>& chrom_names,
                              std::vector<uint32_t>& chrom_sizes, int32_t zoom_levels,
                              size_t buff_size) {
  bigwig::file bw_fp{nullptr, &close_bigwig_file};
  if (bwInit(buff_size)) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to initialize a buffer for file: '{}'"), output_path));
  }

  std::string output_path_tmp = output_path.data();  // This is just a workaround, because bwOpen
  // takes a char* ... not a const char*
  bw_fp.reset(bwOpen(output_path_tmp.data(), nullptr, "w"));
  if (!bw_fp) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("An error occurred while opening file '{}' for writing"), output_path));
  }

  if (bwCreateHdr(bw_fp.get(), zoom_levels)) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to initialize the file header for file '{}'"), output_path));
  }

  bw_fp->cl = bwCreateChromList(chrom_names.data(), chrom_sizes.data(),
                                static_cast<int64_t>(chrom_sizes.size()));
  if (!bw_fp->cl) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to create the chromosome list for file '{}'"), output_path));
  }

  if (bwWriteHdr(bw_fp.get())) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to write file header to file '{}'"), output_path));
  }
  return bw_fp;
}

bigwig::file init_bigwig_file(std::string_view output_path, std::string& chrom_name,
                              uint64_t chrom_size, int32_t zoom_levels, size_t buff_size) {
  // Create the chromosome lists
  std::vector<char*> chrom_names{chrom_name.data()};
  std::vector<uint32_t> chrom_sizes{static_cast<uint32_t>(chrom_size)};
  return init_bigwig_file(output_path, chrom_names, chrom_sizes, zoom_levels, buff_size);
}

void close_bigwig_file(bigwig::file fp) {
  if (fp) {
    bwClose(fp.get());
    bwCleanup();
  }
}

}  // namespace modle::bigwig
