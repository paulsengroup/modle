// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/bigwig/bigwig.hpp"

#include <fmt/format.h>  // for format, FMT_STRING

#include <cassert>     // for assert
#include <cerrno>      // for errno
#include <filesystem>  // for file_size
#include <iosfwd>      // for streamsize
#include <mutex>       // for scoped_lock
#include <stdexcept>   // for runtime_error
#include <string>      // for string
#include <utility>     // for move

#include "libbigwig/bigWig.h"       // for bwCleanup, bwClose, bwAddIntervalSpanSteps, bwCreateChromList
#include "modle/common/common.hpp"  // for u32, u64, i32, i64
#include "modle/common/fmt_helpers.hpp"

namespace modle::io::bigwig {

Writer::Writer(std::filesystem::path name, std::uint_fast8_t zoom_levels, usize buff_size)
    : _fname(std::move(name)), _zoom_levels(zoom_levels), _buff_size(buff_size) {
  std::scoped_lock<std::mutex> l(Writer::_global_state_mutex);
  if (Writer::_global_bigwig_files_opened == 0) {
    if (bwInit(this->_buff_size)) {  // NOLINT(readability-implicit-bool-conversion)
      throw std::runtime_error(fmt::format(
          FMT_STRING("Failed to initialize libBigWig global state while opening file: {}"),
          this->_fname));
    }
    ++Writer::_global_bigwig_files_opened;
  }

  auto tmp_str = this->_fname.string();
  this->_fp = bwOpen(tmp_str.data(), nullptr, "w");
  if (!this->_fp) {
    if (--Writer::_global_bigwig_files_opened == 0) {
      bwCleanup();
    }
    throw fmt::system_error(
        errno, FMT_STRING("An error occurred while opening file {} for writing"), this->_fname);
  }
}

Writer::Writer(Writer&& other) noexcept
    : _fname(std::move(other._fname)),
      _fp(other._fp),
      _zoom_levels(other._zoom_levels),
      _buff_size(other._buff_size),
      _initialized(other._initialized) {
  other._fp = nullptr;
}

Writer::~Writer() {
  if (this->_fp) {
    bwClose(this->_fp);
    std::scoped_lock<std::mutex> l(Writer::_global_state_mutex);
    if (--Writer::_global_bigwig_files_opened == 0) {
      bwCleanup();
    }
  }
}

Writer& Writer::operator=(Writer&& other) noexcept {
  if (this == &other) {
    return *this;
  }

  this->_fname = std::move(other._fname);
  this->_fp = other._fp;
  this->_zoom_levels = other._zoom_levels;
  this->_buff_size = other._buff_size;
  this->_initialized = other._initialized;
  other._fp = nullptr;

  return *this;
}

void Writer::write_chromosomes(const std::vector<std::string>& chrom_names,
                               const std::vector<u32>& chrom_sizes) {
  assert(chrom_names.size() == chrom_sizes.size());
  std::vector<const char*> chrom_names_ptr(chrom_names.size());
  std::transform(chrom_names.begin(), chrom_names.end(), chrom_names_ptr.begin(),
                 [](const std::string& name) { return name.data(); });

  this->write_chromosomes(chrom_names_ptr.data(), chrom_sizes.data(), chrom_names.size());
}

void Writer::write_chromosomes(const char* const* chrom_names, const u32* chrom_sizes,
                               const usize num_chroms) {
  assert(this->_fp);  // NOLINT(hicpp-no-array-decay)
  // GCC 8 and older fails to parse this if constexpr
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 9
  if (utils::ndebug_not_defined()) {
#else
  if constexpr (utils::ndebug_not_defined()) {
#endif
    if (this->_initialized) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("bigwig::Writer::write_chromosomes() was called twice on file {}"),
                      this->_fname));
    }
  }

  if (bwCreateHdr(this->_fp, this->_zoom_levels)) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to initialize the file header for file {}"), this->_fname));
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
  this->_fp->cl = bwCreateChromList(const_cast<char**>(chrom_names), const_cast<u32*>(chrom_sizes),
                                    static_cast<i64>(num_chroms));
  if (!this->_fp->cl) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to create the chromosome list for file {}"), this->_fname));
  }

  if (bwWriteHdr(this->_fp)) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to write file header to file {}"), this->_fname));
  }
  this->_initialized = true;
}

const std::filesystem::path& Writer::path() const noexcept { return this->_fname; }

}  // namespace modle::io::bigwig
