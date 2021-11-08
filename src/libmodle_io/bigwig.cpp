// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/bigwig/bigwig.hpp"

#include <fmt/format.h>  // for format, FMT_STRING
#include <fmt/ostream.h>

#include <boost/filesystem/path.hpp>  // for file_size
#include <cassert>                    // for assert
#include <cerrno>                     // for errno
#include <iosfwd>                     // for streamsize
#include <mutex>                      // for scoped_lock
#include <stdexcept>                  // for runtime_error
#include <string>                     // for string
#include <utility>                    // for move

#include "libBigWig/bigWig.h"  // for bwCleanup, bwClose, bwAddIntervalSpanSteps, bwCreateChromList
#include "modle/common/common.hpp"  // for u32, u64, i32, i64
#include "modle/common/utils.hpp"   // for ndebug_not_defined

namespace modle::io::bigwig {

Writer::Writer(boost::filesystem::path name, uint_fast8_t zoom_levels, usize buff_size)
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
      _initialized(other._initialized),
      _offset(other._offset) {
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
  this->_offset = other._offset;
  other._fp = nullptr;

  return *this;
}

void Writer::write_chromosomes(const char* const* chrom_names, const u32* chrom_sizes,
                               const usize num_chroms) {
  assert(this->_fp);  // NOLINT
  // GCC7 fails to parse this if constexpr
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ == 7
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

  this->_fp->cl = bwCreateChromList(chrom_names, chrom_sizes, static_cast<i64>(num_chroms));
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

const boost::filesystem::path& Writer::path() const noexcept { return this->_fname; }

}  // namespace modle::io::bigwig
