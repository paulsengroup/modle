// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/bigwig/bigwig.hpp"

#include <fmt/format.h>

#include <cassert>
#include <cerrno>
#include <filesystem>
#include <iosfwd>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>

#include "libbigwig/bigWig.h"
#include "modle/common/common.hpp"
#include "modle/common/fmt_helpers.hpp"

#ifndef NOCURL
#error "NOCURL is not defined. Please file an issue on GitHub."
#endif

namespace modle::io::bigwig {

Reader::Reader(Reader&& other) noexcept
    : _fname(std::move(other._fname)), _fp(other._fp), _chroms(std::move(other._chroms)) {
  other._fp = nullptr;
}

[[nodiscard]] static Reader::Chromosomes read_chromosomes(bigWigFile_t* fp) {
  assert(fp);
  const auto num_chroms = static_cast<usize>(fp->cl->nKeys);
  Reader::Chromosomes chroms;

  const auto chrom_names = absl::MakeConstSpan(fp->cl->chrom, num_chroms);
  const auto chrom_sizes = absl::MakeConstSpan(fp->cl->len, num_chroms);

  for (usize i = 0; i < num_chroms; ++i) {
    chroms.emplace(chrom_names[i], utils::conditional_static_cast<bp_t>(chrom_sizes[i]));
  }
  return chroms;
}

Reader::Reader(std::filesystem::path name)
    : _fname(std::move(name)),
      _fp(bwOpen(_fname.c_str(), nullptr, "r")),
      _chroms(_fp ? read_chromosomes(_fp) : Chromosomes{}) {
  if (!_fp) {
    throw fmt::system_error(
        errno, FMT_STRING("an error occurred while opening file {} for reading"), this->_fname);
  }
}

Reader::~Reader() {
  if (this->_fp) {
    bwClose(this->_fp);
  }
}

Reader& Reader::operator=(Reader&& other) noexcept {
  if (this == &other) {
    return *this;
  }

  if (_fp) {
    bwClose(_fp);
  }
  _fp = nullptr;

  _fname = std::move(other._fname);
  std::swap(_fp, other._fp);
  other._fp = nullptr;
  _chroms = std::move(other._chroms);

  return *this;
}

std::vector<float> Reader::read_values(const std::string& chrom, bp_t start, bp_t end) {
  assert(!!*this);
  this->validate_query(chrom, start, end);

  const Intervals intervals{
      bwGetValues(this->_fp, chrom.c_str(), utils::conditional_static_cast<u32>(start),
                  utils::conditional_static_cast<u32>(end), 1),
      &bwDestroyOverlappingIntervals};

  if (!intervals) {
    return {};
  }

  std::vector<float> values(utils::conditional_static_cast<usize>(intervals->l));
  std::copy_n(intervals->value, values.size(), values.begin());
  return values;
}

auto Reader::get_intervals(const std::string& chrom, bp_t start, bp_t end) -> Intervals {
  assert(!!*this);
  this->validate_query(chrom, start, end);

  return Intervals{bwGetOverlappingIntervals(this->_fp, chrom.c_str(),
                                             utils::conditional_static_cast<u32>(start),
                                             utils::conditional_static_cast<u32>(end)),
                   &bwDestroyOverlappingIntervals};
}

void Reader::validate_query(const std::string& chrom, bp_t start, bp_t end) {
  auto it = this->_chroms.find(chrom);
  if (it == this->_chroms.end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("query {}:{}-{}: unable to find chromosome {} in file {}"), chrom,
                    start, end, chrom, this->_fname));
  }

  if (start >= end) {
    throw std::logic_error(
        fmt::format(FMT_STRING("query {}:{}-{}: start position is greater than end position"),
                    chrom, start, end));
  }

  if (end > it->second) {
    throw std::out_of_range(fmt::format(
        FMT_STRING("query {}:{}-{}: query spans past the end of chromosome"), chrom, start, end));
  }
}

Writer::Writer(std::filesystem::path name, std::uint_fast8_t zoom_levels)
    : _fname(std::move(name)),
      _fp(bwOpen(_fname.c_str(), nullptr, "w")),
      _zoom_levels(zoom_levels) {
  if (!_fp) {
    throw fmt::system_error(
        errno, FMT_STRING("an error occurred while opening file {} for writing"), this->_fname);
  }
}

Writer::Writer(Writer&& other) noexcept
    : _fname(std::move(other._fname)),
      _fp(other._fp),
      _zoom_levels(other._zoom_levels),
      _initialized(other._initialized) {
  other._fp = nullptr;
}

Writer::~Writer() {
  if (this->_fp) {
    bwClose(this->_fp);
  }
}

Writer& Writer::operator=(Writer&& other) noexcept {
  if (this == &other) {
    return *this;
  }

  if (this->_fp) {
    bwClose(this->_fp);
  }
  this->_fp = nullptr;

  this->_fname = std::move(other._fname);
  std::swap(this->_fp, other._fp);
  this->_zoom_levels = other._zoom_levels;
  this->_initialized = other._initialized;

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
        fmt::format(FMT_STRING("failed to initialize the file header for file {}"), this->_fname));
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
  this->_fp->cl = bwCreateChromList(const_cast<char**>(chrom_names), const_cast<u32*>(chrom_sizes),
                                    static_cast<i64>(num_chroms));
  if (!this->_fp->cl) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to create the chromosome list for file {}"), this->_fname));
  }

  if (bwWriteHdr(this->_fp)) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to write file header to file {}"), this->_fname));
  }
  this->_initialized = true;
}

}  // namespace modle::io::bigwig
