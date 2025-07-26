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
  const auto num_chroms = static_cast<std::size_t>(fp->cl->nKeys);
  Reader::Chromosomes chroms;

  const std::span<const char*> chrom_names{const_cast<const char**>(fp->cl->chrom), num_chroms};
  const std::span<const std::uint32_t> chrom_sizes{fp->cl->len, num_chroms};

  for (std::size_t i = 0; i < num_chroms; ++i) {
    chroms.emplace(chrom_names[i], utils::conditional_static_cast<bp_t>(chrom_sizes[i]));
  }
  return chroms;
}

Reader::Reader(std::filesystem::path name)
    : _fname(std::move(name)),
      _fp(bwOpen(_fname.c_str(), nullptr, "r")),
      _chroms(_fp ? read_chromosomes(_fp) : Chromosomes{}) {
  if (!_fp) {
    throw fmt::system_error(errno, "an error occurred while opening file {} for reading", _fname);
  }
}

Reader::~Reader() {
  if (_fp) {
    bwClose(_fp);
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
  validate_query(chrom, start, end);

  const Intervals intervals{
      bwGetValues(_fp, chrom.c_str(), utils::conditional_static_cast<std::uint32_t>(start),
                  utils::conditional_static_cast<std::uint32_t>(end), 1),
      &bwDestroyOverlappingIntervals};

  if (!intervals) {
    return {};
  }

  std::vector<float> values(utils::conditional_static_cast<std::size_t>(intervals->l));
  std::copy_n(intervals->value, values.size(), values.begin());
  return values;
}

auto Reader::get_intervals(const std::string& chrom, bp_t start, bp_t end) -> Intervals {
  assert(!!*this);
  validate_query(chrom, start, end);

  return Intervals{bwGetOverlappingIntervals(_fp, chrom.c_str(),
                                             utils::conditional_static_cast<std::uint32_t>(start),
                                             utils::conditional_static_cast<std::uint32_t>(end)),
                   &bwDestroyOverlappingIntervals};
}

void Reader::validate_query(const std::string& chrom, bp_t start, bp_t end) {
  auto it = _chroms.find(chrom);
  if (it == _chroms.end()) {
    throw std::runtime_error(fmt::format("query {}:{}-{}: unable to find chromosome {} in file {}",
                                         chrom, start, end, chrom, _fname));
  }

  if (start >= end) {
    throw std::logic_error(fmt::format(
        "query {}:{}-{}: start position is greater than end position", chrom, start, end));
  }

  if (end > it->second) {
    throw std::out_of_range(
        fmt::format("query {}:{}-{}: query spans past the end of chromosome", chrom, start, end));
  }
}

Writer::Writer(std::filesystem::path name, std::uint_fast8_t zoom_levels)
    : _fname(std::move(name)),
      _fp(bwOpen(_fname.c_str(), nullptr, "w")),
      _zoom_levels(zoom_levels) {
  if (!_fp) {
    throw fmt::system_error(errno, "an error occurred while opening file {} for writing", _fname);
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
  if (_fp) {
    bwClose(_fp);
  }
}

Writer& Writer::operator=(Writer&& other) noexcept {
  if (this == &other) {
    return *this;
  }

  if (_fp) {
    bwClose(_fp);
  }
  _fp = nullptr;

  _fname = std::move(other._fname);
  std::swap(_fp, other._fp);
  _zoom_levels = other._zoom_levels;
  _initialized = other._initialized;

  return *this;
}

void Writer::write_chromosomes(const std::vector<std::string>& chrom_names,
                               const std::vector<std::uint32_t>& chrom_sizes) {
  assert(chrom_names.size() == chrom_sizes.size());
  std::vector<const char*> chrom_names_ptr(chrom_names.size());
  std::transform(chrom_names.begin(), chrom_names.end(), chrom_names_ptr.begin(),
                 [](const std::string& name) { return name.data(); });

  write_chromosomes(chrom_names_ptr.data(), chrom_sizes.data(), chrom_names.size());
}

void Writer::write_chromosomes(const char* const* chrom_names, const std::uint32_t* chrom_sizes,
                               const std::size_t num_chroms) {
  assert(_fp);  // NOLINT(hicpp-no-array-decay)
  // GCC 8 and older fails to parse this if constexpr
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 9
  if (utils::ndebug_not_defined()) {
#else
  if constexpr (utils::ndebug_not_defined()) {
#endif
    if (_initialized) {
      throw std::runtime_error(
          fmt::format("bigwig::Writer::write_chromosomes() was called twice on file {}", _fname));
    }
  }

  if (bwCreateHdr(_fp, _zoom_levels)) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(
        fmt::format("failed to initialize the file header for file {}", _fname));
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
  _fp->cl =
      bwCreateChromList(const_cast<char**>(chrom_names), const_cast<std::uint32_t*>(chrom_sizes),
                        static_cast<std::int64_t>(num_chroms));
  if (!_fp->cl) {
    throw std::runtime_error(
        fmt::format("failed to create the chromosome list for file {}", _fname));
  }

  if (bwWriteHdr(_fp)) {  // NOLINT(readability-implicit-bool-conversion)
    throw std::runtime_error(fmt::format("failed to write file header to file {}", _fname));
  }
  _initialized = true;
}

}  // namespace modle::io::bigwig
