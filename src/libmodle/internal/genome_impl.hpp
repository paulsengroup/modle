// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <limits>
#include <utility>

#include "modle/common/common.hpp"

namespace modle {

namespace internal {
constexpr ContactMatrixLazy::operator bool() const noexcept {
  return this->_matrix->npixels() != 0;
}
constexpr u64 ContactMatrixLazy::nrows() const noexcept { return this->_nrows; }
constexpr u64 ContactMatrixLazy::ncols() const noexcept { return this->_ncols; }
constexpr u64 ContactMatrixLazy::npixels() const noexcept { return this->_ncols * this->_nrows; }

constexpr usize Occupancy1DLazy::size() const noexcept { return this->_size; }
}  // namespace internal

constexpr Chromosome::operator bool() const noexcept {
  return this->_id != (std::numeric_limits<usize>::max)();
}
constexpr usize Chromosome::id() const noexcept { return this->_id; }
constexpr bp_t Chromosome::size() const noexcept { return this->_size; }

constexpr bool Chromosome::operator==(const Chromosome& other) const noexcept {
  return this->_id == other._id;
}
constexpr bool Chromosome::operator!=(const Chromosome& other) const noexcept {
  return !(*this == other);
}
constexpr bool Chromosome::operator<(const Chromosome& other) const noexcept {
  return this->_id < other._id;
}
constexpr bool Chromosome::operator<=(const Chromosome& other) const noexcept {
  return this->_id <= other._id;
}
constexpr bool Chromosome::operator>(const Chromosome& other) const noexcept {
  return this->_id > other._id;
}
constexpr bool Chromosome::operator>=(const Chromosome& other) const noexcept {
  return this->_id >= other._id;
}

template <typename It>
inline GenomicInterval::GenomicInterval(usize id, const std::shared_ptr<const Chromosome>& chrom,
                                        bp_t contact_matrix_resolution, bp_t diagonal_width,
                                        It first_barrier, It last_barrier)
    : GenomicInterval(id, chrom, contact_matrix_resolution, diagonal_width, 0, chrom->size(),
                      first_barrier, last_barrier) {}

template <typename It>
inline GenomicInterval::GenomicInterval(usize id, std::shared_ptr<const Chromosome> chrom,
                                        bp_t start, bp_t end, bp_t contact_matrix_resolution,
                                        bp_t diagonal_width, It first_barrier, It last_barrier)
    : _id(id),
      _chrom(std::move(chrom)),
      _start(start),
      _end(end),
      _barriers(first_barrier, last_barrier),
      _contacts(end - start, diagonal_width, contact_matrix_resolution),
      _lef_1d_occupancy(end - start, contact_matrix_resolution) {
  assert(start <= end);
}

constexpr bool GenomicInterval::operator==(const GenomicInterval& other) const noexcept {
  return this->_id == other._id;
}
constexpr bool GenomicInterval::operator!=(const GenomicInterval& other) const noexcept {
  return !(*this == other);
}
constexpr bool GenomicInterval::operator<(const GenomicInterval& other) const noexcept {
  return this->_id < other._id;
}
constexpr bool GenomicInterval::operator<=(const GenomicInterval& other) const noexcept {
  return this->_id <= other._id;
}
constexpr bool GenomicInterval::operator>(const GenomicInterval& other) const noexcept {
  return this->_id > other._id;
}
constexpr bool GenomicInterval::operator>=(const GenomicInterval& other) const noexcept {
  return this->_id >= other._id;
}

constexpr usize GenomicInterval::id() const noexcept { return this->_id; }
constexpr bp_t GenomicInterval::start() const noexcept { return this->_start; }
constexpr bp_t GenomicInterval::end() const noexcept { return this->_end; }
constexpr bp_t GenomicInterval::size() const noexcept { return this->_end - this->_start; }
constexpr u64 GenomicInterval::npixels() const noexcept { return this->_contacts.npixels(); }

template <typename H>
H AbslHashValue(H h, const GenomicInterval& c) {
  return H::combine(std::move(h), c._id);
}

constexpr const std::vector<std::shared_ptr<const Chromosome>>& Genome::chromosomes()
    const noexcept {
  return this->_chroms;
}
constexpr usize Genome::size() const noexcept { return _size; }
constexpr usize Genome::simulated_size() const noexcept { return _simulated_size; }
constexpr usize Genome::num_barriers() const noexcept { return _num_barriers; }

}  // namespace modle

constexpr auto fmt::formatter<modle::Chromosome>::parse(fmt::format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') throw format_error("invalid format");
  // Return an iterator past the end of the parsed range:
  return ctx.end();
}

template <typename FormatContext>
inline auto fmt::formatter<modle::Chromosome>::format(const modle::Chromosome& chrom,
                                                      FormatContext& ctx) const
    -> decltype(ctx.out()) {
  if (MODLE_UNLIKELY(!chrom)) {
    return fmt::format_to(ctx.out(), FMT_STRING("null"));
  }
  return fmt::format_to(ctx.out(), FMT_STRING("{}:{}"), chrom.name(), chrom.size());
}

constexpr auto fmt::formatter<modle::GenomicInterval>::parse(fmt::format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') throw format_error("invalid format");
  // Return an iterator past the end of the parsed range:
  return ctx.end();
}

template <typename FormatContext>
inline auto fmt::formatter<modle::GenomicInterval>::format(const modle::GenomicInterval& gi,
                                                           FormatContext& ctx) const
    -> decltype(ctx.out()) {
  if (MODLE_UNLIKELY(!gi.chrom())) {
    return fmt::format_to(ctx.out(), FMT_STRING("null"));
  }
  return fmt::format_to(ctx.out(), FMT_STRING("{}:{}-{}"), gi.chrom().name(), gi.start(), gi.end());
}
