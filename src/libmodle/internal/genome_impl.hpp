// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>    // for for_each
#include <cassert>      // for assert
#include <cstddef>      // for size_t
#include <type_traits>  // for is_rvalue_reference_v

#include "modle/chrom_sizes.hpp"
#include "modle/common/common.hpp"  // for bp_t
#include "modle/common/suppress_compiler_warnings.hpp"

namespace modle {

constexpr bp_t Chromosome::start_pos() const { return this->_start; }

constexpr bp_t Chromosome::end_pos() const { return this->_end; }

constexpr bp_t Chromosome::size() const { return this->_size; }
constexpr bp_t Chromosome::simulated_size() const { return this->_end - this->_start; }

constexpr size_t Chromosome::npixels(bp_t diagonal_width, bp_t bin_size) const {
  const auto npix1 = (diagonal_width + bin_size - 1) / bin_size;
  const auto npix2 = (this->simulated_size() + bin_size - 1) / bin_size;

  return npix1 * npix2;
}

template <typename H>
H AbslHashValue(H h, const Chromosome& c) {
  return H::combine(std::move(h), c._name);
}

template <typename Iter, typename>
void Chromosome::add_extrusion_barrier(Iter barriers_begin, Iter barriers_end) {
  std::for_each(barriers_begin, barriers_end, [this](const auto& barrier) {
    this->_barriers.insert(barrier.chrom_start, barrier.chrom_end, barrier);
  });
}
}  // namespace modle
