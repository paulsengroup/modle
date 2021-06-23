#pragma once

#include <algorithm>    // for for_each
#include <cassert>      // for assert
#include <cstddef>      // for size_t
#include <type_traits>  // for is_rvalue_reference_v

#include "modle/chrom_sizes.hpp"
#include "modle/common/common.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"

namespace modle {

DISABLE_WARNING_PUSH
DISABLE_WARNING_SIGN_CONVERSION
DISABLE_WARNING_CONVERSION
template <typename I, typename>
Chromosome::Chromosome(size_t id, std::string_view chrom_name, I chrom_start, I chrom_end,
                       I chrom_size)
    : _name(chrom_name), _start(chrom_start), _end(chrom_end), _size(chrom_size), _id(id) {
  assert(chrom_start <= chrom_end);               // NOLINT
  assert(chrom_end - chrom_start <= chrom_size);  // NOLINT
}

template <typename I, typename>
Chromosome::Chromosome(size_t id, std::string_view chrom_name, I chrom_start, I chrom_end,
                       I chrom_size, const IITree<bp_t, ExtrusionBarrier>& barriers)
    : _name(chrom_name),
      _start(chrom_start),
      _end(chrom_end),
      _size(chrom_size),
      _id(id),
      _barriers(barriers) {
  assert(chrom_start <= chrom_end);               // NOLINT
  assert(chrom_end - chrom_start <= chrom_size);  // NOLINT

  _barriers.make_BST();
}

template <typename I, typename>
Chromosome::Chromosome(size_t id, std::string_view chrom_name, I chrom_start, I chrom_end,
                       I chrom_size, IITree<bp_t, ExtrusionBarrier>&& barriers)
    : _name(chrom_name),
      _start(chrom_start),
      _end(chrom_end),
      _size(chrom_size),
      _id(id),
      _barriers(std::move(barriers)) {
  assert(chrom_start <= chrom_end);               // NOLINT
  assert(chrom_end - chrom_start <= chrom_size);  // NOLINT

  _barriers.make_BST();
}
DISABLE_WARNING_POP

constexpr bp_t Chromosome::start_pos() const { return this->_start; }

constexpr bp_t Chromosome::end_pos() const { return this->_end; }

constexpr bp_t Chromosome::size() const { return this->_size; }
constexpr bp_t Chromosome::simulated_size() const { return this->_end - this->_start; }

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
