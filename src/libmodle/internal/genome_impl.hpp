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
template <typename Iter, typename I, typename, typename>
Chromosome::Chromosome(size_t id, std::string_view chrom_name, I chrom_start, I chrom_end,
                       I chrom_size, Iter barriers_begin, Iter barriers_end)
    : _name(chrom_name), _start(chrom_start), _end(chrom_end), _size(chrom_size), _id(id) {
  assert(chrom_start <= chrom_end);               // NOLINT
  assert(chrom_end - chrom_start <= chrom_size);  // NOLINT
  if (const auto size = std::distance(barriers_begin, barriers_end); size > 0) {
    _barriers.reserve(size);

    using Iter_t = typename std::iterator_traits<Iter>::value_type;
    if constexpr (std::is_rvalue_reference_v<Iter_t>) {
      std::for_each(barriers_begin, barriers_end, [&](auto&& barrier) {
        _barriers.emplace(barrier.chrom_start, barrier.chrom_end, std::move(barrier));
      });
    } else {
      std::for_each(barriers_begin, barriers_end, [&](const auto& barrier) {
        _barriers.insert(barrier.chrom_start, barrier.chrom_end, barrier);
      });
    }
  }
  _barriers.make_BST();
}

template <typename I, typename>
Chromosome::Chromosome(size_t id, std::string_view chrom_name, I chrom_start, I chrom_end,
                       I chrom_size, const interval_tree_value_t& barriers)
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
                       I chrom_size, interval_tree_value_t&& barriers)
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

template <typename Iter, typename>
Chromosome::Chromosome(size_t id, const bed::BED& chrom, Iter barriers_begin, Iter barriers_end)
    : Chromosome(id, chrom.chrom, chrom.thick_start, chrom.thick_end, chrom.size(), barriers_begin,
                 barriers_end){};

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
