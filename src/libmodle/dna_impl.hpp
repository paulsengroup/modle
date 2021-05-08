#pragma once

// IWYU pragma: private, include "modle/dna.hpp"

#include <absl/container/btree_set.h>  // for btree_set
#include <fmt/format.h>                // for format, FMT_STRING

#include <algorithm>    // for min
#include <cassert>      // for assert
#include <cstddef>      // for size_t
#include <cstdint>      // for uint32_t
#include <iterator>     // for move_iterator
#include <memory>       // for shared_ptr, make_shared, __shared_ptr_access
#include <stdexcept>    // for runtime_error
#include <string>       // for char_traits, string, operator==, basic_string<>::_...
#include <string_view>  // for operator<, string_view, operator!=, operator==
#include <utility>      // for move
#include <vector>       // for vector

#include "dna_impl.hpp"           // for FMT_COMPILE_STRING
#include "modle/bed.hpp"          // for BED
#include "modle/chrom_sizes.hpp"  // for ChromSize
#include "modle/common.hpp"       // for bp_t, contacts_t
#include "modle/contacts.hpp"     // for ContactMatrix
#include "modle/utils.hpp"        // for ndebug_defined

namespace modle {
Chromosome::Chromosome(const chrom_sizes::ChromSize &chrom, const std::vector<bed::BED> &barriers)
    : _name(chrom.name),
      _size(chrom.end),
      _start(chrom.start),
      _end(chrom.end),
      _barriers(barriers.begin(), barriers.end()) {}

Chromosome::Chromosome(chrom_sizes::ChromSize &&chrom, std::vector<bed::BED> &&barriers)
    : _name(std::move(chrom.name)),
      _size(chrom.end),
      _start(chrom.start),
      _end(chrom.end),
      _barriers(std::move_iterator(barriers.begin()), std::move_iterator(barriers.end())) {}

bool Chromosome::operator==(const Chromosome &other) const noexcept(utils::ndebug_defined()) {
  return this->name() == other.name() && this->size() == other.size();
}

bool Chromosome::operator==(std::string_view other) const noexcept(utils::ndebug_defined()) {
  return this->name() == other;
}

bool Chromosome::operator<(const Chromosome &other) const noexcept(utils::ndebug_defined()) {
  if (this->name() != other.name()) {
    return this->name() < other.name();
  }
  return this->size() < other.size();
}

bool Chromosome::Comparator::operator()(const Chromosome &c1, const Chromosome &c2) const
    noexcept(utils::ndebug_defined()) {
  return c1 < c2;
}

bool Chromosome::Comparator::operator()(const Chromosome &c1, std::string_view c2) const
    noexcept(utils::ndebug_defined()) {
  return c1.name() < c2;
}

bool Chromosome::Comparator::operator()(std::string_view c1, const Chromosome &c2) const
    noexcept(utils::ndebug_defined()) {
  return c1 < c2.name();
}

bool Chromosome::Comparator::operator()(std::string_view c1, std::string_view c2) const
    noexcept(utils::ndebug_defined()) {
  return c1 < c2;
}

void Chromosome::instantiate_contact_matrix(size_t bin_size, size_t diagonal_width) {
  assert(diagonal_width != 0);  // NOLINT
#ifndef NDEBUG
  if (this->_contacts) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Caught an attempt to re-instantiate the contact matrix for '{}'"),
                    this->_name));
  }
#endif
  const auto nrows = diagonal_width / bin_size;
  const auto ncols = (this->_end - this->_start) / bin_size;
  this->_contacts = std::make_shared<ContactMatrix<contacts_t>>(nrows, ncols);
}

void Chromosome::clear_contacts() {
  assert(this->_contacts);
  this->_contacts->reset();
}

void Chromosome::add_extrusion_barrier(bed::BED barrier) { this->_barriers.emplace(barrier); }

void Chromosome::add_extrusion_barrier(bed::BED &&barrier) { this->_barriers.emplace(barrier); }

std::string_view Chromosome::name() const { return this->_name; }
const char *Chromosome::name_cstr() const { return this->_name.c_str(); }

constexpr bp_t Chromosome::start_pos() const { return this->_start; }

constexpr bp_t Chromosome::end_pos() const { return this->_end; }

constexpr bp_t Chromosome::size() const { return this->_size; }
constexpr bp_t Chromosome::simulated_size() const { return this->_end - this->_start; }
bool Chromosome::ok() const { return !this->_barriers.empty(); }

size_t Chromosome::nbarriers() const { return static_cast<size_t>(this->_barriers.size()); }

const absl::btree_set<bed::BED> &Chromosome::get_barriers() const { return this->_barriers; }

template <typename I>
void Chromosome::increment_contacts(bp_t pos1, bp_t pos2, bp_t bin_size, I n) {
  static_assert(std::is_integral_v<I>, "n should have an integral type.");
  assert(this->_contacts);
  this->_contacts->add((pos1 - this->_start) / bin_size, (pos2 - this->_start) / bin_size, n);
}

void Chromosome::increment_contacts(bp_t pos1, bp_t pos2, bp_t bin_size) {
  this->_contacts->increment((pos1 - this->_start) / bin_size, (pos2 - this->_start) / bin_size);
}

void Chromosome::allocate_contacts(bp_t bin_size, bp_t diagonal_width) {
  const auto ncols = (this->simulated_size() / bin_size) + (this->simulated_size() % bin_size != 0);
  const auto nrows =
      std::min(ncols, (diagonal_width / bin_size) + (diagonal_width % bin_size != 0));
  this->_contacts = std::make_shared<ContactMatrix<contacts_t>>(nrows, ncols);
}

void Chromosome::deallocate_contacts() { this->_contacts = nullptr; }

const ContactMatrix<contacts_t> &Chromosome::contacts() const {
  assert(this->_contacts);  // NOLINT
  return *this->_contacts;
}

template <typename H>
H AbslHashValue(H h, const Chromosome &c) {
  return H::combine(std::move(h), c._name);
}

}  // namespace modle
