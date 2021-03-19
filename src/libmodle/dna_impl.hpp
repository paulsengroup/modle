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

#include "dna_impl.hpp"         // for FMT_COMPILE_STRING
#include "modle/bed.hpp"        // for BED
#include "modle/chr_sizes.hpp"  // for ChrSize
#include "modle/common.hpp"     // for Contacts, Bp
#include "modle/contacts.hpp"   // for ContactMatrix

namespace modle {
Chromosome::Chromosome(const chr_sizes::ChrSize &chrom, const std::vector<bed::BED> &barriers)
    : _name(chrom.name),
      _size(chrom.end),
      _start(chrom.start),
      _end(chrom.end),
      _barriers(barriers.begin(), barriers.end()) {}

Chromosome::Chromosome(chr_sizes::ChrSize &&chrom, std::vector<bed::BED> &&barriers)
    : _name(std::move(chrom.name)),
      _size(chrom.end),
      _start(chrom.start),
      _end(chrom.end),
      _barriers(std::move_iterator(barriers.begin()), std::move_iterator(barriers.end())) {}

bool Chromosome::operator==(const Chromosome &other) const {
  return this->_name == other._name && this->_size == other._size;
}

bool Chromosome::operator==(std::string_view other) const { return this->_name == other; }

bool Chromosome::operator<(const Chromosome &other) const {
  if (this->name() != other.name()) {
    return this->name() < other.name();
  }
  return this->size() < other.size();
}

bool Chromosome::Comparator::operator()(const Chromosome &c1, const Chromosome &c2) const {
  return c1 < c2;
}

bool Chromosome::Comparator::operator()(const Chromosome &c1, std::string_view c2) const {
  return c1.name() < c2;
}

bool Chromosome::Comparator::operator()(std::string_view c1, const Chromosome &c2) const {
  return c1 < c2.name();
}

bool Chromosome::Comparator::operator()(std::string_view c1, std::string_view c2) const {
  return c1 < c2;
}

void Chromosome::instantiate_contact_matrix(std::size_t bin_size, std::size_t diagonal_width) {
  assert(diagonal_width != 0);
  if (this->_contacts) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Caught an attempt to instantiate the contact matrix for '{}'"), this->_name));
  }
  const auto nrows = diagonal_width / bin_size;
  const auto ncols = (this->_end - this->_start) / bin_size;
  this->_contacts = std::make_shared<ContactMatrix<Contacts>>(nrows, ncols);
}

void Chromosome::clear_contacts() {
  assert(this->_contacts);
  this->_contacts->reset();
}

void Chromosome::add_extrusion_barrier(bed::BED barrier) { this->_barriers.emplace(barrier); }

void Chromosome::add_extrusion_barrier(bed::BED &&barrier) { this->_barriers.emplace(barrier); }

std::string_view Chromosome::name() const { return this->_name; }

Bp Chromosome::start_pos() const { return this->_start; }

Bp Chromosome::end_pos() const { return this->_end; }

Bp Chromosome::size() const { return this->_size; }
Bp Chromosome::simulated_size() const { return this->_end - this->_start; }
bool Chromosome::ok() const { return !this->_barriers.empty(); }

std::size_t Chromosome::nbarriers() const { return static_cast<size_t>(this->_barriers.size()); }

const absl::btree_set<bed::BED> &Chromosome::get_barriers() const { return this->_barriers; }

template <typename I>
void Chromosome::increment_contacts(Bp pos1, Bp pos2, Bp bin_size, I n) {
  static_assert(std::is_integral_v<I>, "n should have an integral type.");
  assert(this->_contacts);
  this->_contacts->add((pos1 - this->_start) / bin_size, (pos2 - this->_start) / bin_size, n);
}

void Chromosome::increment_contacts(Bp pos1, Bp pos2, Bp bin_size) {
  this->_contacts->increment((pos1 - this->_start) / bin_size, (pos2 - this->_start) / bin_size);
}

void Chromosome::allocate_contacts(Bp bin_size, Bp diagonal_width) {
  const auto ncols = (this->simulated_size() / bin_size) + (this->simulated_size() % bin_size != 0);
  const auto nrows =
      std::min(ncols, (diagonal_width / bin_size) + (diagonal_width % bin_size != 0));
  this->_contacts = std::make_shared<ContactMatrix<Contacts>>(nrows, ncols);
}

void Chromosome::deallocate_contacts() { this->_contacts = nullptr; }

const ContactMatrix<Contacts> &Chromosome::contacts() const {
  assert(this->_contacts);  // NOLINT
  return *this->_contacts;
}

template <typename H>
H AbslHashValue(H h, const Chromosome &c) {
  return H::combine(std::move(h), c._name);
}

}  // namespace modle
