#pragma once

#include <absl/container/btree_set.h>
#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "modle/bed.hpp"
#include "modle/chr_sizes.hpp"

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

void Chromosome::add_extrusion_barrier(bed::BED barrier) {
  this->_barriers.emplace(std::move(barrier));
}

void Chromosome::add_extrusion_barrier(bed::BED &&barrier) {
  this->_barriers.emplace(std::move(barrier));
}

std::string_view Chromosome::name() const { return this->_name; }

std::size_t Chromosome::start_pos() const { return this->_start; }

std::size_t Chromosome::end_pos() const { return this->_end; }

std::size_t Chromosome::size() const { return this->_size; }
std::size_t Chromosome::simulated_size() const { return this->_end - this->_start; }
bool Chromosome::ok() const { return !this->_barriers.empty(); }

std::size_t Chromosome::nbarriers() const { return static_cast<size_t>(this->_barriers.size()); }

const absl::btree_set<bed::BED> &Chromosome::get_barriers() const { return this->_barriers; }

template <typename H>
H AbslHashValue(H h, const Chromosome &c) {
  return H::combine(std::move(h), c._name);
}

}  // namespace modle
