#include "modle/dna.hpp"

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

#include "modle/bed.hpp"          // for BED
#include "modle/chrom_sizes.hpp"  // for ChromSize
#include "modle/common.hpp"       // for bp_t, contacts_t
#include "modle/contacts.hpp"     // for ContactMatrix
#include "modle/utils.hpp"        // for ndebug_defined

namespace modle {
Chromosome::Chromosome(const chrom_sizes::ChromSize &chrom, size_t id,
                       const std::vector<bed::BED> &barriers)
    : _id(id),
      _name(chrom.name),
      _size(chrom.end),
      _start(chrom.start),
      _end(chrom.end),
      _barriers(barriers.begin(), barriers.end()) {}

Chromosome::Chromosome(chrom_sizes::ChromSize &&chrom, size_t id, std::vector<bed::BED> &&barriers)
    : _id(id),
      _name(std::move(chrom.name)),
      _size(chrom.end),
      _start(chrom.start),
      _end(chrom.end),
      _barriers(std::move_iterator(barriers.begin()), std::move_iterator(barriers.end())) {}

bool Chromosome::operator==(const Chromosome &other) const noexcept(utils::ndebug_defined()) {
  if (this->_id == std::numeric_limits<size_t>::max()) {
    return this->name() == other.name() && this->size() == other.size();
  }
  return this->_id == other._id;
}

bool Chromosome::operator==(std::string_view other) const noexcept(utils::ndebug_defined()) {
  return this->name() == other;
}

bool Chromosome::operator<(const Chromosome &other) const noexcept(utils::ndebug_defined()) {
  if (this->_id != std::numeric_limits<size_t>::max() &&
      other._id != std::numeric_limits<size_t>::max()) {
    return this->_id < other._id;
  }

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

void Chromosome::add_extrusion_barrier(bed::BED barrier) { this->_barriers.emplace(barrier); }

void Chromosome::add_extrusion_barrier(bed::BED &&barrier) { this->_barriers.emplace(barrier); }

size_t Chromosome::id() const { return this->_id; }
std::string_view Chromosome::name() const { return this->_name; }
const char *Chromosome::name_cstr() const { return this->_name.c_str(); }

bool Chromosome::ok() const { return !this->_barriers.empty(); }

size_t Chromosome::nlefs(double nlefs_per_mbp) const {  // NOLINTNEXTLINE
  return static_cast<size_t>((static_cast<double>(this->simulated_size()) / 1.0e6) * nlefs_per_mbp);
}
size_t Chromosome::nbarriers() const { return static_cast<size_t>(this->_barriers.size()); }

size_t Chromosome::num_valid_barriers() const {
  return std::count_if(this->get_barriers().begin(), this->get_barriers().end(),
                       [](const auto &b) { return b.strand == '-' || b.strand == '+'; });
}

const absl::btree_set<bed::BED> &Chromosome::get_barriers() const { return this->_barriers; }

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

uint64_t Chromosome::hash(uint64_t seed, size_t cell_id) {
  auto handle_errors = [&](const auto &status) {
    if (status == XXH_ERROR || !this->_xxh_state) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Failed to hash '{}' for cell #{} using seed {}"), this->name(),
                      cell_id, seed));
    }
  };

  handle_errors(XXH3_64bits_reset_withSeed(this->_xxh_state.get(), seed));
  handle_errors(XXH3_64bits_update(this->_xxh_state.get(), this->_name.data(),
                                   this->_name.size() * sizeof(char)));
  handle_errors(
      XXH3_64bits_update(this->_xxh_state.get(), &this->_size, sizeof(decltype(this->_size))));
  handle_errors(XXH3_64bits_update(this->_xxh_state.get(), &cell_id, sizeof(decltype(cell_id))));

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  return static_cast<uint64_t>(XXH3_64bits_digest(this->_xxh_state.get()));
  DISABLE_WARNING_POP
}

template <typename H>
H AbslHashValue(H h, const Chromosome &c) {
  return H::combine(std::move(h), c._name);
}

}  // namespace modle
