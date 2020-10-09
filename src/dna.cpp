
#include "modle/dna.hpp"

#include <algorithm>

#include "absl/strings/str_format.h"
#include "modle/extr_barrier.hpp"
#include "modle/lefs.hpp"
#include "modle/parsers.hpp"

namespace modle {

DNA::Bin::Bin(uint32_t idx, uint64_t start, uint64_t end, std::vector<ExtrusionBarrier>& barriers)
    : _idx(idx),
      _start(start),
      _end(end),
      _extr_barriers(std::make_unique<std::vector<ExtrusionBarrier>>(std::move(barriers))) {}

DNA::Bin::Bin(uint32_t idx, uint64_t start, uint64_t end)
    : _idx(idx), _start(start), _end(end), _extr_barriers(nullptr) {}

void DNA::Bin::add_barrier(ExtrusionBarrier& b) {
  if (!this->_extr_barriers) {
    this->_extr_barriers = std::make_unique<std::vector<ExtrusionBarrier>>();
  }
  this->_extr_barriers->emplace_back(b);
}

void DNA::Bin::remove_barrier(Direction d) {
  assert(d != Direction::both);
  if (!this->_extr_barriers) {
    throw std::logic_error(
        "Attempt to remove an extrusion barrier from a bin that doesn't have any!");
  }
  auto barrier = std::find_if(this->_extr_barriers->begin(), this->_extr_barriers->end(),
                              [&d](const ExtrusionBarrier& b) { return b.get_direction() == d; });
  if (barrier != this->_extr_barriers->end()) {
    this->_extr_barriers->erase(barrier);
    return;
  }
  throw std::logic_error(absl::StrFormat("Bin doesn't have barriers blocking %s!",
                                         d == Direction::fwd ? "forward" : "reverse"));
}

uint32_t DNA::Bin::get_start() const { return this->_start; }
uint32_t DNA::Bin::get_end() const { return this->_end; }
uint32_t DNA::Bin::size() const { return this->get_end() - this->get_start(); }

uint32_t DNA::Bin::add_extr_unit_binding(ExtrusionUnit* const unit) {
  if (this->_extr_units == nullptr) {
    this->_extr_units = std::make_unique<absl::flat_hash_set<ExtrusionUnit*>>();
  }
  this->_extr_units->insert(unit);
  return this->_extr_units->size();
}

uint32_t DNA::Bin::remove_extr_unit_binding(ExtrusionUnit* const unit) {
  assert(this->_extr_units != nullptr);
  if (this->_extr_units->size() == 1) {
    this->_extr_units = nullptr;
    return 0;
  }
  this->_extr_units->erase(unit);
  return this->_extr_units->size();
}

uint32_t DNA::Bin::get_n_extr_units() const {
  return this->_extr_units == nullptr ? 0 : this->_extr_units->size();
}

absl::flat_hash_set<ExtrusionUnit*>& DNA::Bin::get_extr_units() { return *this->_extr_units; }

uint32_t DNA::Bin::get_index() const { return this->_idx; }

DNA::Bin* DNA::get_ptr_to_bin_from_pos(uint32_t pos) {
  // TODO: Enable these checks only when compiling in Debug
  if (pos > this->length())
    throw std::logic_error(
        absl::StrFormat("DNA::get_ptr_to_bin_from_pos(pos=%lu): pos > this->length(): %lu > %lu\n",
                        pos, pos, this->length()));
  if ((pos / this->get_bin_size()) > this->get_n_bins())
    throw std::logic_error(
        absl::StrFormat("(pos / this->bin_size()) >= this->nbins(): (%lu / %lu) >= %lu\n", pos,
                        this->get_bin_size(), this->get_n_bins()));
  return &this->_bins[pos / this->get_bin_size()];
}

DNA::Bin& DNA::get_bin_from_pos(uint32_t pos) { return *this->get_ptr_to_bin_from_pos(pos); }

std::vector<DNA::Bin>::iterator DNA::begin() { return this->_bins.begin(); }
std::vector<DNA::Bin>::iterator DNA::end() { return this->_bins.end(); }
std::vector<DNA::Bin>::const_iterator DNA::cbegin() const { return this->_bins.cbegin(); }
std::vector<DNA::Bin>::const_iterator DNA::cend() const { return this->_bins.cend(); }

DNA::Bin* DNA::get_ptr_to_prev_bin(const DNA::Bin& current_bin) {
  assert(current_bin._idx > 0);
  return &this->_bins[current_bin._idx - 1];
}

DNA::Bin& DNA::get_prev_bin(const Bin& current_bin) {
  return *this->get_ptr_to_prev_bin(current_bin);
}

DNA::Bin& DNA::get_next_bin(const DNA::Bin& current_bin) {
  return *this->get_ptr_to_next_bin(current_bin);
}

DNA::Bin* DNA::get_ptr_to_next_bin(const DNA::Bin& current_bin) {
  assert(current_bin._idx < this->get_n_bins());
  return &this->_bins[current_bin._idx + 1];
}

DNA::Bin& DNA::get_first_bin() { return *this->get_ptr_to_first_bin(); }

DNA::Bin* DNA::get_ptr_to_first_bin() {
  assert(!this->_bins.empty());
  return &this->_bins.front();
}

DNA::Bin& DNA::get_last_bin() { return *this->get_ptr_to_last_bin(); }

DNA::Bin* DNA::get_ptr_to_last_bin() {
  assert(!this->_bins.empty());
  return &this->_bins.back();
}

void DNA::Bin::remove_barriers() { this->_extr_barriers = nullptr; }

ExtrusionBarrier* DNA::Bin::get_next_extr_barrier(ExtrusionBarrier* b, Direction d) const {
  if (!this->_extr_barriers) return nullptr;
  auto b_itr = b ? this->_extr_barriers->begin() + std::distance(this->_extr_barriers->data(), b)
                 : this->_extr_barriers->begin();  // ptr to iter

  auto barrier =
      std::find_if(b_itr, this->_extr_barriers->end(), [&d](const ExtrusionBarrier& barr) {
        return d == DNA::Direction::both || barr.get_direction() == d;
      });
  if (barrier != this->_extr_barriers->end()) return barrier.operator->();
  return nullptr;
}

std::vector<ExtrusionBarrier>* DNA::Bin::get_all_extr_barriers() const {
  return this->_extr_barriers.get();
}

void DNA::Bin::add_extr_barrier(ExtrusionBarrier b) {
  if (!this->_extr_barriers) {
    this->_extr_barriers =
        std::make_unique<std::vector<ExtrusionBarrier>>(std::vector<ExtrusionBarrier>{b});
  } else {
    this->_extr_barriers->emplace_back(b);
  }
}

void DNA::Bin::add_extr_barrier(double prob_of_barrier_block, DNA::Direction direction) {
  if (!this->_extr_barriers) {
    this->_extr_barriers = std::make_unique<std::vector<ExtrusionBarrier>>(
        std::vector<ExtrusionBarrier>{{prob_of_barrier_block, direction}});
  } else {
    this->_extr_barriers->emplace_back(prob_of_barrier_block, direction);
  }
}

DNA::DNA(uint64_t length, uint32_t bin_size)
    : _bins(make_bins(length, bin_size)), _length(length) {}

void DNA::add_barrier(ExtrusionBarrier& b, uint32_t pos) {
  assert(b.get_prob_of_block() > 0);
  assert(pos <= this->length());
  assert(pos / this->get_bin_size() < this->get_n_bins());
  this->_bins[pos / this->get_bin_size()].add_barrier(b);
}

void DNA::remove_barrier(uint32_t pos, Direction direction) {
  assert(pos <= this->length());
  assert(pos / this->get_bin_size() < this->get_n_bins());
  this->_bins[pos / this->get_bin_size()].remove_barrier(direction);
}

std::vector<DNA::Bin> DNA::make_bins(uint64_t length, uint32_t bin_size) {
  std::vector<DNA::Bin> bins;
  if (length <= bin_size) {
    bins.emplace_back(DNA::Bin{0, 0, length});
    return bins;
  }
  bins.reserve((length / bin_size) + 1);
  uint64_t start = 0;
  for (auto end = bin_size; end <= length; end += bin_size) {
    bins.emplace_back(DNA::Bin{static_cast<uint32_t>(bins.size()), start, end});
    start = end + 1;
  }
  // TODO: Make sure adding a smaller bin at the end when (length % bin_size) != 0 is the right
  // thing to do
  if (bins.back()._end < length) {
    bins.emplace_back(DNA::Bin{static_cast<uint32_t>(bins.size()), start, length});
  }
  return bins;
}

uint32_t DNA::length() const { return this->_length; }

uint32_t DNA::get_n_bins() const { return this->_bins.size(); }

uint32_t DNA::get_bin_size() const { return this->_bins[0]._end - this->_bins[0]._start; }

uint32_t DNA::get_n_barriers() const {
  return std::accumulate(this->_bins.begin(), this->_bins.end(), 0UL,
                         [](uint32_t accumulator, const DNA::Bin& b) {
                           if (b._extr_barriers) accumulator += b._extr_barriers->size();
                           return accumulator;
                         });
}

Chromosome::Chromosome(std::string name, DNA dna)
    : name(std::move(name)),
      dna(std::move(dna)),
      // TODO: Make this a tunable, the first parameter controls the width of the diagonal that we
      // are actually storing
      contacts(500'000 / (this->length() / this->n_bins()), this->dna.get_n_bins()) {}

uint32_t Chromosome::length() const { return this->dna.length(); }
uint32_t Chromosome::n_bins() const { return this->dna.get_n_bins(); }
uint32_t Chromosome::n_barriers() const { return this->barriers.size(); }
void Chromosome::write_contacts_to_tsv(const std::string& path_to_file, bool complete) const {
  if (complete)
    this->contacts.write_full_matrix_to_tsv(path_to_file);
  else
    this->contacts.write_to_tsv(path_to_file);
}

}  // namespace modle