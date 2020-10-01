
#include "modle/dna.hpp"

#include <algorithm>

#include "absl/strings/str_format.h"
#include "modle/extr_barrier.hpp"
#include "modle/lefs.hpp"
#include "modle/parsers.hpp"

namespace modle {

DNA::DNA(uint64_t length, uint32_t bin_size)
    : _bins(make_bins(length, bin_size)), _length(length) {}

void DNA::add_fwd_barrier(const ExtrusionBarrier& barrier, uint32_t pos) {
  assert(pos <= this->length());
  assert(pos / this->bin_size() < this->n_bins());
  this->_bins[pos / this->bin_size()]->add_fwd_barrier(&barrier);
}

void DNA::add_rev_barrier(const ExtrusionBarrier& barrier, uint32_t pos) {
  assert(pos <= this->length());
  assert(pos / this->bin_size() < this->n_bins());
  this->_bins[pos / this->bin_size()]->add_rev_barrier(&barrier);
}

void DNA::remove_fwd_barrier(uint32_t pos) {
  assert(pos <= this->length());
  assert(pos / this->bin_size() < this->n_bins());
  this->_bins[pos / this->bin_size()]->remove_fwd_barrier();
}

void DNA::remove_rev_barrier(uint32_t pos) {
  assert(pos <= this->length());
  assert(pos / this->bin_size() < this->n_bins());
  this->_bins[pos / this->bin_size()]->remove_rev_barrier();
}

std::vector<std::shared_ptr<DNA::Bin>> DNA::make_bins(uint64_t length, uint32_t bin_size) {
  std::vector<std::shared_ptr<DNA::Bin>> bins;
  if (length <= bin_size) {
    bins.emplace_back(std::make_shared<DNA::Bin>(DNA::Bin{0, 0, length}));
    return bins;
  }
  bins.reserve(length / bin_size);
  uint32_t start = 0;
  for (uint32_t end = bin_size; end <= length; end += bin_size) {
    bins.emplace_back(
        std::make_shared<DNA::Bin>(DNA::Bin{static_cast<uint32_t>(bins.size()), start, end}));
    start = end + 1;
  }
  if (bins.back()->_end < length) {
    bins.emplace_back(
        std::make_shared<DNA::Bin>(DNA::Bin{static_cast<uint32_t>(bins.size()), start, length}));
  }
  return bins;
}

uint32_t DNA::length() const { return this->_length; }

uint32_t DNA::n_bins() const { return this->_bins.size(); }

uint32_t DNA::bin_size() const { return this->_bins[0]->_end - this->_bins[0]->_start; }

std::shared_ptr<DNA::Bin> DNA::get_ptr_to_bin_from_pos(uint32_t pos) {
  if (pos > this->length())
    throw std::logic_error(
        absl::StrFormat("DNA::get_ptr_to_bin_from_pos(pos=%lu): pos > this->length(): %lu > %lu\n",
                        pos, pos, this->length()));
  if ((pos / this->bin_size()) >= this->n_bins())
    throw std::logic_error(
        absl::StrFormat("(pos / this->bin_size()) >= this->nbins(): (%lu / %lu) >= %lu\n", pos,
                        this->bin_size(), this->n_bins()));
  return this->_bins[pos / this->bin_size()];
}

std::vector<std::shared_ptr<DNA::Bin>>::iterator DNA::begin() { return this->_bins.begin(); }
std::vector<std::shared_ptr<DNA::Bin>>::iterator DNA::end() { return this->_bins.end(); }
std::vector<std::shared_ptr<DNA::Bin>>::const_iterator DNA::cbegin() const {
  return this->_bins.cbegin();
}
std::vector<std::shared_ptr<DNA::Bin>>::const_iterator DNA::cend() const {
  return this->_bins.cend();
}

std::shared_ptr<DNA::Bin> DNA::get_ptr_to_previous_bin(const DNA::Bin* const current_bin) {
  assert(current_bin);
  assert(current_bin->_idx > 0);
  return this->_bins.at(current_bin->_idx - 1);
}

std::shared_ptr<DNA::Bin> DNA::get_ptr_to_next_bin(const DNA::Bin* const current_bin) {
  assert(current_bin);
  assert(current_bin->_idx < this->n_bins());
  return this->_bins.at(current_bin->_idx + 1);
}

DNA::Bin::Bin(uint32_t idx, uint64_t start, uint64_t end, ExtrusionBarrier const* fwd_barrier,
              ExtrusionBarrier const* rev_barrier)
    : _idx(idx), _start(start), _end(end), _fwd_barrier(fwd_barrier), _rev_barrier(rev_barrier) {}

DNA::Bin::Bin(uint32_t idx, uint64_t start, uint64_t end)
    : _idx(idx), _start(start), _end(end), _fwd_barrier(nullptr), _rev_barrier(nullptr) {}

bool DNA::Bin::has_fwd_barrier() const { return this->_fwd_barrier != nullptr; }
bool DNA::Bin::has_rev_barrier() const { return this->_rev_barrier != nullptr; }
void DNA::Bin::add_fwd_barrier(const ExtrusionBarrier* const barrier) {
  this->_fwd_barrier = barrier;
}
void DNA::Bin::add_rev_barrier(const ExtrusionBarrier* const barrier) {
  this->_rev_barrier = barrier;
}

void DNA::Bin::remove_fwd_barrier() { this->_fwd_barrier = nullptr; }
void DNA::Bin::remove_rev_barrier() { this->_rev_barrier = nullptr; }
void DNA::Bin::remove_barriers() {
  this->remove_fwd_barrier();
  this->remove_rev_barrier();
}

uint32_t DNA::Bin::get_start() const { return this->_start; }
uint32_t DNA::Bin::get_end() const { return this->_end; }
uint32_t DNA::Bin::get_center() const { return (this->get_start() + this->get_end()) / 2; }
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

bool DNA::Bin::pos_is_occupied(uint32_t pos) const {
  return std::any_of(*this->_extr_units->begin(), *this->_extr_units->end(),
                     [&](const auto& extr) { return extr.get_pos() == pos; });
}

uint32_t DNA::Bin::n_extruders() const {
  return this->_extr_units == nullptr ? 0 : this->_extr_units->size();
}

absl::flat_hash_set<ExtrusionUnit*>& DNA::Bin::get_extruders() { return *this->_extr_units; }

uint32_t DNA::Bin::get_index() const { return this->_idx; }

}  // namespace modle