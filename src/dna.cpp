#include "modle/dna.hpp"

#include <algorithm>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_format.h"
#include "modle/extr_barrier.hpp"
#include "modle/lefs.hpp"
#include "modle/parsers.hpp"

namespace modle {

DNA::Bin::Bin(uint32_t idx, uint64_t start, uint64_t end,
              const std::vector<ExtrusionBarrier>& barriers)
    : _idx(idx),
      _start(start),
      _end(end),
      _extr_barriers(std::make_unique<absl::InlinedVector<ExtrusionBarrier, 3>>(barriers.begin(),
                                                                                barriers.end())) {}

DNA::Bin::Bin(uint32_t idx, uint64_t start, uint64_t end)
    : _idx(idx), _start(start), _end(end), _extr_barriers(nullptr) {}

uint32_t DNA::Bin::get_start() const { return this->_start; }
uint32_t DNA::Bin::get_end() const { return this->_end; }
uint32_t DNA::Bin::size() const { return this->get_end() - this->get_start(); }
uint32_t DNA::Bin::get_n_extr_units() const {
  return this->_extr_units ? this->_extr_units->size() : 0;
}
uint32_t DNA::Bin::get_n_extr_barriers() const {
  return this->_extr_barriers ? this->_extr_barriers->size() : 0;
}
uint32_t DNA::Bin::get_index() const { return this->_idx; }

bool DNA::Bin::has_extr_barrier() const { return this->_extr_barriers != nullptr; }

ExtrusionBarrier* DNA::Bin::get_next_extr_barrier(ExtrusionBarrier* b, Direction d) const {
  assert(d != DNA::Direction::none);
  if (!this->_extr_barriers) return nullptr;  // There are no barriers bound to this Bin
  auto* b_itr = b ? this->_extr_barriers->begin() + std::distance(this->_extr_barriers->data(), b)
                  : this->_extr_barriers->begin();  // Convert Bin* to iterator

  auto barrier =
      std::find_if(b_itr, this->_extr_barriers->end(), [&d](const ExtrusionBarrier& barr) {
        return d == DNA::Direction::both || barr.get_direction() == d;
      });
  // Found an ExtrusionBarrier matching the search criteria
  if (barrier != this->_extr_barriers->end()) return barrier;
  return nullptr;  // Unable to find a suitable ExtrusionBarrier
}

ExtrusionBarrier* DNA::Bin::get_prev_extr_barrier(ExtrusionBarrier* b, Direction d) const {
  assert(d != DNA::Direction::none);
  if (!this->_extr_barriers) return nullptr;  // There are no barriers bound to this Bin
  auto b_itr = b ? this->_extr_barriers->rbegin() + std::distance(b, this->_extr_barriers->data())
                 : this->_extr_barriers->rbegin();  // Convert Bin* to iterator

  auto barrier =
      std::find_if(b_itr, this->_extr_barriers->rend(), [&d](const ExtrusionBarrier& barr) {
        return d == DNA::Direction::both || barr.get_direction() == d;
      });
  // Found an ExtrusionBarrier matching the search criteria
  if (barrier != this->_extr_barriers->rend()) return barrier.operator->();
  return nullptr;  // Unable to find a suitable ExtrusionBarrier
}

absl::InlinedVector<ExtrusionBarrier, 3>* DNA::Bin::get_all_extr_barriers() const {
  return this->_extr_barriers.get();
}

absl::InlinedVector<ExtrusionUnit*, 3>& DNA::Bin::get_extr_units() { return *this->_extr_units; }

void DNA::Bin::add_extr_barrier(ExtrusionBarrier b) {
  absl::FPrintF(stderr, "strand=%c\n", b.get_direction() == DNA::Direction::fwd ? '+' : '-');
  //  usleep(5e5);
  if (!this->_extr_barriers) {  // If this is the first ExtrusionBarrier, allocate the std::vector
    this->_extr_barriers = std::make_unique<absl::InlinedVector<ExtrusionBarrier, 3>>(
        absl::InlinedVector<ExtrusionBarrier, 3>{b});
  } else {
    this->_extr_barriers->emplace_back(b);
  }
}

void DNA::Bin::remove_extr_barrier(Direction d) {
  assert(d != Direction::none);
  if (!this->_extr_barriers) {
    throw std::logic_error(
        "Attempt to remove an extrusion barrier from a bin that doesn't have any!");
  }
  auto barrier = std::find_if(this->_extr_barriers->begin(), this->_extr_barriers->end(),
                              [&d](const ExtrusionBarrier& b) {
                                return d == DNA::Direction::both || b.get_direction() == d;
                              });
  if (barrier != this->_extr_barriers->end()) {
    this->_extr_barriers->erase(barrier);
    // Deallocate vector if there are no ExtrusionBarriers left
    if (this->_extr_barriers->empty()) this->_extr_barriers = nullptr;
    return;
  }
  throw std::logic_error(absl::StrFormat("Bin doesn't have barriers blocking %s!",
                                         d == Direction::fwd ? "forward" : "reverse"));
}

uint32_t DNA::Bin::add_extr_unit_binding(ExtrusionUnit* const unit) {
  if (this->_extr_units == nullptr) {  // If if there are no ExtrusionUnits binding to this Bin,
                                       // allocate the vector
    this->_extr_units = std::make_unique<absl::InlinedVector<ExtrusionUnit*, 3>>();
  }
  this->_extr_units->push_back(unit);
  return this->_extr_units->size();
}

uint32_t DNA::Bin::remove_extr_unit_binding(ExtrusionUnit* const unit) {
  assert(this->_extr_units != nullptr);
  if (this->_extr_units->size() == 1) {
    // Deallocate vector if unit is the last ExtrusionUnit binding to this Bin
    this->_extr_units = nullptr;
    return 0;
  }
  this->_extr_units->erase(std::remove(this->_extr_units->begin(), this->_extr_units->end(), unit),
                           this->_extr_units->end());
  return this->_extr_units->size();
}

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

uint32_t DNA::Bin::remove_all_extr_barriers() {
  uint32_t n = this->_extr_barriers->size();
  this->_extr_barriers = nullptr;
  return n;
}

void DNA::Bin::sort_extr_barriers_by_pos() {
  if (this->_extr_barriers && this->_extr_barriers->size() > 1) {
    std::sort(this->_extr_barriers->begin(), this->_extr_barriers->end());
  }
}

void DNA::Bin::add_extr_barrier(uint64_t pos, double prob_of_barrier_block,
                                DNA::Direction direction) {
  if (!this->_extr_barriers) {
    this->_extr_barriers = std::make_unique<absl::InlinedVector<ExtrusionBarrier, 3>>(
        absl::InlinedVector<ExtrusionBarrier, 3>{{pos, prob_of_barrier_block, direction}});
  } else {
    this->_extr_barriers->emplace_back(pos, prob_of_barrier_block, direction);
  }
}

DNA::DNA(uint64_t length, uint32_t bin_size)
    : _bins(make_bins(length, bin_size)), _length(length), _bin_size(bin_size) {}

void DNA::add_extr_barrier(ExtrusionBarrier& b, uint32_t pos) {
  assert(b.get_prob_of_block() >= 0 && b.get_prob_of_block() <= 1);
  assert(pos <= this->length());
  assert(pos / this->get_bin_size() < this->get_n_bins());
  this->_bins[pos / this->get_bin_size()].add_extr_barrier(b);
}

void DNA::add_extr_barrier(const BED& record) {
  assert(record.score >= 0 || record.score <= 1);
  assert(record.chrom_end <= this->length());
  assert(record.chrom_end / this->get_bin_size() < this->get_n_bins());
  assert(record.strand == '+' || record.strand == '-');
  uint64_t pos = (record.chrom_end + record.chrom_start) / 2;
  auto strand = record.strand == '+' ? DNA::Direction::fwd : DNA::Direction::rev;
  this->_bins[pos / this->get_bin_size()].add_extr_barrier(pos, record.score, strand);
}

void DNA::remove_extr_barrier(uint32_t pos, Direction direction) {
  assert(pos <= this->length());
  assert(pos / this->get_bin_size() < this->get_n_bins());
  this->_bins[pos / this->get_bin_size()].remove_extr_barrier(direction);
}

std::vector<DNA::Bin> DNA::make_bins(uint64_t length, uint32_t bin_size) {
  std::vector<DNA::Bin> bins;
  if (length <= bin_size) {  // Deal with short DNA molecules
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

uint32_t DNA::get_bin_size() const { return this->_bin_size; }

uint32_t DNA::get_n_barriers() const {
  return std::accumulate(this->_bins.begin(), this->_bins.end(), 0UL,
                         [](uint32_t accumulator, const DNA::Bin& b) {
                           if (b._extr_barriers) accumulator += b._extr_barriers->size();
                           return accumulator;
                         });
}

Chromosome::Chromosome(std::string name, uint64_t length, uint32_t bin_size,
                       uint32_t diagonal_width)
    : name(std::move(name)),
      dna(length, bin_size),
      contacts(diagonal_width / bin_size, length / bin_size) {}

uint32_t Chromosome::length() const { return this->dna.length(); }
uint32_t Chromosome::get_n_bins() const { return this->dna.get_n_bins(); }
uint32_t Chromosome::get_n_barriers() const { return this->barriers.size(); }

void Chromosome::sort_barriers_by_pos() {
  for (auto& bin : this->dna) {
    bin.sort_extr_barriers_by_pos();
    if (bin.has_extr_barrier()) {
      for (auto& b : *bin.get_all_extr_barriers()) {
        this->barriers.push_back(&b);
      }
    }
  }
}

void Chromosome::write_contacts_to_tsv(std::string_view chr_name, std::string_view output_dir,
                                       bool write_full_matrix) const {
  const auto t0 = absl::Now();
  //  uint32_t bytes_in, bytes_out;
  if (write_full_matrix) {
    auto file = absl::StrFormat("%s/%s.tsv.bz2", output_dir, chr_name);
    absl::FPrintF(stderr, "Writing full contact matrix for '%s' to file '%s'...\n", chr_name, file);
    auto [bytes_in, bytes_out] = this->contacts.write_full_matrix_to_tsv(file);
    absl::FPrintF(stderr,
                  "DONE writing '%s' in %s! Compressed size: %.2f MB (compression ration %.2fx)\n",
                  file, absl::FormatDuration(absl::Now() - t0), bytes_out / 1.0e6,
                  static_cast<double>(bytes_in) / bytes_out);
  } else {
    auto file = absl::StrFormat("%s/%s_raw.tsv.bz2", output_dir, chr_name);
    absl::FPrintF(stderr, "Writing raw contact matrix for '%s' to file '%s'...\n", chr_name, file);
    auto [bytes_in, bytes_out] = this->contacts.write_to_tsv(file);
    absl::FPrintF(stderr,
                  "DONE writing '%s' in %s! Compressed size: %.2f MB (compression ration %.2fx)\n",
                  file, absl::FormatDuration(absl::Now() - t0), bytes_out / 1.0e6,
                  static_cast<double>(bytes_in) / bytes_out);
  }
}

}  // namespace modle