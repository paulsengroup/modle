#pragma once

#include <absl/container/inlined_vector.h>
#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration, Time
#include <fmt/format.h>       // for system_error
#include <fmt/ostream.h>      // for print

#include <algorithm>                         // for min, find_if, sort
#include <boost/iostreams/filter/bzip2.hpp>  // for basic_bzip2_compressor, bzip2_error, bzip2_compressor
#include <boost/iostreams/filtering_stream.hpp>  // for filtering_ostream
#include <cassert>
#include <cerrno>
#include <cmath>       // for llround
#include <cstdint>     // for uint*_t
#include <cstdio>      // for stderr
#include <filesystem>  // for exists, weakly_canonical, path
#include <functional>  // for hash
#include <limits>
#include <numeric>      // for accumulate
#include <sstream>      // streamsize, basic_ios, ofstream, size_t, ios_base
#include <stdexcept>    // for logic_error, runtime_error
#include <string_view>  // for format, print, string_view
#include <type_traits>  // for declval
#include <utility>      // for move

#include "modle/bed.hpp"  // for BED
#include "modle/common.hpp"
#include "modle/extr_barrier.hpp"  // for ExtrusionBarrier
#include "modle/suppress_compiler_warnings.hpp"

namespace modle {

DISABLE_WARNING_PUSH
DISABLE_WARNING_CONVERSION
template <typename I1, typename I2>
DNA::DNA(I1 length, I2 bin_size, bool allocate_bins) : _length(length), _bin_size(bin_size) {
  DISABLE_WARNING_POP
  static_assert(std::is_integral<I1>::value && std::is_integral<I2>::value,
                "I1 and I2 should be an integral numeric type.");

  if (allocate_bins) {
    this->allocate_bins();
  }

#ifndef NDEBUG
  if (length >= std::numeric_limits<decltype(this->_length)>::max()) {
    std::runtime_error(fmt::format(
        FMT_STRING(
            "DNA::DNA(simulated_length={}, bin_size={}): Overflow detected: unable to represent {} "
            "in the range {}-{}"),
        length, bin_size, length, std::numeric_limits<decltype(this->_length)>::min(),
        std::numeric_limits<decltype(this->_length)>::max()));
  }
  if (length >= std::numeric_limits<decltype(this->_length)>::max()) {
    std::runtime_error(fmt::format(
        FMT_STRING(
            "DNA::DNA(simulated_length={}, bin_size={}): Overflow detected: unable to represent {} "
            "in the range {}-{}"),
        length, bin_size, bin_size, std::numeric_limits<decltype(this->_bin_size)>::min(),
        std::numeric_limits<decltype(this->_bin_size)>::max()));
  }
#endif
}

void DNA::allocate_bins() {
  if (this->_bins.empty()) {
    this->_bins = make_bins(this->_length, this->_bin_size);
  }
}

template <typename I1, typename I2>
void validate_params(std::string_view func_signature, I1 start, I2 end, const DNA::Bin* const bin) {
  if (start < std::numeric_limits<decltype(bin->get_start())>::min() ||
      start > std::numeric_limits<decltype(bin->get_start())>::max()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("{}: Overflow detected: start={} cannot be represented in the range {}-{}"),
        func_signature, start, std::numeric_limits<decltype(bin->get_start())>::min(),
        std::numeric_limits<decltype(bin->get_start())>::max()));
  }
  if (end < std::numeric_limits<decltype(bin->get_end())>::min() ||
      end > std::numeric_limits<decltype(bin->get_end())>::max()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("{}: Overflow detected: end={} cannot be represented in the range {}-{}"),
        func_signature, end, std::numeric_limits<decltype(bin->get_end())>::min(),
        std::numeric_limits<decltype(bin->get_end())>::max()));
  }
  if (start >= end) {
    throw std::logic_error(
        fmt::format(FMT_STRING("{}: start coordinate should be smaller than the end one: {} >= {}"),
                    func_signature, start, end));
  }
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_CONVERSION
template <typename I>
DNA::Bin::Bin(std::size_t idx, I start, I end, const std::vector<ExtrusionBarrier>& barriers)
    : _idx(idx),
      _start(start),
      _end(end),
      _extr_barriers(std::make_unique<absl::InlinedVector<ExtrusionBarrier, 3>>(barriers.begin(),
                                                                                barriers.end())) {
  DISABLE_WARNING_POP
  static_assert(std::is_integral_v<I>, "I should be an integral numeric type.");

#ifndef NDEBUG
  validate_params(fmt::format(FMT_STRING("DNA::Bin::Bin(idx={}, start={}, end={}, const "
                                         "std::vector<ExtrusionBarrier>& barriers)"),
                              idx, start, end),
                  start, end, this);
#endif
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_CONVERSION
template <typename I>
DNA::Bin::Bin(std::size_t idx, I start, I end)
    : _idx(idx), _start(start), _end(end), _extr_barriers(nullptr) {
  DISABLE_WARNING_POP
#ifndef NDEBUG
  validate_params(
      fmt::format(FMT_STRING("DNA::Bin::Bin(idx={}, start={}, end={})"), idx, start, end), start,
      end, this);
#endif
}

DNA::Bin::Bin(const Bin& other)
    : _idx(other._idx),
      _start(other._start),
      _end(other._end),
      _lef_affinity(other._lef_affinity),
      _extr_barriers(
          std::make_unique<absl::InlinedVector<ExtrusionBarrier, 3>>(*other._extr_barriers)),
      _extr_units(std::make_unique<absl::InlinedVector<ExtrusionUnit*, 3>>(*other._extr_units)) {}

DNA::Bin& DNA::Bin::operator=(const Bin& other) {
  if (&other == this) {
    return *this;
  }
  _idx = other._idx;
  _start = other._start;
  _end = other._end;
  _lef_affinity = other._lef_affinity;
  _extr_barriers =
      std::make_unique<absl::InlinedVector<ExtrusionBarrier, 3>>(*other._extr_barriers);
  _extr_units = std::make_unique<absl::InlinedVector<ExtrusionUnit*, 3>>(*other._extr_units);
  return *this;
}

uint32_t DNA::Bin::get_start() const { return this->_start; }
uint32_t DNA::Bin::get_end() const { return this->_end; }
uint64_t DNA::Bin::size() const { return this->get_end() - this->get_start(); }
uint64_t DNA::Bin::get_n_extr_units() const {
  return this->_extr_units ? this->_extr_units->size() : 0UL;
}
uint64_t DNA::Bin::get_n_extr_barriers() const {
  return this->_extr_barriers ? this->_extr_barriers->size() : 0UL;
}
std::size_t DNA::Bin::get_index() const { return this->_idx; }

bool DNA::Bin::has_extr_barrier() const { return this->_extr_barriers != nullptr; }

ExtrusionBarrier* DNA::Bin::get_ptr_to_next_extr_barrier(ExtrusionBarrier* b,
                                                         dna::Direction d) const {
  assert(d != dna::Direction::none);  // NOLINT;
  if (!this->_extr_barriers) {
    return nullptr;  // There are no barriers bound to this Bin
  }

  if (b == nullptr) {
    return &this->_extr_barriers->front();  // Return the first extr. barrier
  }
  assert(b >= &this->_extr_barriers->front() && b <= &this->_extr_barriers->back());  // NOLINT;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  for (++b; b < &this->_extr_barriers->back(); ++b) {
    if (d == dna::Direction::both || b->get_direction_of_block() == d) {
      return b;  // Found a suitable extr. barrier
    }
  }
  return nullptr;  // Unable to find a suitable extr. barrier
}

ExtrusionBarrier* DNA::Bin::get_ptr_to_prev_extr_barrier(ExtrusionBarrier* b,
                                                         dna::Direction d) const {
  assert(d != dna::Direction::none);  // NOLINT;
  if (!this->_extr_barriers) {
    return nullptr;  // There are no barriers bound to this Bin
  }

  if (b == nullptr) {
    return &this->_extr_barriers->back();  // Return the last extr. barrier
  }

  assert(b >= &this->_extr_barriers->front() && b <= &this->_extr_barriers->back());  // NOLINT;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  for (--b; b > &this->_extr_barriers->front(); --b) {
    if (d == dna::Direction::both || b->get_direction_of_block() == d) {
      return b;  // Found a suitable extr. barrier
    }
  }
  return nullptr;  // Unable to find a suitable extr. barrier
}

ExtrusionBarrier* DNA::Bin::get_ptr_to_next_extr_barrier(uint64_t pos, dna::Direction d) const {
  if (!this->_extr_barriers) {
    return nullptr;  // There are no barriers bound to this Bin
  }
  switch (d) {
    case dna::Direction::both:
      return &this->_extr_barriers->front();
    case dna::Direction::fwd: {
      const auto& barrier = std::find_if(this->_extr_barriers->begin(), this->_extr_barriers->end(),
                                         [&](const auto& b) { return b.get_pos() >= pos; });
      return barrier != this->_extr_barriers->end() ? &(*barrier) : nullptr;
    }
    case dna::Direction::rev: {
      const auto& barrier = std::find_if(this->_extr_barriers->begin(), this->_extr_barriers->end(),
                                         [&](const auto& b) { return b.get_pos() <= pos; });
      return barrier != this->_extr_barriers->end() ? &(*barrier) : nullptr;
    }
    default:
      throw std::logic_error(
          fmt::format("DNA::Bin::get_ptr_to_next_extr_barrier was called with d = {}. Allowed "
                      "directions are ::fwd, ::rev and ::both",
                      d));
  }
}

absl::InlinedVector<ExtrusionUnit*, 3>& DNA::Bin::get_extr_units() { return *this->_extr_units; }

float DNA::Bin::get_lef_affinity() const { return this->_lef_affinity; }

ExtrusionBarrier* DNA::Bin::add_extr_barrier(ExtrusionBarrier b) {
  if (!this->_extr_barriers) {  // If this is the first ExtrusionBarrier, allocate the std::vector
    this->_extr_barriers = std::make_unique<absl::InlinedVector<ExtrusionBarrier, 3>>(
        absl::InlinedVector<ExtrusionBarrier, 3>{std::move(b)});
  } else {
    this->_extr_barriers->emplace_back(b);
  }
  return &this->_extr_barriers->back();
}

void DNA::Bin::remove_extr_barrier(dna::Direction d) {
  assert(d != dna::Direction::none);  // NOLINT;
  if (!this->_extr_barriers) {
    throw std::logic_error(
        "Attempt to remove an extrusion barrier from a bin that doesn't have any!");
  }
  auto barrier = std::find_if(this->_extr_barriers->begin(), this->_extr_barriers->end(),
                              [&d](const ExtrusionBarrier& b) {
                                return d == dna::Direction::both || b.get_direction_of_block() == d;
                              });
  if (barrier != this->_extr_barriers->end()) {
    this->_extr_barriers->erase(barrier);
    // Deallocate vector if there are no ExtrusionBarriers left
    if (this->_extr_barriers->empty()) {
      this->_extr_barriers = nullptr;
    }
    return;
  }
  throw std::logic_error(fmt::format("Bin #{} doesn't have barriers blocking in {} direction!",
                                     this->get_index(),
                                     d == dna::Direction::fwd ? "forward" : "reverse"));
}

uint64_t DNA::Bin::add_extr_unit_binding(ExtrusionUnit* const unit) {
  assert(unit);
  if (this->_extr_units == nullptr) {  // If if there are no ExtrusionUnits binding to this Bin,
    // allocate the vector
    this->_extr_units = std::make_unique<absl::InlinedVector<ExtrusionUnit*, 3>>();
  }
  this->_extr_units->push_back(unit);
  return this->_extr_units->size();
}

uint64_t DNA::Bin::remove_extr_unit_binding(ExtrusionUnit* const unit) {
  assert(this->_extr_units != nullptr);  // NOLINT;
  if (this->_extr_units->size() == 1) {
    // Deallocate vector if unit is the last ExtrusionUnit binding to this Bin
    this->_extr_units = nullptr;
    return 0;
  }
  this->_extr_units->erase(std::remove(this->_extr_units->begin(), this->_extr_units->end(), unit),
                           this->_extr_units->end());
  return this->_extr_units->size();
}

DNA::Bin* DNA::get_ptr_to_bin_from_pos(uint64_t pos) {
#ifndef NDEBUG
  assert(!this->_bins.empty());
  if (pos > this->length()) {
    throw std::logic_error(fmt::format(
        "DNA::get_ptr_to_bin_from_pos(pos={}): pos > this->simulated_length(): {} > {}\n", pos, pos,
        this->length()));
  }
  if ((pos / this->get_bin_size()) > this->get_n_bins()) {
    throw std::logic_error(
        fmt::format("(pos / this->bin_size()) >= this->nbins(): ({} / {}) >= {}\n", pos,
                    this->get_bin_size(), this->get_n_bins()));
  }
#endif
  return &this->_bins[pos / this->get_bin_size()];
}

DNA::Bin& DNA::get_bin_from_pos(uint64_t pos) {
  assert(!this->_bins.empty());
  return *this->get_ptr_to_bin_from_pos(pos);
}

std::vector<DNA::Bin>::iterator DNA::begin() {
  assert(!this->_bins.empty());
  return this->_bins.begin();
}
std::vector<DNA::Bin>::iterator DNA::end() {
  assert(!this->_bins.empty());
  return this->_bins.end();
}
std::vector<DNA::Bin>::const_iterator DNA::cbegin() const {
  assert(!this->_bins.empty());
  return this->_bins.cbegin();
}
std::vector<DNA::Bin>::const_iterator DNA::cend() const {
  assert(!this->_bins.empty());
  return this->_bins.cend();
}

DNA::Bin* DNA::get_ptr_to_prev_bin(const DNA::Bin& current_bin) {
  assert(!this->_bins.empty());
  assert(current_bin._idx > 0);  // NOLINT;
  return &this->_bins[current_bin._idx - 1];
}

DNA::Bin& DNA::get_prev_bin(const Bin& current_bin) {
  assert(!this->_bins.empty());
  return *this->get_ptr_to_prev_bin(current_bin);
}

DNA::Bin& DNA::get_next_bin(const DNA::Bin& current_bin) {
  assert(!this->_bins.empty());
  return *this->get_ptr_to_next_bin(current_bin);
}

DNA::Bin* DNA::get_ptr_to_next_bin(const DNA::Bin& current_bin) {
  assert(!this->_bins.empty());
  assert(current_bin._idx < this->get_n_bins());  // NOLINT;
  return &this->_bins[current_bin._idx + 1];
}

DNA::Bin& DNA::get_first_bin() {
  assert(!this->_bins.empty());
  return *this->get_ptr_to_first_bin();
}

DNA::Bin* DNA::get_ptr_to_first_bin() {
  assert(!this->_bins.empty());  // NOLINT;
  return &this->_bins.front();
}

DNA::Bin& DNA::get_last_bin() { return *this->get_ptr_to_last_bin(); }

DNA::Bin* DNA::get_ptr_to_last_bin() {
  assert(!this->_bins.empty());  // NOLINT;
  return &this->_bins.back();
}

double DNA::get_total_lef_affinity() const {
  assert(!this->_bins.empty());
  return std::accumulate(this->_bins.begin(), this->_bins.end(), 0.0,
                         [](double accumulator, const auto& bin) {
                           return accumulator + static_cast<double>(bin.get_lef_affinity());
                         });
}

uint64_t DNA::Bin::remove_all_extr_barriers() {
  auto n = this->_extr_barriers->size();
  this->_extr_barriers = nullptr;
  return n;
}

void DNA::Bin::sort_extr_barriers_by_pos() {
  if (this->_extr_barriers && this->_extr_barriers->size() > 1) {
    std::sort(this->_extr_barriers->begin(), this->_extr_barriers->end());
  }
}

ExtrusionBarrier* DNA::Bin::add_extr_barrier(uint64_t pos, double prob_of_barrier_block,
                                             dna::Direction direction) {
  if (!this->_extr_barriers) {
    this->_extr_barriers = std::make_unique<absl::InlinedVector<ExtrusionBarrier, 3>>(
        absl::InlinedVector<ExtrusionBarrier, 3>{{pos, prob_of_barrier_block, direction}});
  } else {
    this->_extr_barriers->emplace_back(pos, prob_of_barrier_block, direction);
  }
  return &this->_extr_barriers->back();
}

ExtrusionBarrier* DNA::add_extr_barrier(ExtrusionBarrier& b, uint32_t pos) {
  assert(!this->_bins.empty());
  assert(b.get_prob_of_block() >= 0 && b.get_prob_of_block() <= 1);  // NOLINT;
  assert(pos <= this->length());                                     // NOLINT;
  assert(pos / this->get_bin_size() < this->get_n_bins());           // NOLINT;
  return this->_bins[pos / this->get_bin_size()].add_extr_barrier(b);
}

ExtrusionBarrier* DNA::add_extr_barrier(const modle::bed::BED& record) {
  assert(!this->_bins.empty());
  assert(record.score >= 0 && record.score <= 1);                        // NOLINT;
  assert(record.chrom_end <= this->length());                            // NOLINT;
  assert(record.chrom_end / this->get_bin_size() < this->get_n_bins());  // NOLINT;
  assert(record.strand == '+' || record.strand == '-');                  // NOLINT;
  const auto pos = static_cast<uint64_t>(
      std::llround(static_cast<double>(record.chrom_end + record.chrom_start) / 2.0));
  dna::Direction strand = record.strand == '+' ? dna::Direction::fwd : dna::Direction::rev;
  return this->_bins[pos / this->get_bin_size()].add_extr_barrier(pos, record.score, strand);
}

void DNA::remove_extr_barrier(uint32_t pos, dna::Direction direction) {
  assert(!this->_bins.empty());
  assert(pos <= this->length());                            // NOLINT;
  assert(pos / this->get_bin_size() < this->get_n_bins());  // NOLINT;
  this->_bins[pos / this->get_bin_size()].remove_extr_barrier(direction);
}

std::vector<DNA::Bin> DNA::make_bins(uint64_t length, uint32_t bin_size) {
  using N = decltype(DNA::Bin::_end);
  assert(length != 0);
  assert(length < std::numeric_limits<N>::max());
  assert(bin_size != 0);
  if (length <= bin_size) {  // Deal with short DNA molecules
    return {{0UL, 0UL, length}};
  }

  std::vector<DNA::Bin> bins((length / bin_size) + (length % bin_size != 0));
  std::generate(bins.begin(), bins.end(), [pos = 0UL, i = 0UL, bin_size]() mutable {
    pos += bin_size;
    return DNA::Bin(i++, pos - bin_size, pos);
  });
  assert(bins.back()._start < length);
  bins.back()._end = static_cast<N>(length);
  return bins;
}

uint64_t DNA::length() const { return this->_length; }

uint64_t DNA::get_n_bins() const {
  assert(!this->_bins.empty());
  return this->_bins.size();
}

uint32_t DNA::get_bin_size() const { return this->_bin_size; }

uint64_t DNA::get_n_barriers() const {
  assert(!this->_bins.empty());
  return std::accumulate(this->_bins.begin(), this->_bins.end(), 0UL,
                         [](uint64_t accumulator, const DNA::Bin& b) {
                           if (b._extr_barriers) {
                             accumulator += b._extr_barriers->size();
                           }
                           return accumulator;
                         });
}

Chromosome::Chromosome(std::string chr_name, uint64_t length, uint32_t bin_size,
                       uint32_t diagonal_width_, uint64_t seed, bool allocate)
    : name(std::move(chr_name)),
      start(0),
      end(length),
      total_length(length),
      diagonal_width(diagonal_width_),
      dna(length, bin_size),
      _seed(seed + std::hash<std::string>{}(this->name) +
            std::hash<uint64_t>{}(this->simulated_length())) {
  if (allocate) {
    this->allocate();
  }

  std::seed_seq seeder{this->_seed};
  this->_rand_eng = std::mt19937(seeder);
}

Chromosome::Chromosome(std::string chr_name, uint64_t chr_start, uint64_t chr_end, uint64_t length,
                       uint32_t bin_size, uint32_t diagonal_width_, uint64_t seed, bool allocate)
    : name(std::move(chr_name)),
      start(chr_start),
      end(chr_end),
      total_length(length),
      diagonal_width(diagonal_width_),
      dna(chr_end - chr_start, bin_size),
      contacts(std::min(static_cast<uint64_t>(diagonal_width_ / bin_size), this->get_n_bins()),
               this->get_n_bins()),
      _seed(seed + std::hash<std::string>{}(this->name) +
            std::hash<uint64_t>{}(this->simulated_length())) {
  if (allocate) {
    this->allocate();
  }

  std::seed_seq seeder{this->_seed};
  this->_rand_eng = std::mt19937(seeder);
}

uint64_t Chromosome::simulated_length() const { return this->end - this->start; }
uint64_t Chromosome::real_length() const { return this->total_length; }
uint64_t Chromosome::get_start_pos() const { return this->start; }
uint64_t Chromosome::get_end_pos() const { return this->end; }
uint64_t Chromosome::get_n_bins() const { return this->dna.get_n_bins(); }
uint32_t Chromosome::get_bin_size() const { return this->dna.get_bin_size(); }
uint64_t Chromosome::get_n_barriers() const { return this->barriers.size(); }
uint64_t Chromosome::get_n_lefs() const { return this->lefs.size(); }
double Chromosome::get_total_lef_affinity() const { return this->dna.get_total_lef_affinity(); }

void Chromosome::sort_barriers_by_pos() {
  for (auto& bin : this->dna) {
    bin.sort_extr_barriers_by_pos();
    if (bin.has_extr_barrier()) {
      for (auto* b = bin.get_ptr_to_next_extr_barrier(nullptr); b != nullptr;
           b = bin.get_ptr_to_next_extr_barrier(b)) {
        this->barriers.push_back(b);
      }
    }
  }
}

void Chromosome::write_contacts_to_tsv(std::string_view output_dir, bool force_overwrite) const {
  auto t0 = absl::Now();
  std::filesystem::create_directories(output_dir);
  std::string path_to_outfile = std::filesystem::weakly_canonical(
      fmt::format("{}/{}_modle_cmatrix.tsv.bz2", output_dir, this->name));
  fmt::print(stderr, "Writing contact matrix for '{}' to file '{}'...\n", this->name,
             path_to_outfile);
  if (!force_overwrite && std::filesystem::exists(path_to_outfile)) {
    fmt::print(stderr, "File '{}' already exists, SKIPPING! Pass --force to overwrite...\n",
               path_to_outfile);
    return;
  }

  const std::string header =
      fmt::format("#{}\t{}\t{}\t{}\t{}\n", this->name, this->get_bin_size(), this->start, this->end,
                  this->contacts.nrows() * this->get_bin_size());

  auto [bytes_in, bytes_out] = this->contacts.write_to_tsv(path_to_outfile, header);
  fmt::print(stderr,
             "DONE writing '{}' in {}! Compressed size: {:.2f} MB (compression ratio {:.2f}x)\n",
             // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
             path_to_outfile, absl::FormatDuration(absl::Now() - t0), bytes_out / 1.0e6,
             static_cast<double>(bytes_in) / bytes_out);
}

void Chromosome::write_barriers_to_tsv(std::string_view output_dir, bool force_overwrite) const {
  auto path_to_outfile = std::filesystem::weakly_canonical(
      fmt::format("{}/{}.extrusion_barriers.tsv.bz2", output_dir, this->name));
  fmt::print(stderr, "Writing extr. barrier coordinates for '{}' to file '{}'\n", this->name,
             path_to_outfile);
  if (!force_overwrite && std::filesystem::exists(path_to_outfile)) {
    fmt::print(stderr, "File '{}' already exists. Pass --force to overwrite... SKIPPING.\n",
               path_to_outfile);
    return;
  }
  try {
    std::ofstream fp(path_to_outfile, std::ios_base::binary);
    if (!fp) {
      throw fmt::system_error(errno, "Unable to open file '{}' for writing", path_to_outfile);
    }
    boost::iostreams::filtering_ostream out;
    std::string buff;
    out.push(boost::iostreams::bzip2_compressor());
    out.push(fp);
    if (!out) {
      throw std::runtime_error(
          fmt::format("An error occurred while initializing the compression stream for file '{}'",
                      path_to_outfile));
    }

    for (const auto& barrier : this->barriers) {
      fmt::print(out, "{}{}\n",
                 barrier->get_direction_of_block() == dna::Direction::fwd ? "+" : "-",
                 barrier->get_pos());
      if (!out || !fp) {
        throw fmt::system_error(errno, "IO error while writing to file '{}'", path_to_outfile);
      }
    }
  } catch (const boost::iostreams::bzip2_error& err) {
    throw std::runtime_error(fmt::format("An error occurred while compressing file '{}': {}",
                                         path_to_outfile, err.what()));
  }
}

void Chromosome::allocate_contacts() {
  using I =
      std::remove_const_t<std::remove_pointer_t<decltype(contacts.get_raw_count_vector().data())>>;
  if (contacts.nrows() == 0 && contacts.nrows() == 0) {
    const auto nrows = std::min(diagonal_width / dna.get_bin_size(), this->get_n_bins());
    const auto ncols = this->get_n_bins();
    contacts = ContactMatrix<I>(nrows, ncols);
  }
}

void Chromosome::allocate() {
  this->allocate_contacts();
  this->dna.allocate_bins();
}

}  // namespace modle
