#pragma once

#include <absl/container/inlined_vector.h>
#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration, Time
#include <absl/types/span.h>
#include <fmt/format.h>   // for system_error
#include <fmt/ostream.h>  // for print

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
#ifdef USE_XOSHIRO
#include <XoshiroCpp.hpp>
#endif

#include "modle/bed.hpp"  // for BED
#include "modle/common.hpp"
#include "modle/extr_barrier.hpp"  // for ExtrusionBarrier
#include "modle/suppress_compiler_warnings.hpp"

namespace modle {

// Pre-declarations
template <typename I1, typename I2>
void validate_params(std::string_view func_signature, I1 start, I2 end, const DNA::Bin* bin);

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

std::size_t DNA::size() const { return this->_length; }

uint32_t DNA::get_bin_size() const { return this->_bin_size; }

double DNA::get_total_lef_affinity() const {
  assert(!this->_bins.empty());
  return std::accumulate(this->_bins.begin(), this->_bins.end(), 0.0,
                         [](double accumulator, const auto& bin) {
                           return accumulator + static_cast<double>(bin.get_lef_affinity());
                         });
}

std::size_t DNA::get_n_bins() const {
  return !this->_bins.empty()
             ? this->_bins.size()
             : (this->_length / this->_bin_size) + (this->_length % this->_bin_size != 0);
}

std::size_t DNA::get_n_barriers() const {
  assert(!this->_bins.empty());
  return std::accumulate(
      this->_bins.begin(), this->_bins.end(), 0UL,
      [](uint64_t accumulator, const auto& b) { return accumulator + b._extr_barriers.size(); });
}

DNA::Bin* DNA::get_ptr_to_bin(std::size_t pos) {
#ifndef NDEBUG
  assert(!this->_bins.empty());
  if (pos > this->size()) {
    throw std::logic_error(fmt::format(
        FMT_STRING("DNA::get_ptr_to_bin(pos={}): pos > this->simulated_length(): {} > {}\n"), pos,
        pos, this->size()));
  }
  if ((pos / this->get_bin_size()) > this->get_n_bins()) {
    throw std::logic_error(
        fmt::format(FMT_STRING("(pos / this->bin_size()) >= this->nbins(): ({} / {}) >= {}\n"), pos,
                    this->get_bin_size(), this->get_n_bins()));
  }
#endif
  return &this->_bins[pos / this->get_bin_size()];
}

DNA::Bin& DNA::get_bin(std::size_t pos) {
  assert(!this->_bins.empty());
  return *this->get_ptr_to_bin(pos);
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

DNA::Bin* DNA::get_ptr_to_next_bin(const DNA::Bin& current_bin) {
  assert(!this->_bins.empty());
  assert(current_bin._idx < this->get_n_bins());  // NOLINT;
  return &this->_bins[current_bin._idx + 1];
}

DNA::Bin& DNA::get_next_bin(const DNA::Bin& current_bin) {
  assert(!this->_bins.empty());
  return *this->get_ptr_to_next_bin(current_bin);
}

DNA::Bin& DNA::get_first_bin() {
  assert(!this->_bins.empty());
  return *this->get_ptr_to_first_bin();
}

DNA::Bin* DNA::get_ptr_to_first_bin() {
  assert(!this->_bins.empty());
  return &this->_bins.front();
}

DNA::Bin& DNA::get_last_bin() { return *this->get_ptr_to_last_bin(); }

DNA::Bin* DNA::get_ptr_to_last_bin() {
  assert(!this->_bins.empty());
  return &this->_bins.back();
}

std::size_t DNA::get_bin_idx(std::size_t pos) const { return pos / this->_bin_size; }
std::size_t DNA::get_bin_idx(const DNA::Bin& bin) const { return bin.get_index(); }

DNA::Bin& DNA::operator[](std::size_t idx) {
  assert(idx < this->_bins.size());
  return this->_bins[idx];
}

const DNA::Bin& DNA::operator[](std::size_t idx) const {
  assert(idx < this->_bins.size());
  return this->_bins[idx];
}

std::vector<DNA::Bin>::iterator DNA::begin() { return this->_bins.begin(); }
std::vector<DNA::Bin>::iterator DNA::end() { return this->_bins.end(); }
std::vector<DNA::Bin>::const_iterator DNA::cbegin() const { return this->_bins.cbegin(); }
std::vector<DNA::Bin>::const_iterator DNA::cend() const { return this->_bins.cend(); }

ExtrusionBarrier* DNA::add_extr_barrier(ExtrusionBarrier& b, std::size_t pos) {
  assert(!this->_bins.empty());
  assert(b.get_prob_of_block() >= 0 && b.get_prob_of_block() <= 1);  // NOLINT;
  assert(pos <= this->size());                                       // NOLINT;
  return this->_bins[pos / this->get_bin_size()].add_extr_barrier(b);
}

ExtrusionBarrier* DNA::add_extr_barrier(const modle::bed::BED& record) {
  assert(!this->_bins.empty());
  assert(record.score >= 0 && record.score <= 1);                        // NOLINT;
  assert(record.chrom_end <= this->size());                              // NOLINT;
  assert(record.chrom_end / this->get_bin_size() < this->get_n_bins());  // NOLINT;
  assert(record.strand == '+' || record.strand == '-');                  // NOLINT;
  const auto pos = static_cast<uint64_t>(
      std::llround(static_cast<double>(record.chrom_end + record.chrom_start) / 2.0));
  dna::Direction strand = record.strand == '+' ? dna::Direction::fwd : dna::Direction::rev;
  return this->_bins[pos / this->get_bin_size()].add_extr_barrier(pos, record.score, strand);
}

void DNA::remove_extr_barrier(std::size_t pos, dna::Direction direction) {
  assert(!this->_bins.empty());
  assert(pos <= this->size());                              // NOLINT;
  assert(pos / this->get_bin_size() < this->get_n_bins());  // NOLINT;
  this->_bins[pos / this->get_bin_size()].remove_extr_barrier(direction);
}

void DNA::allocate_bins() {
  if (this->_bins.empty()) {
    this->_bins = make_bins(this->_length, this->_bin_size);
  }
}

std::vector<DNA::Bin> DNA::make_bins(std::size_t length, uint32_t bin_size) {
  using N = decltype(DNA::Bin::_bin_size);
  assert(length != 0);
  assert(length < std::numeric_limits<N>::max());
  assert(bin_size != 0);
  if (length <= bin_size) {  // Deal with short DNA molecules
    return {{0UL, length}};
  }

  std::vector<DNA::Bin> bins((length / bin_size) + (length % bin_size != 0));
  std::generate(bins.begin(), bins.end(), [pos = 0UL, i = 0UL, bin_size]() mutable {
    pos += bin_size;
    return DNA::Bin(i++, bin_size);
  });
  return bins;
}

// DNA::Bin

template <typename I>
DNA::Bin::Bin(std::size_t idx, I bin_size, absl::Span<ExtrusionBarrier> barriers)
    : _idx(idx),
      _bin_size(static_cast<decltype(_bin_size)>(bin_size)),
      _extr_barriers(barriers.begin(), barriers.end()) {
  static_assert(std::is_integral_v<I>, "bin_size should have an integral numeric type.");

#ifndef NDEBUG
  const auto start = _idx * bin_size;
  const auto end = start + bin_size;
  validate_params(fmt::format(FMT_STRING("DNA::Bin::Bin(idx={}, start={}, end={}, const "
                                         "std::vector<ExtrusionBarrier>& barriers)"),
                              idx, start, end),
                  start, end, this);
#endif
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_CONVERSION
template <typename I>
DNA::Bin::Bin(std::size_t idx, I bin_size) : _idx(idx), _bin_size(bin_size) {
  DISABLE_WARNING_POP
#ifndef NDEBUG
  const auto start = _idx * bin_size;
  const auto end = start + bin_size;
  validate_params(
      fmt::format(FMT_STRING("DNA::Bin::Bin(idx={}, start={}, end={})"), idx, start, end), start,
      end, this);
#endif
}

std::size_t DNA::Bin::size() const { return this->_bin_size; }
std::size_t DNA::Bin::get_n_extr_units() const { return this->_extr_units.size(); }
std::size_t DNA::Bin::get_n_extr_barriers() const { return this->_extr_barriers.size(); }

std::size_t DNA::Bin::get_index() const { return this->_idx; }
std::size_t DNA::Bin::get_start() const { return this->_idx * this->_bin_size; }
std::size_t DNA::Bin::get_end() const { return this->get_start() + this->_bin_size; }

bool DNA::Bin::has_extr_barrier() const { return !this->_extr_barriers.empty(); }
absl::Span<ExtrusionBarrier> DNA::Bin::get_all_extr_barriers() {
  return absl::MakeSpan(this->_extr_barriers);
}

ExtrusionBarrier* DNA::Bin::get_ptr_to_prev_extr_barrier(const ExtrusionBarrier* const old_b,
                                                         dna::Direction d) {
#ifndef NDEBUG
  if (d != dna::Direction::both && d != dna::Direction::fwd && d != dna::Direction::rev) {
    throw std::logic_error(

        fmt::format(
            FMT_STRING("DNA::Bin::get_ptr_to_prev_extr_barrier was called with d = {}. Allowed "
                       "directions are ::fwd, ::rev and ::both"),
            d));
  }
#endif
  if (this->_extr_barriers.empty()) {
    return nullptr;  // There are no barriers bound to this Bin
  }

  auto begin =
      this->_extr_barriers.rbegin() + (old_b ? &this->_extr_barriers.back() - (old_b - 1) : 0);
  auto end = this->_extr_barriers.rend();

  const auto& it = std::find_if(begin, end, [&](const auto& next_b) {
    return d == dna::Direction::both || next_b.get_direction_of_block() == d;
  });

  return it != end ? &(*it) : nullptr;  // Unable to find a suitable extr. barrier
}

ExtrusionBarrier* DNA::Bin::get_ptr_to_prev_extr_barrier(uint64_t pos, dna::Direction d) {
#ifndef NDEBUG
  if (d != dna::Direction::both && d != dna::Direction::fwd && d != dna::Direction::rev) {
    throw std::logic_error(

        fmt::format(
            FMT_STRING("DNA::Bin::get_ptr_to_prev_extr_barrier was called with d = {}. Allowed "
                       "directions are ::fwd, ::rev and ::both"),
            d));
  }
#endif

  if (this->_extr_barriers.empty()) {
    return nullptr;  // There are no barriers bound to this Bin
  }

  const auto& barrier =
      std::find_if(this->_extr_barriers.rbegin(), this->_extr_barriers.rend(), [&](const auto& b) {
        return (d == dna::Direction::both || b.get_direction_of_block() == d) && b.get_pos() <= pos;
      });

  return barrier != this->_extr_barriers.rend() ? &(*barrier) : nullptr;
}

ExtrusionBarrier* DNA::Bin::get_ptr_to_next_extr_barrier(const ExtrusionBarrier* const old_b,
                                                         dna::Direction d) {
#ifndef NDEBUG
  if (d != dna::Direction::both && d != dna::Direction::fwd && d != dna::Direction::rev) {
    throw std::logic_error(

        fmt::format(
            FMT_STRING("DNA::Bin::get_ptr_to_next_extr_barrier was called with d = {}. Allowed "
                       "directions are ::fwd, ::rev and ::both"),
            d));
  }
#endif

  if (this->_extr_barriers.empty()) {
    return nullptr;  // There are no barriers bound to this Bin
  }

  auto begin =
      this->_extr_barriers.begin() + (old_b ? (old_b + 1) - &this->_extr_barriers.front() : 0);
  auto end = this->_extr_barriers.end();

  const auto& it = std::find_if(begin, end, [&](const auto& next_b) {
    return d == dna::Direction::both || next_b.get_direction_of_block() == d;
  });

  return it != end ? &(*it) : nullptr;  // Unable to find a suitable extr. barrier
}

ExtrusionBarrier* DNA::Bin::get_ptr_to_next_extr_barrier(std::size_t pos, dna::Direction d) {
#ifndef NDEBUG
  if (d != dna::Direction::both && d != dna::Direction::fwd && d != dna::Direction::rev) {
    throw std::logic_error(

        fmt::format(
            FMT_STRING("DNA::Bin::get_ptr_to_next_extr_barrier was called with d = {}. Allowed "
                       "directions are ::fwd, ::rev and ::both"),
            d));
  }
#endif

  if (this->_extr_barriers.empty()) {
    return nullptr;  // There are no barriers bound to this Bin
  }

  const auto& barrier =
      std::find_if(this->_extr_barriers.begin(), this->_extr_barriers.end(), [&](const auto& b) {
        return (d == dna::Direction::both || b.get_direction_of_block() == d) && b.get_pos() >= pos;
      });

  return barrier != this->_extr_barriers.end() ? &(*barrier) : nullptr;
}

absl::Span<ExtrusionUnit*> DNA::Bin::get_extr_units() { return absl::MakeSpan(this->_extr_units); }
float DNA::Bin::get_lef_affinity() const { return this->_lef_affinity; }

ExtrusionBarrier* DNA::Bin::add_extr_barrier(ExtrusionBarrier b) {
  auto& ptr = this->_extr_barriers.emplace_back(std::move(b));
  return &ptr;
}
ExtrusionBarrier* DNA::Bin::add_extr_barrier(uint64_t pos, double prob_of_barrier_block,
                                             dna::Direction direction) {
  auto& ptr = this->_extr_barriers.emplace_back(pos, prob_of_barrier_block, direction);
  return &ptr;
}

void DNA::Bin::remove_extr_barrier(dna::Direction d) {
  assert(d != dna::Direction::none);  // NOLINT;
  if (this->_extr_barriers.empty()) {
    throw std::logic_error(
        "Attempt to remove an extrusion barrier from a bin that doesn't have any!");
  }
  auto barrier = std::find_if(this->_extr_barriers.begin(), this->_extr_barriers.end(),
                              [&d](const ExtrusionBarrier& b) {
                                return d == dna::Direction::both || b.get_direction_of_block() == d;
                              });
  if (barrier != this->_extr_barriers.end()) {
    this->_extr_barriers.erase(barrier);
    return;
  }
  throw std::logic_error(
      fmt::format(FMT_STRING("Bin #{} doesn't have barriers blocking in {} direction!"),
                  this->get_index(), d == dna::Direction::fwd ? "forward" : "reverse"));
}

uint64_t DNA::Bin::remove_all_extr_barriers() {
  auto tmp = decltype(this->_extr_barriers){};
  std::swap(this->_extr_barriers, tmp);
  return tmp.size();
}

uint64_t DNA::Bin::register_extr_unit_binding(ExtrusionUnit* const unit) {
  assert(unit);
  this->_extr_units.push_back(unit);
  return this->_extr_units.size();
}

uint64_t DNA::Bin::unregister_extr_unit_binding(ExtrusionUnit* const unit) {
  assert(unit);
  this->_extr_units.erase(std::remove(this->_extr_units.begin(), this->_extr_units.end(), unit),
                          this->_extr_units.end());
  return this->_extr_units.size();
}

void DNA::Bin::sort_extr_barriers_by_pos() {
  std::sort(this->_extr_barriers.begin(), this->_extr_barriers.end());
}

Chromosome::Chromosome(std::string chr_name, std::size_t length, uint32_t bin_size,
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

#ifdef USE_XOSHIRO
  modle::seeder seeder(this->_seed);
  this->_rand_eng = modle::PRNG(seeder.generateSeedSequence<4>());
#else
  modle::seeder seeder_{this->_seed};
  this->_rand_eng = modle::PRNG(seeder_);
#endif
}

Chromosome::Chromosome(std::string chr_name, std::size_t chr_start, std::size_t chr_end,
                       std::size_t length, uint32_t bin_size, uint32_t diagonal_width_,
                       uint64_t seed, bool allocate)
    : name(std::move(chr_name)),
      start(chr_start),
      end(chr_end),
      total_length(length),
      diagonal_width(diagonal_width_),
      dna(chr_end - chr_start, bin_size),
      _seed(seed + std::hash<std::string>{}(this->name) +
            std::hash<uint64_t>{}(this->simulated_length())) {
  if (allocate) {
    this->allocate();
  }

#ifdef USE_XOSHIRO
  modle::seeder seeder(this->_seed);
  this->_rand_eng = modle::PRNG(seeder.generateSeedSequence<4>());
#else
  modle::seeder seeder_{this->_seed};
  this->_rand_eng = modle::PRNG(seeder_);
#endif
}

std::size_t Chromosome::simulated_length() const { return this->end - this->start; }
std::size_t Chromosome::real_length() const { return this->total_length; }
std::size_t Chromosome::get_start_pos() const { return this->start; }
std::size_t Chromosome::get_end_pos() const { return this->end; }
std::size_t Chromosome::get_nbins() const { return this->dna.get_n_bins(); }
uint32_t Chromosome::get_bin_size() const { return this->dna.get_bin_size(); }
std::size_t Chromosome::get_nbarriers() const { return this->barriers.size(); }
std::size_t Chromosome::get_nlefs() const { return this->lefs.size(); }
double Chromosome::get_total_lef_affinity() const { return this->dna.get_total_lef_affinity(); }

/*
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
 */

void Chromosome::allocate_contacts() {
  using I =
      std::remove_const_t<std::remove_pointer_t<decltype(contacts.get_raw_count_vector().data())>>;
  if (contacts.nrows() == 0 && contacts.nrows() == 0) {
    const auto nrows = std::min(diagonal_width / dna.get_bin_size(), this->get_nbins());
    const auto ncols = this->get_nbins();
    contacts = ContactMatrix<I>(nrows, ncols);
  }
}

void Chromosome::allocate() {
  this->allocate_contacts();
  this->dna.allocate_bins();
}

bool Chromosome::ok() const { return this->_ok; }

// Internal functions

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

}  // namespace modle
