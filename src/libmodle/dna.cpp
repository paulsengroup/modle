#include "modle/dna.hpp"

#include <algorithm>
#include <boost/iostreams/filter/bzip2.hpp>
#include <filesystem>

#include "absl/container/flat_hash_set.h"
#include "fmt/printf.h"
#include "modle/extr_barrier.hpp"
#include "modle/lefs.hpp"
#include "modle/parsers.hpp"

namespace modle {

DNA::Bin::Bin(std::size_t idx, uint64_t start, uint64_t end,
              const std::vector<ExtrusionBarrier>& barriers)
    : _idx(idx),
      _start(static_cast<uint32_t>(start)),
      _end(static_cast<uint32_t>(end)),
      _extr_barriers(
          std::make_unique<std::vector<ExtrusionBarrier>>(barriers.begin(), barriers.end())) {
  assert(start <= UINT32_MAX);  // NOLINT;
  assert(end <= UINT32_MAX);    // NOLINT;
  assert(start <= end);         // NOLINT;
}

DNA::Bin::Bin(std::size_t idx, uint64_t start, uint64_t end)
    : _idx(idx),
      _start(static_cast<uint32_t>(start)),
      _end(static_cast<uint32_t>(end)),
      _extr_barriers(nullptr) {
  assert(start <= UINT32_MAX);  // NOLINT;
  assert(end <= UINT32_MAX);    // NOLINT;
  assert(start <= end);         // NOLINT;
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

ExtrusionBarrier* DNA::Bin::get_ptr_to_next_extr_barrier(ExtrusionBarrier* b, Direction d) const {
  assert(d != DNA::Direction::none);  // NOLINT;
  if (!this->_extr_barriers) {
    return nullptr;  // There are no barriers bound to this Bin
  }

  if (b == nullptr) {
    return &this->_extr_barriers->front();  // Return the first extr. barrier
  }
  assert(b >= &this->_extr_barriers->front() && b <= &this->_extr_barriers->back());  // NOLINT;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  for (++b; b < &this->_extr_barriers->back(); ++b) {
    if (d == DNA::Direction::both || b->get_direction_of_block() == d) {
      return b;  // Found a suitable extr. barrier
    }
  }
  return nullptr;  // Unable to find a suitable extr. barrier
}

ExtrusionBarrier* DNA::Bin::get_ptr_to_prev_extr_barrier(ExtrusionBarrier* b, Direction d) const {
  assert(d != DNA::Direction::none);  // NOLINT;
  if (!this->_extr_barriers) {
    return nullptr;  // There are no barriers bound to this Bin
  }

  if (b == nullptr) {
    return &this->_extr_barriers->back();  // Return the last extr. barrier
  }

  assert(b >= &this->_extr_barriers->front() && b <= &this->_extr_barriers->back());  // NOLINT;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  for (--b; b > &this->_extr_barriers->front(); --b) {
    if (d == DNA::Direction::both || b->get_direction_of_block() == d) {
      return b;  // Found a suitable extr. barrier
    }
  }
  return nullptr;  // Unable to find a suitable extr. barrier
}

std::vector<ExtrusionBarrier>& DNA::Bin::get_all_extr_barriers() const {
  return *this->_extr_barriers;
}

absl::InlinedVector<ExtrusionUnit*, 3>& DNA::Bin::get_extr_units() { return *this->_extr_units; }

ExtrusionBarrier* DNA::Bin::add_extr_barrier(ExtrusionBarrier b) {
  if (!this->_extr_barriers) {  // If this is the first ExtrusionBarrier, allocate the std::vector
    this->_extr_barriers =
        std::make_unique<std::vector<ExtrusionBarrier>>(std::vector<ExtrusionBarrier>{b});
  } else {
    this->_extr_barriers->emplace_back(b);
  }
  return &this->_extr_barriers->back();
}

void DNA::Bin::remove_extr_barrier(Direction d) {
  assert(d != Direction::none);  // NOLINT;
  if (!this->_extr_barriers) {
    throw std::logic_error(
        "Attempt to remove an extrusion barrier from a bin that doesn't have any!");
  }
  auto barrier = std::find_if(this->_extr_barriers->begin(), this->_extr_barriers->end(),
                              [&d](const ExtrusionBarrier& b) {
                                return d == DNA::Direction::both || b.get_direction_of_block() == d;
                              });
  if (barrier != this->_extr_barriers->end()) {
    this->_extr_barriers->erase(barrier);
    // Deallocate vector if there are no ExtrusionBarriers left
    if (this->_extr_barriers->empty()) {
      this->_extr_barriers = nullptr;
    }
    return;
  }
  throw std::logic_error(fmt::format("Bin doesn't have barriers blocking %s!",
                                     d == Direction::fwd ? "forward" : "reverse"));
}

uint64_t DNA::Bin::add_extr_unit_binding(ExtrusionUnit* const unit) {
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
  // TODO: Enable these checks only when compiling in Debug
  if (pos > this->length()) {
    throw std::logic_error(
        fmt::format("DNA::get_ptr_to_bin_from_pos(pos=%lu): pos > this->length(): %lu > %lu\n", pos,
                    pos, this->length()));
  }
  if ((pos / this->get_bin_size()) > this->get_n_bins()) {
    throw std::logic_error(
        fmt::format("(pos / this->bin_size()) >= this->nbins(): (%lu / %lu) >= %lu\n", pos,
                    this->get_bin_size(), this->get_n_bins()));
  }
  return &this->_bins[pos / this->get_bin_size()];
}

DNA::Bin& DNA::get_bin_from_pos(uint64_t pos) { return *this->get_ptr_to_bin_from_pos(pos); }

std::vector<DNA::Bin>::iterator DNA::begin() { return this->_bins.begin(); }
std::vector<DNA::Bin>::iterator DNA::end() { return this->_bins.end(); }
std::vector<DNA::Bin>::const_iterator DNA::cbegin() const { return this->_bins.cbegin(); }
std::vector<DNA::Bin>::const_iterator DNA::cend() const { return this->_bins.cend(); }

DNA::Bin* DNA::get_ptr_to_prev_bin(const DNA::Bin& current_bin) {
  assert(current_bin._idx > 0);  // NOLINT;
  return &this->_bins[current_bin._idx - 1];
}

DNA::Bin& DNA::get_prev_bin(const Bin& current_bin) {
  return *this->get_ptr_to_prev_bin(current_bin);
}

DNA::Bin& DNA::get_next_bin(const DNA::Bin& current_bin) {
  return *this->get_ptr_to_next_bin(current_bin);
}

DNA::Bin* DNA::get_ptr_to_next_bin(const DNA::Bin& current_bin) {
  assert(current_bin._idx < this->get_n_bins());  // NOLINT;
  return &this->_bins[current_bin._idx + 1];
}

DNA::Bin& DNA::get_first_bin() { return *this->get_ptr_to_first_bin(); }

DNA::Bin* DNA::get_ptr_to_first_bin() {
  assert(!this->_bins.empty());  // NOLINT;
  return &this->_bins.front();
}

DNA::Bin& DNA::get_last_bin() { return *this->get_ptr_to_last_bin(); }

DNA::Bin* DNA::get_ptr_to_last_bin() {
  assert(!this->_bins.empty());  // NOLINT;
  return &this->_bins.back();
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
                                             DNA::Direction direction) {
  if (!this->_extr_barriers) {
    this->_extr_barriers = std::make_unique<std::vector<ExtrusionBarrier>>(
        std::vector<ExtrusionBarrier>{{pos, prob_of_barrier_block, direction}});
  } else {
    this->_extr_barriers->emplace_back(pos, prob_of_barrier_block, direction);
  }
  return &this->_extr_barriers->back();
}

DNA::DNA(uint64_t length, uint32_t bin_size)
    : _bins(make_bins(length, bin_size)), _length(length), _bin_size(bin_size) {}

ExtrusionBarrier* DNA::add_extr_barrier(ExtrusionBarrier& b, uint32_t pos) {
  assert(b.get_prob_of_block() >= 0 && b.get_prob_of_block() <= 1);  // NOLINT;
  assert(pos <= this->length());                                     // NOLINT;
  assert(pos / this->get_bin_size() < this->get_n_bins());           // NOLINT;
  return this->_bins[pos / this->get_bin_size()].add_extr_barrier(b);
}

ExtrusionBarrier* DNA::add_extr_barrier(const BED& record) {
  assert(record.score >= 0 || record.score <= 1);                        // NOLINT;
  assert(record.chrom_end <= this->length());                            // NOLINT;
  assert(record.chrom_end / this->get_bin_size() < this->get_n_bins());  // NOLINT;
  assert(record.strand == '+' || record.strand == '-');                  // NOLINT;
  uint64_t pos = (record.chrom_end + record.chrom_start) / 2;
  DNA::Direction strand = record.strand == '+' ? DNA::Direction::fwd : DNA::Direction::rev;
  return this->_bins[pos / this->get_bin_size()].add_extr_barrier(pos, record.score, strand);
}

void DNA::remove_extr_barrier(uint32_t pos, Direction direction) {
  assert(pos <= this->length());                            // NOLINT;
  assert(pos / this->get_bin_size() < this->get_n_bins());  // NOLINT;
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
  if (bins.back()._end < length) {
    bins.emplace_back(DNA::Bin{static_cast<uint32_t>(bins.size()), start, length});
  }
  return bins;
}

uint64_t DNA::length() const { return this->_length; }

uint64_t DNA::get_n_bins() const { return this->_bins.size(); }

uint32_t DNA::get_bin_size() const { return this->_bin_size; }

uint64_t DNA::get_n_barriers() const {
  return std::accumulate(this->_bins.begin(), this->_bins.end(), 0UL,
                         [](uint64_t accumulator, const DNA::Bin& b) {
                           if (b._extr_barriers) {
                             accumulator += b._extr_barriers->size();
                           }
                           return accumulator;
                         });
}

Chromosome::Chromosome(std::string name, uint64_t length, uint32_t bin_size,
                       uint32_t diagonal_width)
    : name(std::move(name)),
      dna(length, bin_size),
      contacts(std::min(static_cast<uint64_t>(diagonal_width / bin_size), this->get_n_bins()),
               this->get_n_bins()) {}

uint64_t Chromosome::length() const { return this->dna.length(); }
uint64_t Chromosome::get_n_bins() const { return this->dna.get_n_bins(); }
uint32_t Chromosome::get_bin_size() const { return this->dna.get_bin_size(); }
uint64_t Chromosome::get_n_barriers() const { return this->barriers.size(); }

void Chromosome::sort_barriers_by_pos() {
  for (auto& bin : this->dna) {
    bin.sort_extr_barriers_by_pos();
    if (bin.has_extr_barrier()) {
      for (auto& b : bin.get_all_extr_barriers()) {
        this->barriers.push_back(&b);
      }
    }
  }
}

void Chromosome::write_contacts_to_tsv(std::string_view output_dir, bool force_overwrite) const {
  auto t0 = absl::Now();
  std::string path_to_outfile = std::filesystem::weakly_canonical(
      fmt::format("{}/{}_modle_cmatrix.tsv.bz2", output_dir, this->name));
  fmt::fprintf(stderr, "Writing contact matrix for '%s' to file '%s'...\n", this->name,
               path_to_outfile);
  if (!force_overwrite && std::filesystem::exists(path_to_outfile)) {
    fmt::fprintf(stderr, "File '%s' already exists, SKIPPING! Pass --force to overwrite...\n",
                 path_to_outfile);
    return;
  }

  const std::string header =
      fmt::format("#{}\t{}\t{}\t{}\t{}\n", this->name, this->get_bin_size(), 0, this->length(),
                  this->contacts.n_rows() * this->get_bin_size());

  auto [bytes_in, bytes_out] = this->contacts.write_to_tsv(path_to_outfile, header);
  fmt::fprintf(stderr,
               "DONE writing '%s' in %s! Compressed size: %.2f MB (compression ratio %.2fx)\n",
               // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
               path_to_outfile, absl::FormatDuration(absl::Now() - t0), bytes_out / 1.0e6,
               static_cast<double>(bytes_in) / bytes_out);
}

void Chromosome::write_barriers_to_tsv(std::string_view output_dir, bool force_overwrite) const {
  auto path_to_outfile = std::filesystem::weakly_canonical(
      fmt::format("{}/{}.extrusion_barriers.tsv.bz2", output_dir, this->name));
  if (!force_overwrite && std::filesystem::exists(path_to_outfile)) {
    fmt::fprintf(stderr, "File '%s' already exists. Pass --force to overwrite... SKIPPING.\n",
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
      buff = absl::StrCat(barrier->get_direction_of_block() == DNA::fwd ? "+" : "-",
                          barrier->get_pos(), "\n");
      out.write(buff.data(), buff.size());
      if (!out || !fp) {
        throw fmt::system_error(errno, "IO error while writing to file '{}'", path_to_outfile);
      }
    }
  } catch (const boost::iostreams::bzip2_error& err) {
    throw std::runtime_error(fmt::format("An error occurred while compressing file '{}': {}",
                                         path_to_outfile, err.what()));
  }
}

}  // namespace modle