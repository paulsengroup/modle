// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/genome.hpp"

#include <absl/container/btree_map.h>  // for btree_map
#include <absl/container/btree_set.h>  // for btree_set, btree_iterator
#include <absl/time/clock.h>           // for Now
#include <absl/time/time.h>            // for FormatDuration, operator-, Time
#include <fmt/format.h>                // for format, make_format_args, vformat_to
#include <fmt/std.h>
#include <spdlog/spdlog.h>  // for info

#include <algorithm>  // for max, max_element, find_if
#include <cassert>    // for assert
#include <cmath>      // for round
#include <exception>  // for exception
#include <filesystem>
#include <iosfwd>       // for streamsize
#include <memory>       // for shared_ptr, __shared_ptr_access
#include <numeric>      // for accumulate
#include <string>       // for string, char_traits
#include <string_view>  // for string_view, operator==, operator!=
#include <utility>      // for move, pair, pair<>::second_type
#include <vector>       // for vector

#include "modle/bed/bed.hpp"                  // for BED, Parser, BED_tree, BED_tree::at
#include "modle/chrom_sizes/chrom_sizes.hpp"  // for Parser
#include "modle/common/common.hpp"            // for bp_t, u32, u64, u8
#include "modle/common/fmt_helpers.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_PUSH, DISABLE_WAR...
#include "modle/common/utils.hpp"                       // for XXH3_Deleter, ndebug_defined, XXH...
#include "modle/contact_matrix_dense.hpp"               // for ContactMatrixDense
#include "modle/extrusion_barriers.hpp"                 // for ExtrusionBarrier

namespace modle {

namespace internal {
ContactMatrixLazy::ContactMatrixLazy(bp_t length, bp_t diagonal_width, bp_t bin_size) noexcept
    : _nrows((diagonal_width + bin_size - 1) / bin_size),
      _ncols((length + bin_size - 1) / bin_size) {}

ContactMatrixLazy::ContactMatrixLazy(ContactMatrix matrix) noexcept
    : _matrix(std::move(matrix)), _nrows(_matrix.nrows()), _ncols(_matrix.ncols()) {
  // Signal that _matrix has already been initialized
  std::call_once(this->_alloc_flag, []() {});
}

ContactMatrixLazy::ContactMatrixLazy(ContactMatrixLazy&& other) noexcept
    : _matrix(std::move(other._matrix)), _nrows(other._nrows), _ncols(other._ncols) {
  if (this->_matrix.npixels() == this->npixels()) {
    // Signal that _matrix has already been initialized
    std::call_once(this->_alloc_flag, []() {});
  }
  // This does not handle the case where a matrix was allocated, deallocated and then moved.
  // We don't really care about this case, so the current impl. should be fine
}

ContactMatrixLazy& ContactMatrixLazy::operator=(ContactMatrixLazy&& other) noexcept {
  if (this == &other) {
    return *this;
  }
  _matrix = std::move(other._matrix);
  if (this->_matrix.npixels() == this->npixels()) {  // See comments for move ctor
    std::call_once(this->_alloc_flag, []() {});
  }

  _nrows = other._nrows;
  _ncols = other._ncols;

  return *this;
}

auto ContactMatrixLazy::operator()() noexcept -> ContactMatrix& {
  std::call_once(this->_alloc_flag, [this]() {
    using MatrixT = decltype(this->_matrix);
    this->_matrix = MatrixT(this->_nrows, this->_ncols);
  });
  return this->_matrix;
}

auto ContactMatrixLazy::operator()() const noexcept -> const ContactMatrix& {
  std::call_once(this->_alloc_flag, [this]() {
    using MatrixT = decltype(this->_matrix);
    this->_matrix = MatrixT(this->_nrows, this->_ncols);
  });
  return this->_matrix;
}

void ContactMatrixLazy::deallocate() noexcept {
  std::call_once(this->_dealloc_flag, [this]() {
    using MatrixT = decltype(this->_matrix);
    MatrixT tmp{};
    std::swap(this->_matrix, tmp);
  });
}

Occupancy1DLazy::operator bool() const noexcept { return !this->_buff.empty(); }

Occupancy1DLazy::Occupancy1DLazy(bp_t length, bp_t bin_size) noexcept
    : _size((length + bin_size - 1) / bin_size) {}

Occupancy1DLazy::Occupancy1DLazy(std::vector<std::atomic<u64>> buff) noexcept
    : _buff(std::move(buff)), _size(_buff.size()) {
  // Signal _buff has already been initialized
  std::call_once(this->_alloc_flag, []() {});
}

Occupancy1DLazy::Occupancy1DLazy(Occupancy1DLazy&& other) noexcept
    : _buff(std::move(other._buff)), _size(other._size) {
  if (_buff.size() == this->size()) {
    // Signal _buff has already been initialized
    std::call_once(this->_alloc_flag, []() {});
  }
  // This does not handle the case where buffer was allocated, deallocated and then moved.
  // We don't really care about this case, so the current impl. should be fine
}

Occupancy1DLazy& Occupancy1DLazy::operator=(Occupancy1DLazy&& other) noexcept {
  if (this == &other) {
    return *this;
  }
  _buff = std::move(other._buff);
  if (_buff.size() == this->size()) {  // See comments for move ctor
    std::call_once(this->_alloc_flag, []() {});
  }
  _size = other._size;

  return *this;
}

const std::vector<std::atomic<u64>>& Occupancy1DLazy::operator()() const noexcept {
  std::call_once(this->_alloc_flag, [this]() {
    using BuffT = decltype(this->_buff);
    this->_buff = BuffT(_size);
    std::fill(this->_buff.begin(), this->_buff.end(), 0);
  });
  return this->_buff;
}

std::vector<std::atomic<u64>>& Occupancy1DLazy::operator()() noexcept {
  std::call_once(this->_alloc_flag, [this]() {
    using BuffT = decltype(this->_buff);
    this->_buff = BuffT(_size);
    std::fill(this->_buff.begin(), this->_buff.end(), 0);
  });
  return this->_buff;
}

void Occupancy1DLazy::deallocate() noexcept {
  std::call_once(this->_dealloc_flag, [this]() {
    using BuffT = decltype(this->_buff);
    BuffT tmp{};
    std::swap(this->_buff, tmp);
  });
}

}  // namespace internal

Chromosome::Chromosome(usize id, std::string name, bp_t size) noexcept
    : _name(std::move(name)), _id(id), _size(size) {}

std::string_view Chromosome::name() const noexcept { return this->_name; }
const char* Chromosome::name_cstr() const noexcept { return this->_name.c_str(); }

GenomicInterval::GenomicInterval(usize id, const std::shared_ptr<const Chromosome>& chrom,
                                 bp_t contact_matrix_resolution, bp_t diagonal_width)
    : GenomicInterval(id, chrom, 0, chrom->size(), contact_matrix_resolution, diagonal_width) {}

GenomicInterval::GenomicInterval(usize id, std::shared_ptr<const Chromosome> chrom, bp_t start,
                                 bp_t end, bp_t contact_matrix_resolution, bp_t diagonal_width)
    : GenomicInterval(id, std::move(chrom), start, end, contact_matrix_resolution, diagonal_width,
                      static_cast<const ExtrusionBarrier*>(nullptr),
                      static_cast<const ExtrusionBarrier*>(nullptr)) {}

const Chromosome& GenomicInterval::chrom() const noexcept {
  assert(!!this->_chrom);
  return *this->_chrom;
}

usize GenomicInterval::num_barriers() const { return this->_barriers.size(); }

auto GenomicInterval::barriers() const noexcept -> const std::vector<ExtrusionBarrier>& {
  return this->_barriers;
}
auto GenomicInterval::barriers() noexcept -> std::vector<ExtrusionBarrier>& {
  return this->_barriers;
}

auto GenomicInterval::contacts() const noexcept -> const ContactMatrix& {
  return this->_contacts();
}
auto GenomicInterval::contacts() noexcept -> ContactMatrix& { return this->_contacts(); }

auto GenomicInterval::lef_1d_occupancy() const noexcept -> const std::vector<std::atomic<u64>>& {
  return this->_lef_1d_occupancy();
}

auto GenomicInterval::lef_1d_occupancy() noexcept -> std::vector<std::atomic<u64>>& {
  return this->_lef_1d_occupancy();
}

void GenomicInterval::deallocate() noexcept { this->_contacts.deallocate(); }

struct ComputeBarrierStpResult {
  double stp_active;
  double stp_inactive;
};

[[nodiscard]] static constexpr ComputeBarrierStpResult compute_barrier_stp(
    double score, double default_stp_active, double default_stp_inactive) {
  if (score != 0.0) {
    // When the score field is zero (i.e. when the extr. barrier does not have a custom
    // occupancy), use the default self-transition probabilities
    const auto occupancy = score;
    const auto stp_active =
        ExtrusionBarrier::compute_stp_active_from_occupancy(default_stp_inactive, occupancy);
    return {stp_active, default_stp_inactive};
  }
  return {default_stp_active, default_stp_inactive};
}

void GenomicInterval::add_extrusion_barrier(const bed::BED& record,
                                            const double default_barrier_stp_active,
                                            const double default_barrier_stp_inactive) {
  assert(record.strand == '+' || record.strand == '-' || record.strand == '.');
  const auto pos = (record.chrom_start + record.chrom_end + 1) / 2;
  if (pos < this->start() || pos >= this->end()) {
    return;
  }

  const auto [barrier_stp_active, barrier_stp_inactive] =
      compute_barrier_stp(record.score, default_barrier_stp_active, default_barrier_stp_inactive);

  this->_barriers.emplace_back(pos, barrier_stp_active, barrier_stp_inactive, record.strand);
}

void GenomicInterval::add_extrusion_barriers(std::vector<ExtrusionBarrier> barriers) {
  if constexpr (utils::ndebug_not_defined()) {
    for ([[maybe_unused]] const auto& barrier : barriers) {
      assert(barrier.pos >= this->start());
      assert(barrier.pos < this->end());
    }
  }
  this->_barriers.insert(this->_barriers.end(), std::make_move_iterator(barriers.begin()),
                         std::make_move_iterator(barriers.end()));
}

Genome::Genome(const std::filesystem::path& path_to_chrom_sizes,
               const std::filesystem::path& path_to_extr_barriers,
               const std::filesystem::path& path_to_genomic_intervals,
               bp_t contact_matrix_resolution, bp_t contact_matrix_diagonal_witdh,
               const double default_barrier_pbb, const double default_barrier_puu,
               bool interpret_name_field_as_puu)
    : _chroms(import_chromosomes(path_to_chrom_sizes)),
      _intervals(import_genomic_intervals(path_to_genomic_intervals, _chroms,
                                          contact_matrix_resolution,
                                          contact_matrix_diagonal_witdh)),
      _size(std::accumulate(_chroms.begin(), _chroms.end(), usize(0),
                            [&](auto accumulator, const auto& chrom_ptr) {
                              return accumulator + chrom_ptr->size();
                            })),
      _simulated_size(std::accumulate(
          _intervals.begin(), _intervals.end(), usize(0),
          [&](auto accumulator, const GenomicInterval& gi) { return accumulator + gi.size(); })) {
  assert(!path_to_extr_barriers.empty());
  const auto t0 = absl::Now();
  spdlog::info(FMT_STRING("Importing extrusion barriers from {}..."), path_to_extr_barriers);
  // Parse all the records from the BED file. The parser will throw in case of duplicates.
  const auto barriers_gw =
      bed::Parser(path_to_extr_barriers, bed::BED::BED6).parse_all_in_interval_tree();
  _num_barriers = map_barriers_to_intervals(_intervals, barriers_gw, default_barrier_pbb,
                                            default_barrier_puu, interpret_name_field_as_puu);
  if (_num_barriers == 0) {
    spdlog::warn(FMT_STRING("Imported 0 barriers from {}. Is this intended?"),
                 path_to_extr_barriers);
  } else {
    spdlog::info(FMT_STRING("Imported {} barriers from {} in {}."), _num_barriers,
                 path_to_extr_barriers, absl::FormatDuration(absl::Now() - t0));
  }
}

std::vector<std::shared_ptr<const Chromosome>> Genome::import_chromosomes(
    const std::filesystem::path& path_to_chrom_sizes) {
  assert(!path_to_chrom_sizes.empty());

  const auto t0 = absl::Now();
  spdlog::info(FMT_STRING("Importing chromosomes from {}..."), path_to_chrom_sizes);

  absl::btree_map<std::string, usize> chrom_names;

  usize id = 0;
  std::vector<std::shared_ptr<const Chromosome>> buffer{};
  for (auto&& record : chrom_sizes::Parser(path_to_chrom_sizes).parse_all()) {
    if (auto it = chrom_names.find(record.chrom); it != chrom_names.end()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING(
              "Found duplicate entry for {} at line {} of file {}! First entry was at line {}"),
          record.chrom, id, path_to_chrom_sizes, it->second));
    }
    chrom_names.emplace(record.chrom, id);
    buffer.emplace_back(
        std::make_shared<const Chromosome>(id++, std::move(record.chrom), record.size()));
  }

  if (buffer.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to import any chromosome from {}!"), path_to_chrom_sizes));
  }

  spdlog::info(FMT_STRING("Imported {} chromosomes in {}."), buffer.size(),
               absl::FormatDuration(absl::Now() - t0));
  return buffer;
}

absl::btree_set<GenomicInterval> Genome::import_genomic_intervals(
    const std::filesystem::path& path_to_bed,
    const std::vector<std::shared_ptr<const Chromosome>>& chromosomes,
    bp_t contact_matrix_resolution, bp_t diagonal_width) {
  assert(!chromosomes.empty());
  absl::btree_set<GenomicInterval> buffer{};
  if (path_to_bed.empty()) {
    spdlog::debug(
        FMT_STRING("Path to genomic regions to simulate is empty. Assuming whole chromosomes are "
                   "to be simulated!"));
    std::transform(chromosomes.begin(), chromosomes.end(), std::inserter(buffer, buffer.begin()),
                   [&](const auto& chrom_ptr) {
                     return GenomicInterval{chrom_ptr->id(), chrom_ptr, contact_matrix_resolution,
                                            diagonal_width};
                   });
    return buffer;
  }

  const auto t0 = absl::Now();
  spdlog::info(FMT_STRING("Importing genomic intervals from {}..."), path_to_bed);

  absl::btree_map<std::string_view, std::shared_ptr<const Chromosome>> chrom_names;
  std::transform(
      chromosomes.begin(), chromosomes.end(), std::inserter(chrom_names, chrom_names.begin()),
      [](auto chrom_ptr) { return std::make_pair(chrom_ptr->name(), std::move(chrom_ptr)); });

  usize id = 0;
  // Parse all the records from the BED file. The parser will throw in case of duplicates.
  for (auto&& [chrom_name, intervals] :
       bed::Parser(path_to_bed, bed::BED::BED3).parse_all_in_interval_tree()) {
    auto it = chrom_names.find(chrom_name);
    if (it == chrom_names.end()) {
      spdlog::warn(
          FMT_STRING(
              "Skipping {} intervals from {}. Reason: {} was not present in the .chrom.sizes file"),
          intervals.size(), chrom_name, chrom_name);
      continue;
    }

    auto chrom_ptr = it->second;
    for (const auto& interval : intervals.data()) {
      buffer.emplace(GenomicInterval{id++, chrom_ptr, interval.chrom_start, interval.chrom_end,
                                     contact_matrix_resolution, diagonal_width});
    }
  }

  if (buffer.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to import any interval from {}!"), path_to_bed));
  }

  spdlog::info(FMT_STRING("Imported {} intervals in {}."), buffer.size(),
               absl::FormatDuration(absl::Now() - t0));

  return buffer;
}

[[nodiscard]] static std::vector<ExtrusionBarrier> generate_barriers_from_bed_records(
    const absl::Span<const bed::BED> records, double default_barrier_pbb,
    double default_barrier_puu, bool interpret_name_field_as_puu) {
  std::vector<ExtrusionBarrier> buff;
  buff.reserve(records.size());

  for (const auto& record : records) {
    try {
      if (record.strand == '.') {
        continue;
      }
      if (record.strand != '-' && record.strand != '+') {
        throw std::runtime_error(fmt::format(FMT_STRING("invalid strand '{}'"), record.strand));
      }

      if (record.score < 0 || record.score > 1) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("invalid score field: expected a score between 0 and 1, found {:.4g}."),
            record.score));
      }

      auto puu = default_barrier_puu;
      try {
        if (interpret_name_field_as_puu) {
          utils::parse_numeric_or_throw(record.name, puu);
          if (puu < 0 || puu > 1) {
            throw std::runtime_error("");
          }
        }
      } catch (const std::exception&) {
        throw std::runtime_error(fmt::format(
            FMT_STRING(
                "invalid name field: expected name to be a number between 0 and 1, found {}."),
            record.name));
      }

      const auto pos = (record.chrom_start + record.chrom_end + 1) / 2;
      const auto [barrier_stp_active, barrier_stp_inactive] =
          compute_barrier_stp(record.score, default_barrier_pbb, default_barrier_puu);

      buff.emplace_back(pos, barrier_stp_active, barrier_stp_inactive, record.strand);
    } catch (const std::exception& e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Found invalid extrusion barrier {:bed3}: {}"), record, e.what()));
    }
  }

  return buff;
}

usize Genome::map_barriers_to_intervals(absl::btree_set<GenomicInterval>& intervals,
                                        const bed::BED_tree<>& barriers_bed,
                                        double default_barrier_pbb, double default_barrier_puu,
                                        bool interpret_name_field_as_puu) {
  assert(!intervals.empty());
  usize tot_num_barriers = 0;
  for (auto& interval : intervals) {
    auto barriers = generate_barriers_from_bed_records(
        barriers_bed.find_overlaps(std::string{interval.chrom().name()},
                                   utils::conditional_static_cast<u64>(interval.start()),
                                   utils::conditional_static_cast<u64>(interval.end())),
        default_barrier_pbb, default_barrier_puu, interpret_name_field_as_puu);

    tot_num_barriers += barriers.size();
    interval.add_extrusion_barriers(std::move(barriers));
  }
  return tot_num_barriers;
}

auto Genome::begin() -> iterator { return this->_intervals.begin(); }
auto Genome::end() -> iterator { return this->_intervals.end(); }

auto Genome::begin() const -> const_iterator { return this->_intervals.cbegin(); }
auto Genome::end() const -> const_iterator { return this->_intervals.cend(); }

auto Genome::cbegin() const -> const_iterator { return this->_intervals.cbegin(); }
auto Genome::cend() const -> const_iterator { return this->_intervals.cend(); }

auto Genome::find(const GenomicInterval& query) -> iterator { return this->_intervals.find(query); }
auto Genome::find(const GenomicInterval& query) const -> const_iterator {
  return this->_intervals.find(query);
}

auto Genome::find(usize query) -> iterator {
  return std::find_if(this->_intervals.begin(), this->_intervals.end(),
                      [&](const auto& gi) { return gi.id() == query; });
}
auto Genome::find(usize query) const -> const_iterator {
  return std::find_if(this->_intervals.begin(), this->_intervals.end(),
                      [&](const auto& gi) { return gi.id() == query; });
}

auto Genome::find(const Chromosome& query) -> iterator {
  return std::find_if(this->_intervals.begin(), this->_intervals.end(),
                      [&](const auto& gi) { return gi.chrom() == query; });
}
auto Genome::find(const Chromosome& query) const -> const_iterator {
  return std::find_if(this->_intervals.begin(), this->_intervals.end(),
                      [&](const auto& gi) { return gi.chrom() == query; });
}

auto Genome::find(std::string_view query) -> iterator {
  return std::find_if(this->_intervals.begin(), this->_intervals.end(),
                      [&](const auto& gi) { return gi.chrom().name() == query; });
}
auto Genome::find(std::string_view query) const -> const_iterator {
  return std::find_if(this->_intervals.begin(), this->_intervals.end(),
                      [&](const auto& gi) { return gi.chrom().name() == query; });
}

bool Genome::contains(const GenomicInterval& query) const {
  return this->_intervals.contains(query);
}
bool Genome::contains(usize query) const { return query < this->_chroms.size(); }
bool Genome::contains(const Chromosome& query) const { return this->contains(query.id()); }
bool Genome::contains(std::string_view query) const {
  return std::find_if(this->_chroms.begin(), this->_chroms.end(), [&](const auto& chrom_ptr) {
           return chrom_ptr->name() == query;
         }) != this->_chroms.end();
}

usize Genome::num_chromosomes() const noexcept {
  return utils::conditional_static_cast<usize>(this->_chroms.size());
}

const Chromosome& Genome::chromosome_with_longest_name() const noexcept {
  assert(!this->_chroms.empty());
  return **std::max_element(this->_chroms.begin(), this->_chroms.end(),
                            [&](const auto& chrom_ptr1, const auto& chrom_ptr2) {
                              return chrom_ptr1->name().size() < chrom_ptr2->name().size();
                            });
}

const Chromosome& Genome::longest_chromosome() const noexcept {
  assert(!this->_chroms.empty());
  return **std::max_element(this->_chroms.begin(), this->_chroms.end(),
                            [&](const auto& chrom_ptr1, const auto& chrom_ptr2) {
                              return chrom_ptr1->size() < chrom_ptr2->size();
                            });
}

const GenomicInterval& Genome::longest_interval() const noexcept {
  assert(!this->_intervals.empty());
  return *std::max_element(
      this->_intervals.begin(), this->_intervals.end(),
      [&](const auto& gi1, const auto& gi2) { return gi1.size() < gi2.size(); });
}

const GenomicInterval& Genome::interval_with_most_barriers() const noexcept {
  assert(!this->_intervals.empty());
  return *std::max_element(this->_intervals.begin(), this->_intervals.end(),
                           [&](const auto& gi1, const auto& gi2) {
                             return gi1.barriers().size() < gi2.barriers().size();
                           });
}

usize Genome::max_target_contacts(usize bin_size, usize diagonal_width,
                                  double target_contact_density, usize simulation_iterations,
                                  double lef_fraction_contact_sampling, double nlefs_per_mbp,
                                  usize ncells) const {
  const auto& interval = this->longest_interval();
  if (target_contact_density == 0.0) {
    const auto nlefs = static_cast<usize>(
        std::round(nlefs_per_mbp * (static_cast<double>(interval.size()) / 1.0e6)));
    return static_cast<usize>(
        (static_cast<double>(simulation_iterations * nlefs) * lef_fraction_contact_sampling) /
        static_cast<double>(ncells));
  }

  const auto npixels =
      ((interval.size() + bin_size - 1) / bin_size) * ((diagonal_width + bin_size - 1) / bin_size);

  return static_cast<usize>(std::round((static_cast<double>(npixels) * target_contact_density) /
                                       static_cast<double>(ncells)));
}

}  // namespace modle
