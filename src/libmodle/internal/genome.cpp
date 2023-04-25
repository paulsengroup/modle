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

Chromosome::Chromosome(usize id, const bed::BED& chrom,
                       const IITree<bp_t, ExtrusionBarrier>& barriers)
    : Chromosome(id, chrom.chrom, chrom.thick_start, chrom.thick_end, chrom.size(), barriers) {}

Chromosome::Chromosome(usize id, const bed::BED& chrom)
    : _name(chrom.chrom),
      _start(chrom.chrom_start),
      _end(chrom.chrom_end),
      _size(chrom.size()),
      _id(id) {}

Chromosome::Chromosome(usize id, const bed::BED& chrom, IITree<bp_t, ExtrusionBarrier>&& barriers)
    : Chromosome(id, chrom.chrom, chrom.thick_start, chrom.thick_end, chrom.size(), barriers) {}

DISABLE_WARNING_PUSH
DISABLE_WARNING_SIGN_CONVERSION
DISABLE_WARNING_CONVERSION
Chromosome::Chromosome(usize id, std::string_view chrom_name, bp_t chrom_start, bp_t chrom_end,
                       bp_t chrom_size)
    : _name(chrom_name), _start(chrom_start), _end(chrom_end), _size(chrom_size), _id(id) {
  assert(chrom_start <= chrom_end);
  assert(chrom_end - chrom_start <= chrom_size);
}

Chromosome::Chromosome(usize id, std::string_view chrom_name, bp_t chrom_start, bp_t chrom_end,
                       bp_t chrom_size, const IITree<bp_t, ExtrusionBarrier>& barriers)
    : _name(chrom_name),
      _start(chrom_start),
      _end(chrom_end),
      _size(chrom_size),
      _id(id),
      _barriers(barriers) {
  assert(chrom_start <= chrom_end);
  assert(chrom_end - chrom_start <= chrom_size);

  _barriers.make_BST();
}

Chromosome::Chromosome(usize id, std::string_view chrom_name, bp_t chrom_start, bp_t chrom_end,
                       bp_t chrom_size, IITree<bp_t, ExtrusionBarrier>&& barriers)
    : _name(chrom_name),
      _start(chrom_start),
      _end(chrom_end),
      _size(chrom_size),
      _id(id),
      _barriers(std::move(barriers)) {
  assert(chrom_start <= chrom_end);
  assert(chrom_end - chrom_start <= chrom_size);

  _barriers.make_BST();
}
DISABLE_WARNING_POP

Chromosome::Chromosome(const Chromosome& other)
    : _name(other._name),
      _start(other._start),
      _end(other._end),
      _size(other._size),
      _id(other._id),
      _barriers(other._barriers),
      _contacts(other._contacts) {
  _barriers.make_BST();
}

Chromosome::Chromosome(Chromosome&& other) noexcept
    : _name(std::move(other._name)),
      _start(other._start),
      _end(other._end),
      _size(other._size),
      _id(other._id),
      _barriers(std::move(other._barriers)),
      _contacts(std::move(other._contacts)) {
  _barriers.make_BST();
}

Chromosome& Chromosome::operator=(const Chromosome& other) {
  if (this == &other) {
    return *this;
  }

  _id = other._id;
  _name = other._name;
  _size = other._size;
  _start = other._start;
  _end = other._end;
  _barriers = other._barriers;
  _contacts = other._contacts;

  _barriers.make_BST();

  return *this;
}

Chromosome& Chromosome::operator=(Chromosome&& other) noexcept {
  if (this == &other) {
    return *this;
  }

  _id = other._id;
  _name = std::move(other._name);
  _size = other._size;
  _start = other._start;
  _end = other._end;
  _barriers = std::move(other._barriers);
  _contacts = std::move(other._contacts);

  _barriers.make_BST();

  return *this;
}

bool Chromosome::operator==(const Chromosome& other) const noexcept(utils::ndebug_defined()) {
  if (this->_id == (std::numeric_limits<usize>::max)()) {
    return this->name() == other.name() && this->size() == other.size();
  }
  return this->_id == other._id;
}

bool Chromosome::operator==(std::string_view other_name) const noexcept {
  return this->name() == other_name;
}

bool Chromosome::operator<(const Chromosome& other) const noexcept(utils::ndebug_defined()) {
  if (this->_id != (std::numeric_limits<usize>::max)() &&
      other._id != (std::numeric_limits<usize>::max)()) {
    return this->_id < other._id;
  }

  if (this->name() != other.name()) {
    return this->name() < other.name();
  }
  return this->size() < other.size();
}

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

void Chromosome::add_extrusion_barrier(const bed::BED& record,
                                       const double default_barrier_stp_active,
                                       const double default_barrier_stp_inactive) {
  assert(record.strand == '+' || record.strand == '-' || record.strand == '.');
  const auto pos = (record.chrom_start + record.chrom_end + 1) / 2;
  if (pos < this->start_pos() || pos >= this->end_pos()) {
    // Barrier lies outside of the genomic regions to be simulated
    return;
  }

  const auto [barrier_stp_active, barrier_stp_inactive] =
      compute_barrier_stp(record.score, default_barrier_stp_active, default_barrier_stp_inactive);

  this->_barriers.emplace(
      record.chrom_start, record.chrom_end,
      ExtrusionBarrier{pos, barrier_stp_active, barrier_stp_inactive, record.strand});
}

usize Chromosome::id() const { return this->_id; }
std::string_view Chromosome::name() const { return this->_name; }
const char* Chromosome::name_cstr() const { return this->_name.c_str(); }

usize Chromosome::num_lefs(double nlefs_per_mbp) const {
  return static_cast<usize>((static_cast<double>(this->simulated_size()) / 1.0e6) * nlefs_per_mbp);
}
usize Chromosome::num_barriers() const { return this->_barriers.size(); }

const IITree<bp_t, ExtrusionBarrier>& Chromosome::barriers() const { return this->_barriers; }
IITree<bp_t, ExtrusionBarrier>& Chromosome::barriers() { return this->_barriers; }

bool Chromosome::allocate_contact_matrix(bp_t bin_size, bp_t diagonal_width) {
  if (std::scoped_lock lck(this->_buff_mtx); this->_contacts.npixels() == 0) {
    this->_contacts.unsafe_resize(this->simulated_size(), diagonal_width, bin_size);
    return true;
  }
  return false;
}

bool Chromosome::allocate_lef_occupancy_buffer(modle::bp_t bin_size) {
  if (std::scoped_lock lck(this->_buff_mtx); this->_lef_1d_occupancy.empty()) {
    // Can't resize a vector of atomics
    using BuffT = decltype(this->_lef_1d_occupancy);
    this->_lef_1d_occupancy = BuffT((this->size() + bin_size - 1) / bin_size);
    std::fill(this->_lef_1d_occupancy.begin(), this->_lef_1d_occupancy.end(), 0);
    return true;
  }
  return false;
}

bool Chromosome::deallocate_contact_matrix() {
  if (std::scoped_lock lck(this->_buff_mtx); this->_contacts.npixels() != 0) {
    using MatrixT = decltype(this->_contacts);
    MatrixT tmp{};
    std::swap(this->_contacts, tmp);
    return true;
  }
  return false;
}

bool Chromosome::deallocate_lef_occupancy_buffer() {
  if (std::scoped_lock lck(this->_buff_mtx); !this->_lef_1d_occupancy.empty()) {
    using BuffT = decltype(this->_lef_1d_occupancy);
    BuffT tmp{};
    std::swap(this->_lef_1d_occupancy, tmp);
    return true;
  }
  return false;
}

usize Chromosome::npixels() const { return this->contacts().npixels(); }

const Chromosome::contact_matrix_t& Chromosome::contacts() const noexcept {
  return this->_contacts;
}

Chromosome::contact_matrix_t& Chromosome::contacts() noexcept { return this->_contacts; }

const std::vector<std::atomic<u64>>& Chromosome::lef_1d_occupancy() const noexcept {
  return this->_lef_1d_occupancy;
}

std::vector<std::atomic<u64>>& Chromosome::lef_1d_occupancy() noexcept {
  return this->_lef_1d_occupancy;
}

u64 Chromosome::hash(XXH3_state_t* const xxh_state, u64 seed, usize cell_id) const {
  auto handle_errors = [&](const auto& status) {
    if (MODLE_UNLIKELY(status == XXH_ERROR || !xxh_state)) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Failed to hash \"{}\" for cell #{} using seed {}"), this->name(),
                      cell_id, seed));
    }
  };

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USED_BUT_MARKED_UNUSED
  handle_errors(XXH3_64bits_reset_withSeed(xxh_state, seed));
  handle_errors(
      XXH3_64bits_update(xxh_state, this->_name.data(), this->_name.size() * sizeof(char)));
  handle_errors(XXH3_64bits_update(xxh_state, &this->_size, sizeof(decltype(this->_size))));
  handle_errors(XXH3_64bits_update(xxh_state, &cell_id, sizeof(decltype(cell_id))));

  return utils::conditional_static_cast<u64>(XXH3_64bits_digest(xxh_state));
  DISABLE_WARNING_POP
}

u64 Chromosome::hash(u64 seed, usize cell_id) const {
  std::unique_ptr<XXH3_state_t, utils::XXH3_Deleter> xxh_state{XXH3_createState()};
  return this->hash(xxh_state.get(), seed, cell_id);
}

Genome::Genome(const std::filesystem::path& path_to_chrom_sizes,
               const std::filesystem::path& path_to_extr_barriers,
               const std::filesystem::path& path_to_chrom_subranges,
               const double default_barrier_pbb, const double default_barrier_puu,
               bool interpret_name_field_as_puu)
    : _chromosomes(instantiate_genome(path_to_chrom_sizes, path_to_extr_barriers,
                                      path_to_chrom_subranges, default_barrier_pbb,
                                      default_barrier_puu, interpret_name_field_as_puu)) {}

absl::btree_set<Chromosome> Genome::import_chromosomes(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_chrom_subranges) {
  assert(!path_to_chrom_sizes.empty());

  const auto t0 = absl::Now();
  if (path_to_chrom_subranges.empty()) {
    spdlog::info(FMT_STRING("Importing chromosomes from file {}..."), path_to_chrom_sizes);
  } else {
    spdlog::info(FMT_STRING("Importing chromosomes from files {} and {}..."), path_to_chrom_sizes,
                 path_to_chrom_subranges);
  }
  auto print_status_update_on_return = [&](auto num_chromosomes) {
    spdlog::info(FMT_STRING("Imported {} chromosomes in {}."), num_chromosomes,
                 absl::FormatDuration(absl::Now() - t0));
  };

  // Parse chrom. sizes and build the set of chromosome to be simulated.
  // When the BED file with the chrom. subranges is not provided, all the chromosome in the
  // chrom.sizes file will be selected and returned. When a BED file with the chrom. subranges is
  // available, then only chromosomes that are present in both files will be selected. Furthermore
  // we are also checking that the subrange lies within the genomic coordinates specified in the
  // chrom. sizes file
  usize id = 0;
  absl::btree_set<Chromosome> chromosomes;
  if (path_to_chrom_subranges.empty()) {
    for (auto record : chrom_sizes::Parser(path_to_chrom_sizes).parse_all()) {
      record.thick_start = record.chrom_start;
      record.thick_end = record.chrom_end;
      chromosomes.emplace(id++, std::move(record));
    }
    print_status_update_on_return(chromosomes.size());
    return chromosomes;
  }

  // Parse chrom subranges from BED. We parse everything at once to deal with duplicate entries
  const auto chrom_ranges = [&]() {
    absl::btree_map<std::string, bed::BED> ranges;
    for (const auto& record : bed::Parser(path_to_chrom_subranges, bed::BED::BED3).parse_all()) {
      // passing std::string{record.chrom} is required in order to make GCC 8.4 happy
      const auto new_insertion = ranges.try_emplace(std::string{record.chrom}, record).second;
      if (!new_insertion) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Found more than one entry for chromosome \"{}\" in file {}"),
                        record.chrom, path_to_chrom_subranges));
      }
    }
    return ranges;
  }();

  for (auto record : chrom_sizes::Parser(path_to_chrom_sizes).parse_all()) {
    const auto match = chrom_ranges.find(record.chrom);
    if (match != chrom_ranges.end()) {
      const auto& begin_pos = match->second.chrom_start;
      const auto& end_pos = match->second.chrom_end;
      if (begin_pos < record.chrom_start || end_pos > record.chrom_end) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("According to the chrom.sizes file {}, chromosome \"{}\" should have a size "
                       "of \"{}\", but the subrange specified in BED file {} extends past this "
                       "region: range {}:{}-{} does not fit in range {}:{}-{}"),
            path_to_chrom_sizes, record.chrom, record.chrom_end, path_to_chrom_subranges,
            record.chrom, begin_pos, end_pos, record.chrom, record.chrom_start, record.chrom_end));
      }

      record.thick_start = begin_pos;
      record.thick_end = end_pos;
    } else {
      record.thick_start = record.chrom_start;
      record.thick_end = record.chrom_end;
    }

    chromosomes.emplace(id++, record.chrom, record.thick_start, record.thick_end, record.size(),
                        match != chrom_ranges.end());
  }

  print_status_update_on_return(chromosomes.size());
  return chromosomes;
}

usize Genome::import_barriers(absl::btree_set<Chromosome>& chromosomes,
                              const std::filesystem::path& path_to_extr_barriers,
                              double default_barrier_pbb, double default_barrier_puu,
                              bool interpret_name_field_as_puu) {
  assert(!chromosomes.empty());
  assert(!path_to_extr_barriers.empty());

  const auto t0 = absl::Now();
  spdlog::info(FMT_STRING("Importing extrusion barriers from file {}..."), path_to_extr_barriers);

  usize tot_num_barriers = 0;

  // Parse all the records from the BED file. The parser will throw in case of duplicates.
  const auto barriers =
      bed::Parser(path_to_extr_barriers, bed::BED::BED6).parse_all_in_interval_tree();

  for (auto& chrom : chromosomes) {
    if (const auto chrom_name = std::string{chrom.name()}; barriers.contains(chrom_name)) {
      for (const auto& record : barriers.at(chrom_name).data()) {
        if (record.score < 0 || record.score > 1) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("Invalid score field detected for record {}[{}-{}]: expected a score "
                         "between 0 and 1, got {:.4g}."),
              record.chrom, record.chrom_start, record.chrom_end, record.score));
        }

        if (interpret_name_field_as_puu) {
          const auto puu = utils::parse_numeric_or_throw<double>(record.name);
          if (puu < 0 || puu > 1) {
            throw std::runtime_error(fmt::format(
                FMT_STRING(
                    "Invalid score field detected for record {}[{}-{}]: expected name field to be "
                    "between 0 and 1, got {:.4g}."),
                record.chrom, record.chrom_start, record.chrom_end, puu));
          }
          chrom.add_extrusion_barrier(record, default_barrier_pbb, puu);
        } else {
          chrom.add_extrusion_barrier(record, default_barrier_pbb, default_barrier_puu);
        }
        ++tot_num_barriers;
      }
    }
    chrom.barriers().make_BST();
  }
  spdlog::info(FMT_STRING("Imported {} barriers in {}."), tot_num_barriers,
               absl::FormatDuration(absl::Now() - t0));
  return tot_num_barriers;
}

absl::btree_set<Chromosome> Genome::instantiate_genome(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_extr_barriers,
    const std::filesystem::path& path_to_chrom_subranges, double default_barrier_pbb,
    double default_barrier_puu, bool interpret_name_field_as_puu) {
  auto chroms = import_chromosomes(path_to_chrom_sizes, path_to_chrom_subranges);

  if (chroms.empty()) {
    if (path_to_chrom_subranges.empty()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Unable to import any chromosome from file {}. Is the file empty?"),
          path_to_chrom_sizes));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to import any chromosome from file {} using the ranges "
                               "specified by file {}. Please make sure that neither of the file is "
                               "empty and that both files refer to the same genome (i.e. they use "
                               "the same genomic coordinates)."),
                    path_to_chrom_sizes, path_to_chrom_subranges));
  }

  const auto tot_barriers_imported =
      import_barriers(chroms, path_to_extr_barriers, default_barrier_pbb, default_barrier_puu,
                      interpret_name_field_as_puu);
  if (tot_barriers_imported == 0) {
    spdlog::warn(
        FMT_STRING(
            "Unable to import any barrier from file {}. Please make sure this was not a mistake."),
        path_to_extr_barriers);
  }

  return chroms;
}

Genome::iterator Genome::begin() { return this->_chromosomes.begin(); }
Genome::iterator Genome::end() { return this->_chromosomes.end(); }

Genome::const_iterator Genome::begin() const { return this->_chromosomes.cbegin(); }
Genome::const_iterator Genome::end() const { return this->_chromosomes.cend(); }

Genome::const_iterator Genome::cbegin() const { return this->_chromosomes.cbegin(); }
Genome::const_iterator Genome::cend() const { return this->_chromosomes.cend(); }

Genome::iterator Genome::find(const Chromosome& other_chrom) {
  return this->_chromosomes.find(other_chrom);
}

Genome::const_iterator Genome::find(const Chromosome& other_chrom) const {
  return this->_chromosomes.find(other_chrom);
}

Genome::iterator Genome::find(std::string_view other_chrom_name) {
  return std::find_if(this->_chromosomes.begin(), this->_chromosomes.end(),
                      [&](const auto& chrom) { return chrom.name() == other_chrom_name; });
}

Genome::const_iterator Genome::find(std::string_view other_chrom_name) const {
  return std::find_if(this->_chromosomes.begin(), this->_chromosomes.end(),
                      [&](const auto& chrom) { return chrom.name() == other_chrom_name; });
}

bool Genome::contains(const Chromosome& other_chromosome) const {
  return this->_chromosomes.contains(other_chromosome);
}

bool Genome::contains(std::string_view other_chrom_name) const {
  return this->find(other_chrom_name) != this->end();
}

usize Genome::size() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
      [](auto accumulator, const auto& chrom) { return accumulator + chrom.size(); });
}

usize Genome::number_of_chromosomes() const {
  return utils::conditional_static_cast<usize>(this->_chromosomes.size());
}

usize Genome::simulated_size() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), usize(0),
                         [](auto accumulator, const auto& chrom) {
                           return accumulator + (chrom.end_pos() - chrom.start_pos());
                         });
}

const Chromosome& Genome::chromosome_with_longest_name() const {
  assert(!this->_chromosomes.empty());
  return *std::max_element(
      this->_chromosomes.begin(), this->_chromosomes.end(),
      [&](const auto& c1, const auto& c2) { return c1.name().size() < c2.name().size(); });
}

const Chromosome& Genome::longest_chromosome() const {
  assert(!this->_chromosomes.empty());
  return *std::max_element(
      this->_chromosomes.begin(), this->_chromosomes.end(),
      [&](const auto& c1, const auto& c2) { return c1.simulated_size() < c2.simulated_size(); });
}

const Chromosome& Genome::chromosome_with_max_nbarriers() const {
  assert(!this->_chromosomes.empty());
  return *std::max_element(
      this->_chromosomes.begin(), this->_chromosomes.end(),
      [&](const auto& c1, const auto& c2) { return c1.num_barriers() < c2.num_barriers(); });
}

usize Genome::max_target_contacts(usize bin_size, usize diagonal_width,
                                  double target_contact_density, usize simulation_iterations,
                                  double lef_fraction_contact_sampling, double nlefs_per_mbp,
                                  usize ncells) const {
  const auto& chrom = this->longest_chromosome();
  if (target_contact_density == 0.0) {
    const auto nlefs = static_cast<usize>(
        std::round(nlefs_per_mbp * (static_cast<double>(chrom.simulated_size()) / 1.0e6)));
    return static_cast<usize>(
        (static_cast<double>(simulation_iterations * nlefs) * lef_fraction_contact_sampling) /
        static_cast<double>(ncells));
  }

  const auto npixels = ((chrom.simulated_size() + bin_size - 1) / bin_size) *
                       ((diagonal_width + bin_size - 1) / bin_size);

  return static_cast<usize>(std::round((static_cast<double>(npixels) * target_contact_density) /
                                       static_cast<double>(ncells)));
}

}  // namespace modle
