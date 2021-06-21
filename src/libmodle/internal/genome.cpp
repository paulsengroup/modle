#include "modle/genome.hpp"

#include <absl/container/btree_set.h>  // for btree_set
#include <fmt/format.h>                // for format, FMT_STRING
#include <fmt/ostream.h>

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

#include "modle/bed.hpp"                                // for BED
#include "modle/chrom_sizes.hpp"                        // for ChromSize
#include "modle/common/common.hpp"                      // for bp_t, contacts_t
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNING_PUSH
#include "modle/common/utils.hpp"                       // for ndebug_defined
#include "modle/contacts.hpp"                           // for ContactMatrix
#include "modle/extrusion_barriers.hpp"                 // for ExtrusionBarrier

namespace modle {

Chromosome::Chromosome(size_t id, const bed::BED& chrom, absl::Span<const bed::BED> barriers)
    : Chromosome(id, chrom.chrom, chrom.thick_start, chrom.thick_end, chrom.size(),
                 barriers.begin(), barriers.end()) {}

Chromosome::Chromosome(size_t id, const bed::BED& chrom, const interval_tree_value_t& barriers)
    : Chromosome(id, chrom.chrom, chrom.thick_start, chrom.thick_end, chrom.size(), barriers) {}

Chromosome::Chromosome(size_t id, const bed::BED& chrom, interval_tree_value_t&& barriers)
    : Chromosome(id, chrom.chrom, chrom.thick_start, chrom.thick_end, chrom.size(), barriers) {}

Chromosome::Chromosome(const Chromosome& other)
    : _name(other._name),
      _start(other._start),
      _end(other._end),
      _size(other._size),
      _id(other._id),
      _barriers(other._barriers),
      _contacts(std::make_unique<contact_matrix_t>(*other._contacts)) {
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
  _contacts = std::make_unique<contact_matrix_t>(*other._contacts);

  _barriers.make_BST();

  return *this;
}

bool Chromosome::operator==(const Chromosome& other) const noexcept(utils::ndebug_defined()) {
  if (this->_id == std::numeric_limits<size_t>::max()) {
    return this->name() == other.name() && this->size() == other.size();
  }
  return this->_id == other._id;
}

bool Chromosome::operator==(std::string_view other_name) const noexcept {
  return this->name() == other_name;
}

bool Chromosome::operator<(const Chromosome& other) const noexcept(utils::ndebug_defined()) {
  if (this->_id != std::numeric_limits<size_t>::max() &&
      other._id != std::numeric_limits<size_t>::max()) {
    return this->_id < other._id;
  }

  if (this->name() != other.name()) {
    return this->name() < other.name();
  }
  return this->size() < other.size();
}

void Chromosome::add_extrusion_barrier(const bed::BED& barrier) {
  this->_barriers.insert(barrier.chrom_start, barrier.chrom_end, barrier);
}

void Chromosome::add_extrusion_barrier(const absl::Span<const bed::BED> barriers) {
  this->add_extrusion_barrier(barriers.begin(), barriers.end());
}

void Chromosome::add_extrusion_barrier(bed::BED&& barrier) {
  this->_barriers.emplace(barrier.chrom_start, barrier.chrom_end, std::move(barrier));
}

size_t Chromosome::id() const { return this->_id; }
std::string_view Chromosome::name() const { return this->_name; }
const char* Chromosome::name_cstr() const { return this->_name.c_str(); }

bool Chromosome::ok() const { return !this->_barriers.empty(); }

size_t Chromosome::nlefs(double nlefs_per_mbp) const {  // NOLINTNEXTLINE
  return static_cast<size_t>((static_cast<double>(this->simulated_size()) / 1.0e6) * nlefs_per_mbp);
}
size_t Chromosome::nbarriers() const { return this->_barriers.size(); }

size_t Chromosome::num_valid_barriers() const {
  return static_cast<size_t>(
      std::count_if(this->get_barriers().begin(), this->get_barriers().end(),
                    [](const auto& b) { return b.strand == '-' || b.strand == '+'; }));
}

absl::Span<const bed::BED> Chromosome::get_barriers() const { return this->_barriers.data(); }

void Chromosome::increment_contacts(bp_t pos1, bp_t pos2, bp_t bin_size) {
  this->_contacts->increment((pos1 - this->_start) / bin_size, (pos2 - this->_start) / bin_size);
}

void Chromosome::increment_contacts(bp_t bin1, bp_t bin2) {
  this->_contacts->increment(bin1, bin2);
}

void Chromosome::allocate_contacts(bp_t bin_size, bp_t diagonal_width) {
  const auto ncols = (this->simulated_size() / bin_size) + (this->simulated_size() % bin_size != 0);
  const auto nrows =
      std::min(ncols, (diagonal_width / bin_size) + (diagonal_width % bin_size != 0));
  this->_contacts = std::make_shared<ContactMatrix<contacts_t>>(nrows, ncols);
}

void Chromosome::deallocate_contacts() { this->_contacts = nullptr; }

const Chromosome::contact_matrix_t& Chromosome::contacts() const {
  assert(this->_contacts);  // NOLINT
  return *this->_contacts;
}

uint64_t Chromosome::hash(uint64_t seed, size_t cell_id) {
  auto handle_errors = [&](const auto& status) {
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

Genome::Genome(const std::filesystem::path& path_to_chrom_sizes,
               const std::filesystem::path& path_to_extr_barriers,
               const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms)
    : _chromosomes(instantiate_genome(path_to_chrom_sizes, path_to_extr_barriers,
                                      path_to_chrom_subranges, keep_all_chroms)) {}

absl::btree_set<Chromosome> Genome::import_chromosomes(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms) {
  assert(!path_to_chrom_sizes.empty());  // NOLINT

  // Parse chrom subranges from BED. We parse everything at once to deal with duplicate entries
  const auto chrom_ranges = [&]() {
    absl::btree_map<std::string, bed::BED> ranges;
    if (!path_to_chrom_subranges.empty()) {
      for (auto&& record : bed::Parser(path_to_chrom_subranges).parse_all()) {
        const auto chrom_name = record.chrom;
        const auto [node, new_insertion] = ranges.try_emplace(chrom_name, std::move(record));
        if (!new_insertion) {
          throw std::runtime_error(
              fmt::format(FMT_STRING("Found more than one entry for chromosome '{}' in file {}"),
                          chrom_name, path_to_chrom_subranges));
        }
      }
    }
    return ranges;
  }();

  // Parse chrom. sizes and build the set of chromosome to be simulated.
  // When the BED file with the chrom. subranges is not provided, all the chromosome in the
  // chrom.sizes file will be selected and returned. When a BED file with the chrom. subranges is
  // available, then only chromosomes that are present in both files will be selected. Furthermore
  // we are also checking that the subrange lies within the genomic coordinates specified in the
  // chrom. sizes file
  auto id = 0UL;
  absl::btree_set<Chromosome> chromosomes;
  if (chrom_ranges.empty()) {
    for (auto record : chrom_sizes::Parser(path_to_chrom_sizes).parse_all()) {
      record.thick_start = record.chrom_start;
      record.thick_end = record.chrom_end;
      chromosomes.emplace(id++, std::move(record));
    }
    return chromosomes;
  }

  for (auto record : chrom_sizes::Parser(path_to_chrom_sizes).parse_all()) {
    const auto match = chrom_ranges.find(record.chrom);
    if (match != chrom_ranges.end()) {
      const auto& begin_pos = match->second.chrom_start;
      const auto& end_pos = match->second.chrom_end;
      if (begin_pos < record.chrom_start || end_pos > record.chrom_end) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("According to the chrom.sizes file {}, chromosome '{}' should have a size "
                       "of '{}', but the subrange specified in BED file {} extends past this "
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

    if (keep_all_chroms || match != chrom_ranges.end()) {
      chromosomes.emplace(id++, record.chrom, record.thick_start, record.thick_end, record.size());
    }
  }

  return chromosomes;
}

size_t Genome::import_barriers(absl::btree_set<Chromosome>& chromosomes,
                               const std::filesystem::path& path_to_extr_barriers) {
  assert(!chromosomes.empty());            // NOLINT
  assert(!path_to_extr_barriers.empty());  // NOLINT
  size_t nbarriers = 0;
  /*
  absl::flat_hash_map<std::string_view, Chromosome*> tmp_chrom_names(
      static_cast<size_t>(chromosomes.size()));
  for (auto& chrom : chromosomes) {
    tmp_chrom_names.emplace(chrom.name(), &chrom);
  }
   */

  // Parse all the records from the BED file. The parser will throw in case of duplicates.
  const auto barriers = bed::Parser(path_to_extr_barriers).parse_all_in_interval_tree();

  for (auto& chrom : chromosomes) {
    if (const auto chrom_name = std::string{chrom.name()}; barriers.contains(chrom_name)) {
      for (const auto& record : barriers.at(chrom_name).data()) {
        if (record.score < 0 || record.score > 1) {
          throw std::runtime_error(
              fmt::format("Invalid score field detected for record {}[{}-{}]: expected a score "
                          "between 0 and 1, got {:.4g}.",
                          record.chrom, record.chrom_start, record.chrom_end, record.score));
        }
        chrom.add_extrusion_barrier(record);
        ++nbarriers;
      }
    }
  }
  return nbarriers;
}

std::vector<ExtrusionBarrier> Genome::generate_vect_of_barriers(std::string_view chrom_name,
                                                                double ctcf_prob_occ_to_occ,
                                                                double ctcf_prob_nocc_to_nocc) {
  std::vector<ExtrusionBarrier> barriers;
  size_t barriers_skipped = 0;

  const auto* chrom = [&]() -> const Chromosome* {
    for (const auto& c : this->_chromosomes) {
      if (c.name() == chrom_name) {
        return &c;
      }
    }
    return nullptr;
  }();

  assert(chrom);  // NOLINT
  for (const auto& record : chrom->get_barriers()) {
    // Only instantiate barriers with a known motif direction.
    if (record.strand != '+' && record.strand != '-') {
      ++barriers_skipped;
      continue;
    }
    const auto pos = (record.chrom_start + record.chrom_end + 1) / 2;
    if (pos < chrom->start_pos() || pos >= chrom->end_pos()) {
      // Barrier lies outside of the genomic regions to be simulated
      ++barriers_skipped;
      continue;
    }

    if (record.score != 0) {
      // When the score field is zero (i.e. when the extr. barrier does not have a custom
      // occupancy), use the occupancy specified through the CLI
      const auto pblock = record.score;
      const auto pno = ctcf_prob_nocc_to_nocc;
      const auto poo =
          ExtrusionBarrier::compute_blocking_to_blocking_transition_probabilities_from_pblock(
              pblock, pno);

      barriers.emplace_back(pos, poo, pno, record.strand);
    } else {
      barriers.emplace_back(pos, ctcf_prob_occ_to_occ, ctcf_prob_nocc_to_nocc, record.strand);
    }
  }

  fmt::print(stderr,
             FMT_STRING("Instantiated {} extr. barriers for '{}' ({} barriers were skipped).\n"),
             barriers.size(), chrom->name(), barriers_skipped);

  return barriers;
}

absl::btree_set<Chromosome> Genome::instantiate_genome(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_extr_barriers,
    const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms) {
  auto chroms = import_chromosomes(path_to_chrom_sizes, path_to_chrom_subranges, keep_all_chroms);
  import_barriers(chroms, path_to_extr_barriers);

  return chroms;
}

Genome::iterator Genome::begin() { return this->_chromosomes.begin(); }
Genome::iterator Genome::end() { return this->_chromosomes.end(); }

Genome::const_iterator Genome::begin() const { return this->_chromosomes.cbegin(); }
Genome::const_iterator Genome::end() const { return this->_chromosomes.cend(); }

size_t Genome::size() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
      [](auto accumulator, const auto& chrom) { return accumulator + chrom.size(); });
}

size_t Genome::simulated_size() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
                         [](auto accumulator, const auto& chrom) {
                           return accumulator + (chrom.end_pos() - chrom.start_pos());
                         });
}

const Chromosome& Genome::chromosome_with_longest_name() const {
  assert(!this->_chromosomes.empty());  // NOLINT
  return *std::max_element(
      this->_chromosomes.begin(), this->_chromosomes.end(),
      [&](const auto& c1, const auto& c2) { return c1.name().size() < c2.name().size(); });
}

const Chromosome& Genome::longest_chromosome() const {
  assert(!this->_chromosomes.empty());  // NOLINT
  return *std::max_element(
      this->_chromosomes.begin(), this->_chromosomes.end(),
      [&](const auto& c1, const auto& c2) { return c1.simulated_size() < c2.simulated_size(); });
}

const Chromosome& Genome::chromosome_with_max_nbarriers() const {
  assert(!this->_chromosomes.empty());  // NOLINT
  return *std::max_element(this->_chromosomes.begin(), this->_chromosomes.end(),
                           [&](const auto& c1, const auto& c2) {
                             return c1.num_valid_barriers() < c2.num_valid_barriers();
                           });
}

size_t Genome::max_target_contacts(size_t bin_size, size_t diagonal_width,
                                   double target_contact_density, size_t simulation_iterations,
                                   double lef_fraction_contact_sampling, double nlefs_per_mbp,
                                   size_t ncells) const {
  const auto& chrom = this->longest_chromosome();
  if (target_contact_density == 0.0) {
    const auto nlefs = static_cast<size_t>(
        std::round(nlefs_per_mbp * (static_cast<double>(chrom.simulated_size()) / 1.0e6)));
    return static_cast<size_t>(
        (static_cast<double>(simulation_iterations * nlefs) * lef_fraction_contact_sampling) /
        static_cast<double>(ncells));
  }

  const auto npixels = ((chrom.simulated_size() + bin_size - 1) / bin_size) *
                       ((diagonal_width + bin_size - 1) / bin_size);

  return static_cast<size_t>(std::round((static_cast<double>(npixels) * target_contact_density) /
                                        static_cast<double>(ncells)));
}

}  // namespace modle
