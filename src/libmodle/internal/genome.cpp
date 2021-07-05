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

Chromosome::Chromosome(size_t id, const bed::BED& chrom,
                       const IITree<bp_t, ExtrusionBarrier>& barriers, bool ok_)
    : Chromosome(id, chrom.chrom, chrom.thick_start, chrom.thick_end, chrom.size(), barriers, ok_) {
}

Chromosome::Chromosome(size_t id, const bed::BED& chrom, bool ok_)
    : _name(chrom.chrom),
      _start(chrom.chrom_start),
      _end(chrom.chrom_end),
      _size(chrom.size()),
      _id(id),
      _ok(ok_) {}

Chromosome::Chromosome(size_t id, const bed::BED& chrom, IITree<bp_t, ExtrusionBarrier>&& barriers,
                       bool ok_)
    : Chromosome(id, chrom.chrom, chrom.thick_start, chrom.thick_end, chrom.size(), barriers, ok_) {
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_SIGN_CONVERSION
DISABLE_WARNING_CONVERSION
Chromosome::Chromosome(size_t id, std::string_view chrom_name, bp_t chrom_start, bp_t chrom_end,
                       bp_t chrom_size, bool ok_)
    : _name(chrom_name),
      _start(chrom_start),
      _end(chrom_end),
      _size(chrom_size),
      _id(id),
      _ok(ok_) {
  assert(chrom_start <= chrom_end);               // NOLINT
  assert(chrom_end - chrom_start <= chrom_size);  // NOLINT
}

Chromosome::Chromosome(size_t id, std::string_view chrom_name, bp_t chrom_start, bp_t chrom_end,
                       bp_t chrom_size, const IITree<bp_t, ExtrusionBarrier>& barriers, bool ok_)
    : _name(chrom_name),
      _start(chrom_start),
      _end(chrom_end),
      _size(chrom_size),
      _id(id),
      _barriers(barriers),
      _ok(ok_) {
  assert(chrom_start <= chrom_end);               // NOLINT
  assert(chrom_end - chrom_start <= chrom_size);  // NOLINT

  _barriers.make_BST();
}

Chromosome::Chromosome(size_t id, std::string_view chrom_name, bp_t chrom_start, bp_t chrom_end,
                       bp_t chrom_size, IITree<bp_t, ExtrusionBarrier>&& barriers, bool ok_)
    : _name(chrom_name),
      _start(chrom_start),
      _end(chrom_end),
      _size(chrom_size),
      _id(id),
      _barriers(std::move(barriers)),
      _ok(ok_) {
  assert(chrom_start <= chrom_end);               // NOLINT
  assert(chrom_end - chrom_start <= chrom_size);  // NOLINT

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
      _contacts(std::make_unique<contact_matrix_t>(*other._contacts)),
      _features(other._features),
      _ok(other._ok) {
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
  _features = other._features;
  _ok = other._ok;

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

void Chromosome::add_extrusion_barrier(const bed::BED& record, const double ctcf_prob_occ_to_occ,
                                       const double ctcf_prob_nocc_to_nocc) {
  if (record.strand != '+' && record.strand != '-') {
    return;
  }
  const auto pos = (record.chrom_start + record.chrom_end + 1) / 2;
  if (pos < this->start_pos() || pos >= this->end_pos()) {
    // Barrier lies outside of the genomic regions to be simulated
    return;
  }

  if (record.score != 0) {
    // When the score field is zero (i.e. when the extr. barrier does not have a custom
    // occupancy), use the occupancy specified through the CLI
    const auto pblock = record.score;
    const auto pno = ctcf_prob_nocc_to_nocc;
    const auto poo =
        ExtrusionBarrier::compute_blocking_to_blocking_transition_probabilities_from_pblock(pblock,
                                                                                            pno);

    this->_barriers.emplace(record.chrom_start, record.chrom_end,
                            ExtrusionBarrier{pos, poo, pno, record.strand});
  } else {
    this->_barriers.emplace(
        record.chrom_start, record.chrom_end,
        ExtrusionBarrier{pos, ctcf_prob_occ_to_occ, ctcf_prob_nocc_to_nocc, record.strand});
  }
}

size_t Chromosome::id() const { return this->_id; }
std::string_view Chromosome::name() const { return this->_name; }
const char* Chromosome::name_cstr() const { return this->_name.c_str(); }

bool Chromosome::ok() const { return this->_ok; }

size_t Chromosome::num_lefs(double nlefs_per_mbp) const {  // NOLINTNEXTLINE
  return static_cast<size_t>((static_cast<double>(this->simulated_size()) / 1.0e6) * nlefs_per_mbp);
}
size_t Chromosome::num_barriers() const { return this->_barriers.size(); }

const IITree<bp_t, ExtrusionBarrier>& Chromosome::barriers() const { return this->_barriers; }
IITree<bp_t, ExtrusionBarrier>& Chromosome::barriers() { return this->_barriers; }
absl::Span<const Chromosome::bed_tree_value_t> Chromosome::get_features() const {
  return this->_features;
}

void Chromosome::increment_contacts(bp_t pos1, bp_t pos2, bp_t bin_size) {
  this->_contacts->increment((pos1 - this->_start) / bin_size, (pos2 - this->_start) / bin_size);
}

void Chromosome::increment_contacts(bp_t bin1, bp_t bin2) {
  this->_contacts->increment(bin1, bin2);
}

void Chromosome::allocate_contacts(bp_t bin_size, bp_t diagonal_width) {
  this->_contacts =
      std::make_shared<ContactMatrix<contacts_t>>(this->simulated_size(), diagonal_width, bin_size);
}

void Chromosome::deallocate_contacts() { this->_contacts = nullptr; }

const Chromosome::contact_matrix_t& Chromosome::contacts() const {
  assert(this->_contacts);  // NOLINT
  return *this->_contacts;
}

Chromosome::contact_matrix_t& Chromosome::contacts() {
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

Genome::Genome(const boost::filesystem::path& path_to_chrom_sizes,
               const boost::filesystem::path& path_to_extr_barriers,
               const boost::filesystem::path& path_to_chrom_subranges,
               const absl::Span<const boost::filesystem::path> paths_to_extra_features,
               const double ctcf_prob_occ_to_occ, const double ctcf_prob_nocc_to_nocc,
               bool keep_all_chroms)
    : _chromosomes(instantiate_genome(
          path_to_chrom_sizes, path_to_extr_barriers, path_to_chrom_subranges,
          paths_to_extra_features, ctcf_prob_occ_to_occ, ctcf_prob_nocc_to_nocc, keep_all_chroms)) {
}

absl::btree_set<Chromosome> Genome::import_chromosomes(
    const boost::filesystem::path& path_to_chrom_sizes,
    const boost::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms) {
  assert(!path_to_chrom_sizes.empty());  // NOLINT

  const auto t0 = absl::Now();
  if (path_to_chrom_subranges.empty()) {
    fmt::print(stderr, FMT_STRING("Importing chromosomes from file {}..."), path_to_chrom_sizes);
  } else {
    fmt::print(stderr, FMT_STRING("Importing chromosomes from files {} and {}..."),
               path_to_chrom_sizes, path_to_chrom_subranges);
  }
  auto print_status_update_on_return = [&](auto num_chromosomes) {
    fmt::print(stderr, FMT_STRING(" DONE!\nImported {} chromosomes in {}.\n"), num_chromosomes,
               absl::FormatDuration(absl::Now() - t0));
  };

  // Parse chrom subranges from BED. We parse everything at once to deal with duplicate entries
  const auto chrom_ranges = [&]() {
    absl::btree_map<std::string, bed::BED> ranges;
    if (!path_to_chrom_subranges.empty()) {
      for (const auto& record : bed::Parser(path_to_chrom_subranges, bed::BED::BED3).parse_all()) {
        // passing std::string{record.chrom} is required in order to make GCC 8.4 happy
        const auto [node, new_insertion] = ranges.try_emplace(std::string{record.chrom}, record);
        (void)node;
        if (!new_insertion) {
          throw std::runtime_error(
              fmt::format(FMT_STRING("Found more than one entry for chromosome '{}' in file {}"),
                          record.chrom, path_to_chrom_subranges));
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
    print_status_update_on_return(chromosomes.size());
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
      chromosomes.emplace(id++, record.chrom, record.thick_start, record.thick_end, record.size(),
                          match != chrom_ranges.end());
    }
  }

  print_status_update_on_return(chromosomes.size());
  return chromosomes;
}

size_t Genome::import_barriers(absl::btree_set<Chromosome>& chromosomes,
                               const boost::filesystem::path& path_to_extr_barriers,
                               const double ctcf_prob_occ_to_occ,
                               const double ctcf_prob_nocc_to_nocc) {
  assert(!chromosomes.empty());            // NOLINT
  assert(!path_to_extr_barriers.empty());  // NOLINT

  const auto t0 = absl::Now();
  fmt::print(stderr, FMT_STRING("Importing extrusion barriers from file {}..."),
             path_to_extr_barriers);

  size_t tot_num_barriers = 0;

  // Parse all the records from the BED file. The parser will throw in case of duplicates.
  const auto barriers =
      bed::Parser(path_to_extr_barriers, bed::BED::BED6).parse_all_in_interval_tree();

  for (auto& chrom : chromosomes) {
    if (const auto chrom_name = std::string{chrom.name()}; barriers.contains(chrom_name)) {
      for (const auto& record : barriers.at(chrom_name).data()) {
        if (record.score < 0 || record.score > 1) {
          throw std::runtime_error(
              fmt::format("Invalid score field detected for record {}[{}-{}]: expected a score "
                          "between 0 and 1, got {:.4g}.",
                          record.chrom, record.chrom_start, record.chrom_end, record.score));
        }
        chrom.add_extrusion_barrier(record, ctcf_prob_occ_to_occ, ctcf_prob_nocc_to_nocc);
        ++tot_num_barriers;
      }
    }
    chrom.barriers().make_BST();
  }
  fmt::print(stderr, FMT_STRING(" DONE!\nImported {} barriers in {}.\n"), tot_num_barriers,
             absl::FormatDuration(absl::Now() - t0));
  return tot_num_barriers;
}

size_t Genome::import_extra_features(absl::btree_set<Chromosome>& chromosomes,
                                     const boost::filesystem::path& path_to_extra_features) {
  assert(!chromosomes.empty());             // NOLINT
  assert(!path_to_extra_features.empty());  // NOLINT

  const auto t0 = absl::Now();
  fmt::print(stderr, FMT_STRING("Importing features from the following file: {}..."),
             path_to_extra_features);
  size_t num_features = 0;

  // Parse all the records from the BED file. The parser will throw in case of duplicates.
  const auto features =
      bed::Parser(path_to_extra_features, bed::BED::BED6).parse_all_in_interval_tree();

  for (auto& chrom : chromosomes) {
    if (const auto chrom_name = std::string{chrom.name()}; features.contains(chrom_name)) {
      const auto element = chrom._features.emplace_back(features.at(chrom_name));
      num_features += element.size();
    }
  }
  fmt::print(stderr, FMT_STRING(" DONE!\nImported {} features in {}.\n"), num_features,
             absl::FormatDuration(absl::Now() - t0));
  return num_features;
}

absl::btree_set<Chromosome> Genome::instantiate_genome(
    const boost::filesystem::path& path_to_chrom_sizes,
    const boost::filesystem::path& path_to_extr_barriers,
    const boost::filesystem::path& path_to_chrom_subranges,
    const absl::Span<const boost::filesystem::path> paths_to_extra_features,
    const double ctcf_prob_occ_to_occ, const double ctcf_prob_nocc_to_nocc, bool keep_all_chroms) {
  auto chroms = import_chromosomes(path_to_chrom_sizes, path_to_chrom_subranges, keep_all_chroms);
  import_barriers(chroms, path_to_extr_barriers, ctcf_prob_occ_to_occ, ctcf_prob_nocc_to_nocc);
  for (const auto& path_to_feature_bed : paths_to_extra_features) {
    import_extra_features(chroms, path_to_feature_bed);
  }

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
  return *std::max_element(
      this->_chromosomes.begin(), this->_chromosomes.end(),
      [&](const auto& c1, const auto& c2) { return c1.num_barriers() < c2.num_barriers(); });
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
