#include "modle/setup.hpp"

#include <absl/container/flat_hash_map.h>
#include <fmt/format.h>
#include <fmt/ostream.h>  // required to format std::filesystem::path

#include <cassert>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>

#include "modle/chrom_sizes.hpp"
#include "modle/dna.hpp"

namespace modle {

/// Import chromosomes from a chrom.sizes file

//! When \p path_to_extr_barriers is non-empty, import the intersection of the chromosomes present
//! in the chrom.sizes and BED files. The optional BED file can be used to instruct ModLE to
//! simulate loop extrusion on a sub-region of a chromosome from the chrom.sizes file.
[[nodiscard]] absl::btree_set<Chromosome, Chromosome::Comparator> import_chromosomes(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms);

/// Parse a BED file containing the genomic coordinates of extrusion barriers and add them to the
/// Genome
size_t import_barriers(absl::btree_set<Chromosome, Chromosome::Comparator>& chromosomes,
                       const std::filesystem::path& path_to_extr_barriers);

Genome::Genome(const std::filesystem::path& path_to_chrom_sizes,
               const std::filesystem::path& path_to_extr_barriers,
               const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms)
    : _chromosomes(instantiate_genome(path_to_chrom_sizes, path_to_extr_barriers,
                                      path_to_chrom_subranges, keep_all_chroms)) {}

absl::btree_set<Chromosome, Chromosome::Comparator> import_chromosomes(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms) {
  assert(!path_to_chrom_sizes.empty());  // NOLINT

  // Parse chrom subranges from BED. We parse everything at once to deal with duplicate entries
  absl::flat_hash_map<std::string, std::pair<uint64_t, uint64_t>> chrom_ranges;
  if (!path_to_chrom_subranges.empty()) {
    for (auto&& record : modle::bed::Parser(path_to_chrom_subranges).parse_all()) {
      chrom_ranges.emplace(std::move(record.chrom),
                           std::make_pair(record.chrom_start, record.chrom_end));
    }
  }

  // Parse chrom. sizes and build the set of chromosome to be simulated.
  // When the BED file with the chrom. subranges is not provided, all the chromosome in the
  // chrom.sizes file will be selected and returned. When a BED file with the chrom. subranges is
  // available, then only chromosomes that are present in both files will be selected. Furthermore
  // we are also checking that the subrange lies within the genomic coordinates specified in the
  // chrom. sizes file
  auto chrom_id = 0UL;
  absl::btree_set<Chromosome, Chromosome::Comparator> chromosomes;
  for (auto&& chrom : chrom_sizes::Parser(path_to_chrom_sizes).parse_all()) {
    if (auto match = chrom_ranges.find(chrom.name); match != chrom_ranges.end()) {
      const auto& range_start = match->second.first;
      const auto& range_end = match->second.second;
      if (range_start < chrom.start || range_end > chrom.end) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("According to the chrom.sizes file {}, chromosome '{}' should have a size "
                       "of '{}', but the subrange specified in BED file {} extends past this "
                       "region: range {}:{}-{} does not fit in range {}:{}-{}"),
            path_to_chrom_sizes, chrom.name, chrom.end, path_to_chrom_subranges, chrom.name,
            range_start, range_end, chrom.name, chrom.start, chrom.end));
      }
      chrom.start = range_start;
      chrom.end = range_end;
      chromosomes.emplace(std::move(chrom), chrom_id++);
    } else if (chrom_ranges.empty() || keep_all_chroms) {
      chromosomes.emplace(std::move(chrom), chrom_id++);
    }
  }
  return chromosomes;
}

size_t import_barriers(absl::btree_set<Chromosome, Chromosome::Comparator>& chromosomes,
                       const std::filesystem::path& path_to_extr_barriers) {
  assert(!chromosomes.empty());
  assert(!path_to_extr_barriers.empty());
  size_t nbarriers = 0;
  // Parse all the records from the BED file. parse_all() will throw in case of duplicates.
  for (auto&& record : modle::bed::Parser(path_to_extr_barriers).parse_all()) {
    if (record.score < 0 || record.score > 1) {
      throw std::runtime_error(
          fmt::format("Invalid score field detected for record {}[{}-{}]: expected a score "
                      "between 0 and 1, got {:.4g}.",
                      record.chrom, record.chrom_start, record.chrom_end, record.score));
    }

    if (auto match = chromosomes.find(record.chrom); match != chromosomes.end()) {
      match->add_extrusion_barrier(record);
      ++nbarriers;
    }
  }
  return nbarriers;
}

std::vector<ExtrusionBarrier> Genome::generate_vect_of_barriers(std::string_view chrom_name,
                                                                double ctcf_prob_occ_to_occ,
                                                                double ctcf_prob_nocc_to_nocc) {
  std::vector<ExtrusionBarrier> barriers;
  size_t barriers_skipped = 0;

  const auto& chrom = this->_chromosomes.find(chrom_name);
  assert(chrom != this->_chromosomes.end());  // NOLINT
  for (const auto& b : chrom->get_barriers()) {
    // Only instantiate barriers with a known motif direction.
    if (b.strand != '+' && b.strand != '-') MODLE_UNLIKELY {
        ++barriers_skipped;
        continue;
      }
    const auto pos = (b.chrom_start + b.chrom_end + 1) / 2;
    if (pos < chrom->start_pos() || pos >= chrom->end_pos()) {
      // Barrier lies outside of the genomic regions to be simulated
      ++barriers_skipped;
      continue;
    }

    if (b.score != 0) {
      // When the score field is zero (i.e. when the extr. barrier does not have a custom
      // occupancy), use the occupancy specified through the CLI
      const auto pblock = b.score;
      const auto pno = ctcf_prob_nocc_to_nocc;
      const auto poo =
          ExtrusionBarrier::compute_blocking_to_blocking_transition_probabilities_from_pblock(
              pblock, pno);

      barriers.emplace_back(pos, poo, pno, b.strand);
    } else {
      barriers.emplace_back(pos, ctcf_prob_occ_to_occ, ctcf_prob_nocc_to_nocc, b.strand);
    }
  }

  fmt::print(stderr,
             FMT_STRING("Instantiated {} extr. barriers for '{}' ({} barriers were skipped).\n"),
             barriers.size(), chrom->name(), barriers_skipped);

  return barriers;
}

absl::btree_set<Chromosome, Chromosome::Comparator> Genome::instantiate_genome(
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
