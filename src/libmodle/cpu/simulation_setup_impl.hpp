#pragma once

#include <absl/container/btree_set.h>      // for btree_iterator, btree_set
#include <absl/container/flat_hash_map.h>  // for flat_hash_map, BitMask, raw_hash_set<>::iterator
#include <fmt/format.h>                    // for format, print, FMT_STRING
#include <fmt/ostream.h>                   // for formatbuf<>::int_type

#include <algorithm>   // for max
#include <cassert>     // for assert
#include <cstddef>     // for size_t
#include <cstdint>     // for uint64_t
#include <cstdio>      // for stderr
#include <filesystem>  // for path, operator<<
#include <sstream>     // for basic_stringbuf<>::int_type, basic_stringbuf<>...
#include <stdexcept>   // for runtime_error
#include <string>      // for string, basic_string
#include <utility>     // for pair, move, make_pair
#include <vector>      // for vector

#include "modle/bed.hpp"                 // for BED, Parser
#include "modle/chrom_sizes.hpp"         // for ChromSize, Parser
#include "modle/common.hpp"              // for MODLE_LIKELY
#include "modle/dna.hpp"                 // for Chromosome
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier

namespace modle {

// TODO Add flag to import skip chrom without barriers
Simulation::Genome Simulation::instantiate_genome(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_extr_barriers,
    const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms) {
  auto genome =
      Simulation::import_chromosomes(path_to_chrom_sizes, path_to_chrom_subranges, keep_all_chroms);
  Simulation::import_barriers(genome, path_to_extr_barriers);
  return genome;
}

Simulation::Genome Simulation::import_chromosomes(
    const std::filesystem::path& path_to_chrom_sizes,
    const std::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms) {
  assert(!path_to_chrom_sizes.empty());  // NOLINT
  Genome chroms;

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
      chroms.emplace(std::move(chrom));
    } else if (chrom_ranges.empty() || keep_all_chroms) {
      chroms.emplace(std::move(chrom));
    }
  }
  return chroms;
}

size_t Simulation::import_barriers(Genome& genome,
                                   const std::filesystem::path& path_to_extr_barriers) {
  assert(!genome.empty());
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

    if (auto match = genome.find(record.chrom); match != genome.end()) {
      match->add_extrusion_barrier(record);
      ++nbarriers;
    }
  }
  return nbarriers;
}

std::vector<ExtrusionBarrier> Simulation::allocate_barriers(const Chromosome* const chrom) const {
  std::vector<ExtrusionBarrier> barriers;
  size_t barriers_skipped = 0;

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
      const auto pno = this->ctcf_not_occupied_self_prob;
      const auto poo =
          ExtrusionBarrier::compute_blocking_to_blocking_transition_probabilities_from_pblock(
              pblock, pno);

      barriers.emplace_back(pos, poo, pno, b.strand);
    } else {
      barriers.emplace_back(pos, this->ctcf_occupied_self_prob, this->ctcf_not_occupied_self_prob,
                            b.strand);
    }
  }

  fmt::print(stderr,
             FMT_STRING("Instantiated {} extr. barriers for '{}' ({} barriers were skipped).\n"),
             barriers.size(), chrom->name(), barriers_skipped);

  return barriers;
}
}  // namespace modle
