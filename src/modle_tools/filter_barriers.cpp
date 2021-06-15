#include <absl/container/flat_hash_map.h>  // for flat_hash_map, BitMask, operator!=
#include <absl/strings/ascii.h>            // for AsciiStrToLower
#include <fmt/format.h>                    // for print, FMT_STRING, format

#include <cassert>  // for assert

#include "modle/common/suppress_compiler_warnings.hpp"

DISABLE_WARNING_PUSH
DISABLE_WARNING_SIGN_CONVERSION
DISABLE_WARNING_SIGN_COMPARE
#include <cgranges/IITree.hpp>  // for IITree
DISABLE_WARNING_POP
#include <cstddef>  // for size_t
#include <cstdint>  // for uint64_t
#include <vector>   // for vector

#include "modle/bed.hpp"           // for Parser
#include "modle_tools/config.hpp"  // for config
#include "modle_tools/tools.hpp"

namespace modle::tools {

void filter_barriers_intersection(const modle::tools::config& c) {
  using IITree_t = IITree<uint64_t, uint8_t>;
  assert(!c.path_to_bed_files_for_filtering.empty());  // NOLINT
  absl::flat_hash_map<std::string, std::vector<IITree_t>> intervals(
      c.path_to_bed_files_for_filtering.size());

  for (auto i = 0UL; i < c.path_to_bed_files_for_filtering.size(); ++i) {
    for (const auto& record : bed::Parser(c.path_to_bed_files_for_filtering[i].string(),
                                          c.bed_dialect, c.strict_bed_validation)
                                  .parse_all(false)) {
      auto [node, new_insertion] = intervals.try_emplace(record.chrom, std::vector<IITree_t>{});
      if (new_insertion) {
        node->second.resize(c.path_to_bed_files_for_filtering.size());
      }
      node->second[i].add(record.chrom_start, record.chrom_end, uint8_t(0));
    }

    for (auto& [_, trees] : intervals) {
      trees[i].index();
    }
  }

  std::vector<size_t> buff;
  for (const auto& record :
       bed::Parser(c.path_to_extrusion_barrier_motifs_bed, c.bed_dialect, c.strict_bed_validation)
           .parse_all(false)) {
    if (auto node = intervals.find(record.chrom); node != intervals.end()) {
      auto print = true;
      for (const auto& tree : node->second) {
        if (tree.size() == 0) {
          print = false;
          break;
        }
        tree.overlap(record.chrom_start, record.chrom_end, buff);
        if (buff.empty()) {
          print = false;
          break;
        }
      }
      if (print) {
        fmt::print(stdout, "{}\n", record.to_string());
      }
    }
  }
}

// TODO: Find a better name than pairwise-intersection
void filter_barriers_pairwise_intersection(const modle::tools::config& c) {
  using IITree_t = IITree<uint64_t, uint64_t>;
  assert(!c.path_to_bed_files_for_filtering.empty());  // NOLINT
  absl::flat_hash_map<std::string, std::vector<IITree_t>> intervals(
      c.path_to_bed_files_for_filtering.size());

  for (auto i = 0UL; i < c.path_to_bed_files_for_filtering.size(); ++i) {
    for (const auto& record : bed::Parser(c.path_to_bed_files_for_filtering[i].string(),
                                          c.bed_dialect, c.strict_bed_validation)
                                  .parse_all(false)) {
      auto [node, new_insertion] = intervals.try_emplace(record.chrom, std::vector<IITree_t>{});
      if (new_insertion) {
        node->second.resize(c.path_to_bed_files_for_filtering.size());
      }
      node->second[i].add(record.chrom_start, record.chrom_end, uint8_t(0));
    }

    for (auto& [_, trees] : intervals) {
      trees[i].index();
    }
  }

  std::vector<size_t> buff;
  for (const auto& record :
       bed::Parser(c.path_to_extrusion_barrier_motifs_bed, c.bed_dialect, c.strict_bed_validation)
           .parse_all(false)) {
    if (auto node = intervals.find(record.chrom); node != intervals.end()) {
      for (const auto& tree : node->second) {
        if (tree.size() != 0) {
          tree.overlap(record.chrom_start, record.chrom_end, buff);
          if (!buff.empty()) {
            fmt::print(stdout, "{}\n", record.to_string());
          }
        }
      }
    }
  }
}

void filter_barriers_subcmd(const modle::tools::config& c) {
  if (absl::AsciiStrToLower(c.filtering_criterion) == "intersection") {
    filter_barriers_intersection(c);
  } else if (absl::AsciiStrToLower(c.filtering_criterion) == "pairwise-intersection") {
    filter_barriers_pairwise_intersection(c);
  } else {
    assert(false);
  }
}

}  // namespace modle::tools
