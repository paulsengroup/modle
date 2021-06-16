#include <absl/strings/ascii.h>  // for AsciiStrToLower
#include <fmt/format.h>          // for print, FMT_STRING, format

#include <cassert>  // for assert
#include <cstddef>  // for size_t
#include <cstdint>  // for uint64_t
#include <vector>   // for vector

#include "modle/bed.hpp"           // for Parser
#include "modle_tools/config.hpp"  // for config
#include "modle_tools/tools.hpp"

namespace modle::tools {

void filter_barriers_intersection(const modle::tools::config& c) {
  assert(!c.path_to_bed_files_for_filtering.empty());  // NOLINT
  std::vector<bed::BED_tree<std::string, uint64_t>> intervals(
      c.path_to_bed_files_for_filtering.size());

  for (auto i = 0UL; i < c.path_to_bed_files_for_filtering.size(); ++i) {
    intervals[i].emplace(
        bed::Parser(c.path_to_bed_files_for_filtering[i], c.bed_dialect, c.strict_bed_validation)
            .parse_all());
  }
  for (auto& i : intervals) {
    i.index();
  }

  for (auto&& record :
       bed::Parser(c.path_to_extrusion_barrier_motifs_bed, c.bed_dialect, c.strict_bed_validation)
           .parse_all()) {
    auto print = true;
    for (const auto& i : intervals) {
      if (!i.contains_overlap(record)) {
        print = false;
        break;
      }
    }
    if (print) {
      fmt::print(stdout, "{}\n", record.to_string());
    }
  }
}

// TODO: Find a better name than pairwise-intersection
void filter_barriers_pairwise_intersection(const modle::tools::config& c) {
  assert(!c.path_to_bed_files_for_filtering.empty());  // NOLINT
  std::vector<bed::BED_tree<std::string, uint64_t>> intervals(
      c.path_to_bed_files_for_filtering.size());

  for (auto i = 0UL; i < c.path_to_bed_files_for_filtering.size(); ++i) {
    intervals[i].emplace(
        bed::Parser(c.path_to_bed_files_for_filtering[i], c.bed_dialect, c.strict_bed_validation)
            .parse_all());
  }

  for (auto& i : intervals) {
    i.index();
  }

  for (auto&& record :
       bed::Parser(c.path_to_extrusion_barrier_motifs_bed, c.bed_dialect, c.strict_bed_validation)
           .parse_all()) {
    for (const auto& i : intervals) {
      if (!i.contains_overlap(record)) {
        fmt::print(stdout, "{}\n", record.to_string());
        break;
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
