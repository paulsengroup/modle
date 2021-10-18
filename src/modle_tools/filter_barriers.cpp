// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/strings/ascii.h>  // for AsciiStrToLower
#include <fmt/format.h>          // for print, FMT_STRING

#include <cassert>    // for assert
#include <cstdio>     // for stdout
#include <memory>     // for allocator_traits<>::value_type
#include <stdexcept>  // for logic_error
#include <string>     // for string, operator==
#include <vector>     // for vector

#include "modle/bed.hpp"            // for BED_tree, Parser, formatter<>::format, formatter<>::parse
#include "modle/common/common.hpp"  // for u64, u32, u8
#include "modle_tools/config.hpp"   // for filter_barrier_config
#include "modle_tools/tools.hpp"    // for filter_barriers_subcmd

namespace modle::tools {

void filter_barriers_intersection(const modle::tools::filter_barrier_config& c) {
  assert(!c.path_to_bed_files_for_filtering.empty());  // NOLINT
  std::vector<bed::BED_tree<std::string, u64>> intervals(c.path_to_bed_files_for_filtering.size());

  for (usize i = 0; i < c.path_to_bed_files_for_filtering.size(); ++i) {
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
      fmt::print(stdout, FMT_STRING("{}\n"), record);
    }
  }
}

// TODO: Find a better name than pairwise-intersection
void filter_barriers_pairwise_intersection(const modle::tools::filter_barrier_config& c) {
  assert(!c.path_to_bed_files_for_filtering.empty());  // NOLINT
  std::vector<bed::BED_tree<std::string, u64>> intervals(c.path_to_bed_files_for_filtering.size());

  for (usize i = 0; i < c.path_to_bed_files_for_filtering.size(); ++i) {
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
        fmt::print(stdout, FMT_STRING("{}\n"), record);
        break;
      }
    }
  }
}

void filter_barriers_subcmd(const modle::tools::filter_barrier_config& c) {
  if (absl::AsciiStrToLower(c.filtering_criterion) == "intersection") {
    filter_barriers_intersection(c);
  } else if (absl::AsciiStrToLower(c.filtering_criterion) == "pairwise-intersection") {
    filter_barriers_pairwise_intersection(c);
  } else {
    throw std::logic_error("Unreachable code");
  }
}

}  // namespace modle::tools
