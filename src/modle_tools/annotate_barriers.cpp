// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>

#include <cassert>
#include <ctgmath>

#include "modle/bed/bed.hpp"
#include "modle/bigwig/bigwig.hpp"
#include "modle_tools/modle_tools_config.hpp"
#include "modle_tools/tools.hpp"

namespace modle::tools {

[[nodiscard]] static constexpr double logit(double x, double scaling_factor) noexcept {
  assert(scaling_factor != 0);
  return 1.0 / (1.0 + std::exp(-x / scaling_factor));
}

static void process_interval(io::bigwig::Reader& bw, bed::BED interval, double mu, double lb,
                             double ub, bool clamp_occupancies) {
  assert(lb < ub);
  const auto peak =
      bw.stats<bwStatsType::max>(interval.chrom, interval.chrom_start, interval.chrom_end);
  interval.score = logit(peak, mu);
  if (clamp_occupancies) {
    interval.score = std::clamp(interval.score, lb, ub);
  }

  if (interval.score >= lb && interval.score <= ub) {
    fmt::print(FMT_COMPILE("{:bed6}\n"), interval);
  }
}

void annotate_barriers_subcmd(const annotate_barriers_config& c) {
  auto bw = io::bigwig::Reader(c.path_to_bigwig);
  auto bed_parser = bed::Parser(c.path_to_bed, bed::BED::BED6);
  while (true) {
    auto interval = bed_parser.parse_next();
    if (!interval) {
      break;
    }

    process_interval(bw, std::move(interval), c.scaling_factor, c.occupancy_lb, c.occupancy_ub,
                     c.clamp_occupancy);
  }
}

}  // namespace modle::tools
