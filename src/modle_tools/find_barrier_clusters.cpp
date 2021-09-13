// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/container/btree_map.h>  // for btree_iterator
#include <absl/strings/str_format.h>   // for StrAppendFormat
#include <absl/types/span.h>           // for Span, MakeConstSpan
#include <fmt/format.h>                // for FMT_STRING, format, make_format_args, print, vform...
#include <fmt/os.h>                    // for ostream, output_file
#include <spdlog/spdlog.h>             // for warn

#include <boost/filesystem/path.hpp>  // for path
#include <cstdint>                    // for uint32_t, uint8_t
#include <cstdio>                     // for size_t, stdout
#include <exception>                  // for exception
#include <limits>                     // for numeric_limits
#include <memory>                     // for unique_ptr, make_unique
#include <string>                     // for string, basic_string
#include <type_traits>                // for add_const<>::type
#include <vector>                     // for vector

#include "modle/bed.hpp"            // for BED, BED_tree, formatter<>::format, formatter<>::p...
#include "modle/common/common.hpp"  // for bp_t
#include "modle/interval_tree.hpp"  // for IITree, IITree::IITree<I, T>, IITree::data_begin
#include "modle_tools/config.hpp"   // for find_barrier_clusters_config
#include "modle_tools/tools.hpp"    // for find_barrier_clusters_subcmd

namespace modle::tools {

void print_warnings(const bed::BED& cluster, const bed::BED& barrier1,
                    const find_barrier_clusters_config& c, const size_t cluster_span,
                    const size_t cluster_size, std::string& warning_buffer) {
  warning_buffer.clear();
  if (cluster_span < c.min_cluster_span) {
    absl::StrAppendFormat(&warning_buffer, " Cluster span is too small (%d < %d);", cluster_span,
                          c.min_cluster_span);
  } else if (cluster_span >= c.max_cluster_span) {
    absl::StrAppendFormat(&warning_buffer, " Cluster span is too large (%d >= %d);", cluster_span,
                          c.max_cluster_span);
  }

  if (cluster_size < c.min_cluster_size) {
    absl::StrAppendFormat(&warning_buffer, " Cluster does not contain enough barriers (%d < %d);",
                          cluster_size, c.min_cluster_size);
  }
  if (cluster_size >= c.max_cluster_size) {
    absl::StrAppendFormat(&warning_buffer, " Cluster contains too many barriers (%d >= %d);",
                          cluster_size, c.max_cluster_size);
  }
  spdlog::warn(FMT_STRING("Warning: Skipping a cluster {}:{}-{}. Reason: {}"), barrier1.chrom,
               cluster.chrom_start, cluster.chrom_end, warning_buffer);
}

void write_cluster(bed::BED& cluster, const size_t cluster_id, const size_t cluster_size,
                   std::unique_ptr<fmt::ostream>& fp) {
  cluster.name = fmt::format(FMT_STRING("cluster_{:07d}"), cluster_id);
  cluster.score = static_cast<double>(cluster_size);
  if (fp) {
    fp->print(FMT_STRING("{:bed5}\n"), cluster);
  } else {
    fmt::print(stdout, FMT_STRING("{:bed5}\n"), cluster);
  }
}

void write_single_barrier(bed::BED& buff, const bed::BED& barrier, const size_t barrier_id,
                          std::unique_ptr<fmt::ostream>& fp) {
  buff.chrom = barrier.chrom;
  buff.chrom_start = barrier.chrom_start;
  buff.chrom_end = barrier.chrom_end;
  buff.name = fmt::format(FMT_STRING("barrier_{:07d}"), barrier_id);
  buff.score = 1;
  if (fp) {
    fp->print(FMT_STRING("{:bed5}\n"), buff);
  } else {
    fmt::print(stdout, FMT_STRING("{:bed5}\n"), buff);
  }
}

void find_barrier_clusters_subcmd(const find_barrier_clusters_config& c) {
  const bed::BED_tree<> barrier_intervals_gw(c.path_to_input_barriers, bed::BED::BED3);
  const bed::BED_tree<> breaking_intervals_gw(c.path_to_breaking_points, bed::BED::BED3);
  auto fp = [c]() {
    if (c.path_to_output.empty()) {
      return std::unique_ptr<fmt::ostream>{nullptr};
    }
    return std::make_unique<fmt::ostream>(fmt::output_file(c.path_to_output.string()));
  }();

  const auto min_cluster_size = c.min_cluster_size;
  const auto max_cluster_size =
      c.max_cluster_size == 0 ? (std::numeric_limits<size_t>::max)() : c.max_cluster_size;
  const auto min_cluster_span = c.min_cluster_span;
  const auto max_cluster_span =
      c.max_cluster_span == 0 ? (std::numeric_limits<bp_t>::max)() : c.max_cluster_span;
  std::string warning_buffer;

  bed::BED cluster;  //{bed::BED::BED5};
  size_t cluster_id = 0;

  for (const auto& [chrom, barriers_itree] : barrier_intervals_gw) {
    // Create two spans corresponding to the barriers and breaking points mapping on the chrom that
    // is being processed
    const auto barriers =
        absl::MakeConstSpan(&(*barriers_itree.data_begin()), &(*barriers_itree.data_end()));
    const auto breaking_intervals = [&, chrom = chrom]() {
      if (breaking_intervals_gw.contains(chrom)) {
        return breaking_intervals_gw.at(chrom);
      }
      return bed::BED_tree<>::value_type{};
    }();

    cluster.chrom = chrom;
    cluster.chrom_start = 0;
    cluster.chrom_end = 0;
    size_t cluster_size = 0;

    auto process_cluster = [&](bed::BED& cluster_, const bed::BED& barrier,
                               const size_t barrier_id) {
      const auto cluster_span = cluster_.chrom_end - cluster_.chrom_start;
      if (cluster_span != 0 && cluster_span >= min_cluster_span &&
          cluster_span < max_cluster_span && cluster_size >= min_cluster_size &&
          cluster_size < max_cluster_size) {
        // Print the cluster in BED format if cluster size and span constraints are satisfied
        write_cluster(cluster_, cluster_id, cluster_size, fp);
      } else if (const auto barrier_span = barrier.chrom_end - barrier.chrom_start;
                 c.min_cluster_size == 1 && cluster_size <= 1 &&
                 barrier_span >= c.min_cluster_span && barrier_span < c.max_cluster_span) {
        // Handle case where --min-cluster-size=1
        write_single_barrier(cluster_, barrier, barrier_id, fp);
      } else if (!c.quiet && cluster_span != 0) {
        print_warnings(cluster_, barrier, c, cluster_span, cluster_size, warning_buffer);
      }
    };

    // Look for uninterrupted clusters of barriers
    // Clusters are identified by initially considering a window starting at the position of the
    // first extrusion barrier that has yet to be processed, and extends for a width equal to
    // c.extension_window.
    // This window will be extended by c.extension_window bp until extending the window does not
    // yield an increase in cluster size (i.e. number of barriers belonging to a cluster).
    // Cluster extension is also stopped when overlapping with one or more breaking intervals (e.g.
    // genes, promoters, enhancer etc.)
    for (size_t i = 1; i < barriers.size(); ++i) {
      const auto& b1 = barriers[i - 1];
      const auto& b2 = barriers[i];
      if (b2.chrom_end - b1.chrom_start < c.extension_window &&
          !breaking_intervals.overlaps_with(b1.chrom_start, b2.chrom_end)) {
        if (cluster_size == 0) {
          cluster.chrom_start = b1.chrom_start;
          cluster_size = 1;
        }
        cluster.chrom_end = b2.chrom_end + 1;
        ++cluster_size;
        continue;
      }

      process_cluster(cluster, b1, i - 1);
      ++cluster_id;
      cluster_size = 0;
      cluster.chrom_start = cluster.chrom_end;
    }
    // Handle last barrier mapping on a chromosome
    process_cluster(cluster, barriers.back(), barriers.size() - 1);
  }
}
}  // namespace modle::tools
