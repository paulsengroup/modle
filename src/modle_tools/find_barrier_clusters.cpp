#include <fmt/format.h>
#include <fmt/os.h>

#include "modle/bed.hpp"
#include "modle_tools/config.hpp"
#include "modle_tools/tools.hpp"

namespace modle::tools {
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
      c.max_cluster_size == 0 ? std::numeric_limits<size_t>::max() : c.max_cluster_size;
  const auto min_cluster_span = c.min_cluster_span;
  const auto max_cluster_span =
      c.max_cluster_span == 0 ? std::numeric_limits<bp_t>::max() : c.max_cluster_span;

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

      const auto cluster_span = cluster.chrom_end - cluster.chrom_start;
      // Print the cluster in BED format
      if (cluster_span != 0 && cluster_span >= min_cluster_span &&
          cluster_span < max_cluster_span && cluster_size >= min_cluster_size &&
          cluster_size < max_cluster_size) {
        cluster.name = fmt::format(FMT_STRING("cluster_{:07d}"), cluster_id);
        cluster.score = static_cast<double>(cluster_size);
        if (fp) {
          fp->print(FMT_STRING("{:bed5}\n"), cluster);
        } else {
          fmt::print(stdout, FMT_STRING("{:bed5}\n"), cluster);
        }
      } else if (!c.quiet && cluster_span != 0) {
        fmt::print(stderr, FMT_STRING("Warning: Skipping a cluster {}:{}-{}. Reason:"), b1.chrom,
                   cluster.chrom_start, cluster.chrom_end);
        if (cluster_span < min_cluster_span) {
          fmt::print(stderr, FMT_STRING(" Cluster span is too small ({} < {});"), cluster_span,
                     min_cluster_span);
        } else if (cluster_span >= max_cluster_span) {
          fmt::print(stderr, FMT_STRING(" Cluster span is too large ({} >= {});"), cluster_span,
                     max_cluster_span);
        }

        if (cluster_size < min_cluster_size) {
          fmt::print(stderr, FMT_STRING(" Cluster does not contain enough barriers ({} < {});"),
                     cluster_size, min_cluster_size);
        }
        if (cluster_size >= max_cluster_size) {
          fmt::print(stderr, FMT_STRING(" Cluster contains too many barriers ({} >= {});"),
                     cluster_size, max_cluster_size);
        }
        fmt::print(stderr, "\n");
      }
      ++cluster_id;
      cluster_size = 0;
      cluster.chrom_start = cluster.chrom_end;
    }
  }
}
}  // namespace modle::tools
