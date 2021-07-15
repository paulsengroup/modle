#include <fmt/format.h>
#include <fmt/os.h>

#include "modle/bed.hpp"
#include "modle_tools/config.hpp"
#include "modle_tools/tools.hpp"

namespace modle::tools {
void find_barrier_clusters_subcmd(const find_barrier_clusters_config& c) {
  const bed::BED_tree<> intervals(c.path_to_input_barriers);
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

  for (const auto& [chrom, itree] : intervals) {
    const auto barriers = absl::MakeConstSpan(&(*itree.data_begin()), &(*itree.data_end()));
    bp_t start_pos{};
    bp_t end_pos{};
    size_t cluster_size = 0;
    for (size_t i = 1; i < barriers.size(); ++i) {
      const auto& b1 = barriers[i - 1];
      const auto& b2 = barriers[i];
      if (b2.chrom_end - b1.chrom_start < c.extension_window) {
        if (cluster_size == 0) {
          start_pos = b1.chrom_start;
          cluster_size = 1;
        }
        end_pos = b2.chrom_end + 1;
        ++cluster_size;
        continue;
      }

      const auto cluster_span = end_pos - start_pos;
      if (cluster_span != 0 && cluster_span >= min_cluster_span &&
          cluster_span < max_cluster_span && cluster_size >= min_cluster_size &&
          cluster_size < max_cluster_size) {
        if (fp) {
          fp->print(FMT_STRING("{}\t{}\t{}\t{}\n"), b1.chrom, start_pos, end_pos, cluster_size);
        } else {
          fmt::print(stdout, FMT_STRING("{}\t{}\t{}\t{}\n"), b1.chrom, start_pos, end_pos,
                     cluster_size);
        }
      } else if (cluster_span != 0) {
        fmt::print(stderr, FMT_STRING("Warning: Skipping a cluster located at {}:{}-{}. Reason:"),
                   b1.chrom, start_pos, end_pos);
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
      cluster_size = 0;
      start_pos = end_pos;
    }
  }
}
}  // namespace modle::tools
