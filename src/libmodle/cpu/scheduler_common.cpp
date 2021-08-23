// clang-format off
#include "modle/simulation.hpp"
// clang-format on

#include <absl/container/flat_hash_set.h>  // for flat_hash_set, BitMask
#include <absl/strings/match.h>            // for StartsWith
#include <absl/strings/str_cat.h>          // for StrAppend
#include <absl/strings/str_format.h>       // for StrAppendFormat
#include <absl/time/clock.h>               // for Now
#include <absl/time/time.h>                // for FormatDuration, operator-, Time
#include <absl/types/span.h>               // for Span, MakeConstSpan
#include <fmt/format.h>                    // for format, make_format_args, vformat_to, FMT_STRING
#include <fmt/ostream.h>                   // for formatbuf<>::int_type
#include <spdlog/spdlog.h>                 // for info, error

#include <algorithm>                    // for min, max
#include <boost/filesystem/path.hpp>    // for path
#include <cassert>                      // for assert
#include <cstddef>                      // IWYU pragma: keep for size_t
#include <cstdint>                      // for uint32_t, uint8_t
#include <exception>                    // for exception, rethrow_exception
#include <stdexcept>                    // for runtime_error
#include <string>                       // for string
#include <thread_pool/thread_pool.hpp>  // for thread_pool
#include <vector>                       // for vector

#include "modle/bed.hpp"                 // for BED_tree, BED_tree::size, BED, BED::BED3, BED_...
#include "modle/common/utils.hpp"        // for parse_numeric_or_throw
#include "modle/compressed_io.hpp"       // for Reader
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier
#include "modle/genome.hpp"              // for Chromosome, Genome
#include "modle/interval_tree.hpp"       // for IITree, IITree::data_end, IITree::equal_range

namespace modle {

bed::BED_tree<> Simulation::import_deletions() const {
  spdlog::info(FMT_STRING("Importing deletions from file {}..."), this->path_to_deletion_bed);
  const auto t0 = absl::Now();
  auto deletions = bed::BED_tree<>{this->path_to_deletion_bed, bed::BED::BED3};
  spdlog::info(FMT_STRING("Imported {} deletions in {}"), deletions.size(),
               absl::FormatDuration(absl::Now() - t0));
  return deletions;
}

bed::BED_tree<> Simulation::generate_deletions() const {
  spdlog::info(FMT_STRING("Generating deletions for {} chromosomes..."), this->_genome.size());
  const auto t0 = absl::Now();
  bed::BED_tree<> deletions;
  for (const auto& chrom : this->_genome) {
    for (const auto& barrier : chrom.barriers().data()) {
      // Compute the span of the deletion
      const auto deletion_size_ = std::min(barrier.pos() + 1, this->deletion_size);
      const auto deletion_begin = barrier.pos() + 1 - deletion_size_;
      const auto deletion_end = deletion_begin + deletion_size_;
      deletions.insert(std::string{chrom.name()}, deletion_begin, deletion_end);
    }
  }
  deletions.index();
  spdlog::info(FMT_STRING("Generated {} deletions in {}"), deletions.size(),
               absl::FormatDuration(absl::Now() - t0));
  return deletions;
}

bool Simulation::advance_window(TaskPW& base_task, const Chromosome& chrom) const {
  base_task.active_window_start += this->diagonal_width;
  base_task.active_window_end =
      std::min(base_task.active_window_start + (2 * this->diagonal_width), chrom.end_pos());

  base_task.window_start = base_task.active_window_start - this->diagonal_width;
  base_task.window_end =
      std::min(base_task.active_window_end + this->diagonal_width, chrom.end_pos());

  return base_task.active_window_start < chrom.end_pos();
}

bool Simulation::map_features_to_window(TaskPW& base_task, const Chromosome& chrom) {
  const auto chrom_name = std::string{chrom.name()};
  const auto& features = chrom.get_features();
  auto print_status_update = [](const auto& t) {
    spdlog::info(FMT_STRING("Skipping {}[{}-{}]..."), t.chrom->name(), t.active_window_start,
                 t.active_window_end);
  };

  // Find all features of type 1 falling within the outer window. Skip over windows with 0
  // features
  auto [first_feat1, last_feat1] =
      features[0].equal_range(base_task.active_window_start, base_task.active_window_end);
  if (first_feat1 == features[0].data_end()) {
    print_status_update(base_task);
    return false;
  }

  // Find all features of type 2 falling within the outer window. Skip over windows with 0
  // features
  auto [first_feat2, last_feat2] =
      features[1].equal_range(base_task.active_window_start, base_task.active_window_end);
  if (first_feat2 == features[1].data_end()) {
    print_status_update(base_task);
    return false;
  }

  base_task.feats1 = absl::MakeConstSpan(&(*first_feat1), &(*last_feat1));
  base_task.feats2 = absl::MakeConstSpan(&(*first_feat2), &(*last_feat2));
  return true;
}

bool Simulation::map_barriers_to_window(TaskPW& base_task, const Chromosome& chrom) {
  const auto& barriers = chrom.barriers();

  // Find all barriers falling within the outer window. Skip over windows with 0 barriers
  auto [first_barrier, last_barrier] =
      barriers.equal_range(base_task.window_start, base_task.window_end);
  if (first_barrier == barriers.data_end()) {
    spdlog::info(FMT_STRING("Skipping {}[{}-{}]..."), base_task.chrom->name(),
                 base_task.active_window_start, base_task.active_window_end);
    return false;
  }

  base_task.barriers = absl::MakeConstSpan(&(*first_barrier), &(*last_barrier));
  return true;
}

absl::flat_hash_set<size_t> Simulation::import_task_filter(
    const boost::filesystem::path& path_to_task_filter) {
  if (path_to_task_filter.empty()) {
    return absl::flat_hash_set<size_t>{};
  }

  compressed_io::Reader r(path_to_task_filter);
  std::string buff;
  size_t i = 0;
  absl::flat_hash_set<size_t> tasks;
  try {
    size_t id;  // NOLINT
    for (; r.getline(buff); ++i) {
      utils::parse_numeric_or_throw(buff, id);
      if (tasks.contains(id)) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Found a duplicate entry for task #{}"), id));
      }
      tasks.insert(id);
    }
  } catch (const std::exception& e) {
    if (absl::StartsWith(e.what(), "Unable to convert field")) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Detected a malformed entry at line {} of file {}: {}"), i,
                      path_to_task_filter, e.what()));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while parsing line {} of file {}: {}"),
                    i, path_to_task_filter, e.what()));
  }
  return tasks;
}

void Simulation::rethrow_exceptions() const {
  assert(!this->ok());                 // NOLINT
  assert(!this->_exceptions.empty());  // NOLINT

  std::string error_msg = "The following error(s) occurred while simulating loop extrusion:";
  for (const auto& exc_ptr : this->_exceptions) {
    try {
      std::rethrow_exception(exc_ptr);
    } catch (const std::exception& e) {
      absl::StrAppendFormat(&error_msg, "\n  %s", e.what());
    } catch (...) {
      absl::StrAppend(&error_msg,
                      "\n  An unhandled exception was caught! This should never happen! If you see "
                      "this message, please file an issue on GitHub.");
    }
  }
  throw std::runtime_error(error_msg);
}

void Simulation::handle_exceptions() {
  assert(!this->ok());  // NOLINT
  spdlog::error("MoDLE encountered an exception. Shutting down worker threads...");
  this->_tpool.wait_for_tasks();  // Wait on simulate_worker threads
  this->rethrow_exceptions();
}

}  // namespace modle
