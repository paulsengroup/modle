#pragma once

#include <absl/container/flat_hash_set.h>  // for BitMask, flat_hash_set, raw_hash_set

#include <cstddef>     // IWYU pragma: keep for size_t
#include <cstdint>     // for uint64_t
#include <filesystem>  // for path, temp_directory_path
#include <string>      // for string
#include <thread>      // for thread
#include <vector>      // for vector

#include "modle/bed.hpp"  // for BED::Dialect

namespace modle::tools {

struct config {
  std::filesystem::path path_to_input_matrix;
  std::filesystem::path path_to_output_matrix;
  std::filesystem::path path_to_reference_matrix;
  std::filesystem::path path_to_chrom_subranges;
  std::filesystem::path output_base_name;
  std::filesystem::path chrom_sizes;
  std::filesystem::path tmp_dir{std::filesystem::temp_directory_path()};

  bool force{false};
  bool keep_tmp_files{false};
  size_t nthreads{std::thread::hardware_concurrency()};

  size_t bin_size{0};
  size_t diagonal_width;
  bed::BED::Dialect bed_dialect{bed::BED::Dialect::BED6};
  bool strict_bed_validation{true};

  // eval
  bool compute_spearman{true};
  bool compute_pearson{true};
  bool compute_edist{true};
  size_t sliding_window_size{0};
  size_t sliding_window_overlap{0};
  bool deplete_contacts_from_reference{true};

  // filter barriers
  std::filesystem::path path_to_extrusion_barrier_motifs_bed;
  std::vector<std::filesystem::path> path_to_bed_files_for_filtering;
  std::string filtering_criterion{"intersection"};

  // noisify_contacts
  double genextreme_mu{0};
  double genextreme_sigma{7'500};
  double genextreme_xi{0.001};
  uint64_t seed{0};

  // stats
  bool dump_depleted_matrices{false};
  std::string output_path_for_histograms{};
  std::vector<std::string> chromosomes_excluded_vect{};
  double depletion_multiplier{1.0};

  // Internal
  absl::flat_hash_set<std::string> chromosomes_excluded{};
};
}  // namespace modle::tools
