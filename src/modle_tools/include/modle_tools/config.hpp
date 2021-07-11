#pragma once

#include <absl/container/flat_hash_set.h>  // for BitMask, flat_hash_set, raw_hash_set
#include <absl/types/variant.h>

#include <boost/filesystem/operations.hpp>  // for temp_directory_path
#include <boost/filesystem/path.hpp>        // for path
#include <cstddef>                          // IWYU pragma: keep for size_t
#include <cstdint>                          // for uint64_t
#include <string>                           // for string
#include <thread>                           // for thread
#include <vector>                           // for vector

#include "modle/bed.hpp"  // for BED::Dialect

namespace modle::tools {

struct eval_config {  // NOLINT
  // IO
  boost::filesystem::path path_to_input_matrix;
  boost::filesystem::path output_base_name;
  boost::filesystem::path path_to_reference_matrix;
  boost::filesystem::path path_to_chrom_subranges;
  boost::filesystem::path tmp_dir{boost::filesystem::temp_directory_path()};
  bool keep_tmp_files{false};
  bool force{false};

  // Correlation methods
  bool compute_spearman{true};
  bool compute_pearson{true};
  bool compute_edist{true};

  // Reference contacts
  size_t bin_size{0};
  size_t diagonal_width;
  double depletion_multiplier{1.0};
  bool deplete_contacts_from_reference{true};

  // Other
  size_t nthreads{std::thread::hardware_concurrency()};
  size_t sliding_window_size{0};
  size_t sliding_window_overlap{0};
};

struct filter_barrier_config {  // NOLINT
  // IO
  boost::filesystem::path path_to_extrusion_barrier_motifs_bed;
  std::vector<boost::filesystem::path> path_to_bed_files_for_filtering;

  // Other
  std::string filtering_criterion{"intersection"};
  bed::BED::Dialect bed_dialect{bed::BED::Dialect::BED6};
  bool strict_bed_validation{true};
};

struct noisify_config {  // NOLINT
  // IO
  boost::filesystem::path path_to_input_matrix;
  boost::filesystem::path path_to_output_matrix;
  bool force{false};

  // Generic
  size_t diagonal_width;
  size_t bin_size{0};

  // Noise properties
  double genextreme_mu{0};
  double genextreme_sigma{7'500};  // NOLINT
  double genextreme_xi{0.001};     // NOLINT
  uint64_t seed{0};
};

struct stats_config {  // NOLINT
  // IO
  boost::filesystem::path path_to_input_matrix;
  boost::filesystem::path path_to_chrom_subranges;
  std::string output_path_for_histograms{};
  bool dump_depleted_matrices{false};
  bool force{false};

  // Contact matrix
  size_t bin_size{0};
  size_t diagonal_width;
  std::vector<std::string> chromosomes_excluded_vect{};
  double depletion_multiplier{1.0};
};

// clang-format off
using config = absl::variant<absl::monostate,
                             eval_config,
                             filter_barrier_config,
                             noisify_config,
                             stats_config>;
// clang-format on

}  // namespace modle::tools
