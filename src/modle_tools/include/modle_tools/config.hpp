#pragma once

#include <absl/container/flat_hash_set.h>

#include <cstdint>
#include <filesystem>
#include <string>
#include <thread>
#include <vector>

namespace modle::tools {

struct config {
  std::filesystem::path path_to_input_matrix;
  std::filesystem::path path_to_reference_matrix;
  std::filesystem::path path_to_chr_subranges;
  std::filesystem::path output_base_name;
  std::filesystem::path chr_sizes;
  std::filesystem::path tmp_dir{std::filesystem::temp_directory_path()};

  bool force{false};
  bool keep_tmp_files{false};
  std::size_t nthreads{std::thread::hardware_concurrency()};

  std::size_t bin_size{0};
  std::size_t diagonal_width;

  // eval
  bool compute_spearman{true};
  bool compute_pearson{true};
  std::size_t sliding_window_size{0};
  std::size_t sliding_window_overlap{0};
  bool deplete_contacts_from_reference{true};

  // stats
  bool dump_depleted_matrices{false};
  std::string output_path_for_histograms{};
  std::vector<std::string> chromosomes_excluded_vect{};
  double depletion_multiplier{1.0};

  // Internal
  absl::flat_hash_set<std::string> chromosomes_excluded{};
};
}  // namespace modle::tools
