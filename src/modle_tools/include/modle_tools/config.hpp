#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <thread>
#include <vector>

namespace modle::tools {

struct config {
  std::vector<std::string> path_to_input_matrices;
  std::string output_base_name;
  std::string chr_sizes;
  bool convert_to_hic{true};
  bool convert_to_tsv{true};
  bool compress{true};
  bool add_noise{false};
  bool compute_spearman{true};
  bool compute_pearson{true};
  bool force{false};
  bool keep_tmp_files{false};
  uint64_t sliding_window_size{0};
  uint64_t sliding_window_overlap{0};
  std::string path_to_juicer_tools;
  std::string tmp_dir{std::filesystem::temp_directory_path()};
  int exit_code{0};
  uint64_t seed{0};
  double noise_stdev{0};
  uint64_t noise_range{50'000};
  uint64_t juicer_tools_mem{2 * 1024 * 1024 * 1024ULL};  // 2GiB
  std::string path_to_reference_matrix;
  std::string chr_name_hic{};
  uint64_t chr_offset_hic{UINT64_MAX};
  std::size_t nthreads{std::thread::hardware_concurrency()};
};
}  // namespace modle::tools