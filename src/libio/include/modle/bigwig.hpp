#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"

namespace modle::bigwig {
template <typename N>
uint64_t write_range(
    absl::flat_hash_map<std::pair<std::string, N>, const std::vector<double>*>& data,
    uint64_t offset, uint64_t span, uint64_t step, std::string_view output_file);

uint64_t write_range(const std::string& chr_name, uint64_t chr_length,
                     const std::vector<double>& vals, uint64_t offset, uint64_t span, uint64_t step,
                     std::string_view output_file);
}  // namespace modle::bigwig
