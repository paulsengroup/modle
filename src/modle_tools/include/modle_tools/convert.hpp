#pragma once

#include <cstdint>
#include <random>

#include "modle_tools/cli.hpp"

namespace modle::tools {

std::normal_distribution<double> init_noise_generator(uint64_t range);

std::string prepare_files_for_juicer_tools(const modle::tools::config& c);

void convert_to_hic(const modle::tools::config& c, std::string& argv);

void convert_to_tsv(const modle::tools::config& c);

}  // namespace modle::tools