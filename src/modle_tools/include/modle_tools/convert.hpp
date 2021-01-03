#pragma once

#include <boost/asio/thread_pool.hpp>
#include <cstdint>  // for uint64_t
#include <random>   // for normal_distribution
#include <string_view>

namespace modle::tools {

struct config;  // Pre-declaration

inline std::normal_distribution<double> init_noise_generator(uint64_t noise_range,
                                                             double noise_stddev);

inline void convert_to_cooler(const modle::tools::config& c);

inline void convert_to_hic(const modle::tools::config& c, std::string_view argv);

inline void convert_to_tsv(boost::asio::thread_pool& tpool, const modle::tools::config& c);

}  // namespace modle::tools

#include "../../convert_impl.hpp"