#pragma once

#include <boost/asio/thread_pool.hpp>
#include <cstdint>  // for uint64_t
#include <random>   // for normal_distribution
#include <string_view>

namespace modle::tools {

struct config;  // Pre-declaration

std::normal_distribution<double> init_noise_generator(uint64_t range);

std::string prepare_files_for_juicer_tools(boost::asio::thread_pool& tpool, std::size_t i,
                                           const modle::tools::config& c);

void convert_to_hic(boost::asio::thread_pool& tpool, const modle::tools::config& c,
                    std::string_view argv);

void convert_to_tsv(boost::asio::thread_pool& tpool, const modle::tools::config& c);

}  // namespace modle::tools