# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Boost QUIET REQUIRED)
find_package(cpp-sort QUIET REQUIRED)

add_library(modle_correlation INTERFACE)
add_library(Modle::correlation ALIAS modle_correlation)

target_sources(cmatrix INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/correlation_impl.hpp)

target_include_directories(modle_correlation INTERFACE include/)
target_link_libraries(modle_correlation INTERFACE project_warnings project_options cpp-sort::cpp-sort)

if(ENABLE_PCH)
  target_precompile_headers(
    modle_correlation
    INTERFACE
    <algorithm>
    <cstddef>
    <cstdint>
    <type_traits>
    <utility>
    <vector>
    <absl/types/span.h>
    <boost/dynamic_bitset/dynamic_bitset.hpp>
    [["modle/correlation.hpp"]])
endif()
