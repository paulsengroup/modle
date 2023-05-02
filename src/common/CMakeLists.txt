# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(bitflags CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS headers random)
find_package(CLI11 CONFIG REQUIRED)
find_package(FastFloat CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(range-v3 CONFIG REQUIRED)
find_package(xoshiro-cpp CONFIG REQUIRED)
find_package(xxHash CONFIG REQUIRED)

add_library(modle_common INTERFACE)

add_library(Modle::common ALIAS modle_common)

target_sources(
  modle_common
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/cli_utils.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/common.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/const_map.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/dna.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/fmt_helpers.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/genextreme_value_distribution.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/numeric_utils.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/random.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/random_sampling.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/simulation_config.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/suppress_compiler_warnings.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/utils.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cli_utils_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/const_map_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/numeric_utils_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_impl.hpp)

target_link_libraries(modle_common INTERFACE project_options project_warnings fmt::fmt)

target_link_system_libraries(
  modle_common
  INTERFACE
  absl::span
  absl::strings
  bitflags::bitflags
  Boost::headers
  Boost::random
  CLI11::CLI11
  FastFloat::fast_float
  range-v3::range-v3
  xoshiro-cpp::xoshiro-cpp
  xxHash::xxhash)

target_include_directories(modle_common INTERFACE include/)
