# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED)
find_package(CLI11 CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(range-v3 CONFIG REQUIRED)
find_package(xxHash CONFIG REQUIRED)

add_library(modle_common INTERFACE)

add_library(Modle::common ALIAS modle_common)

target_sources(
  modle_common
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/common.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/simulation_config.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/genextreme_value_distribution.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/random.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/random_sampling.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/suppress_compiler_warnings.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/utils.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_impl.hpp)

target_link_libraries(modle_common INTERFACE project_options project_warnings)

target_link_system_libraries(
  modle_common
  INTERFACE
  absl::span
  absl::strings
  bitflags::bitflags
  Boost::headers
  Boost::filesystem
  Boost::random
  CLI11::CLI11
  fmt::fmt
  range-v3::range-v3
  $<$<NOT:$<BOOL:${MODLE_USE_MERSENNE_TWISTER}>>:XoshiroCpp::XoshiroCpp>
  xxHash::xxhash
  mscharconv)

target_include_directories(modle_common INTERFACE include/)
