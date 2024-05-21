# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS headers)
find_package(fmt CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(xxHash CONFIG REQUIRED)

add_library(libmodle_internal)
add_library(Modle::libmodle_internal ALIAS libmodle_internal)
set_target_properties(libmodle_internal PROPERTIES OUTPUT_NAME modle_internal)

target_sources(
  libmodle_internal
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_barriers.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/genome.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_barriers_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/genome_impl.hpp)

target_link_libraries(
  libmodle_internal
  PRIVATE project_options project_warnings Modle::io_chrom_sizes
  PUBLIC Modle::cmatrix
         Modle::common
         Modle::io_bed
         Modle::interval_tree
         fmt::fmt-header-only)

target_link_system_libraries(
  libmodle_internal
  PRIVATE
  absl::time
  PUBLIC
  absl::btree
  absl::span
  Boost::headers
  spdlog::spdlog_header_only
  xxHash::xxhash)

target_include_directories(libmodle_internal PUBLIC include/)
