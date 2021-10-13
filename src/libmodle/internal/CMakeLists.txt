# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl QUIET REQUIRED)
find_package(Boost QUIET REQUIRED)
find_package(fmt QUIET REQUIRED)
find_package(spdlog QUIET REQUIRED)
find_package(xxHash QUIET REQUIRED)

add_library(libmodle_internal)
add_library(Modle::libmodle_internal ALIAS libmodle_internal)
set_target_properties(libmodle_internal PROPERTIES OUTPUT_NAME modle_internal)

target_sources(
  libmodle_internal
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/genome.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/genome_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_barriers.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors_impl.hpp)
target_link_libraries(
  libmodle_internal
  PRIVATE project_options
          project_warnings
          absl::span
          Modle::common
          Modle::cmatrix
  PUBLIC absl::btree
         absl::optional
         Boost::filesystem
         fmt::fmt
         xxHash::xxhash
         spdlog::spdlog
         Modle::libmodle_io)

target_include_directories(libmodle_internal PUBLIC include/)
