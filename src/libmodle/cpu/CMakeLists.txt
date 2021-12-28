# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl QUIET REQUIRED)
find_package(Boost QUIET REQUIRED)
find_package(concurrentqueue QUIET REQUIRED)
find_package(cpp-sort QUIET REQUIRED)
find_package(spdlog QUIET REQUIRED)
find_package(xxHash QUIET REQUIRED)

add_library(libmodle_cpu)
add_library(Modle::libmodle_cpu ALIAS libmodle_cpu)
set_target_properties(libmodle_cpu PROPERTIES OUTPUT_NAME modle_cpu)

target_sources(
  libmodle_cpu
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/collision_encoding_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_correct_moves.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_detect_collisions.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_common.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_simulate.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_perturbate.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_replay.cpp)

target_include_directories(libmodle_cpu PUBLIC include/)

target_link_libraries(
  libmodle_cpu
  PRIVATE project_warnings
          project_options
          absl::btree
          absl::fixed_array
          absl::hash
          absl::strings
          absl::str_format
          absl::time
          cpp-sort::cpp-sort
          spdlog::spdlog
          Modle::interval_tree
          Modle::io_compressed
          Modle::io_cooler
          Modle::stats
  PUBLIC absl::flat_hash_set
         absl::span
         Boost::filesystem
         Boost::headers
         concurrentqueue::concurrentqueue
         fmt::fmt
         xxHash::xxHash
         Modle::cmatrix
         Modle::common
         Modle::libmodle_internal
         Modle::io_bed
         thread-pool)

if(ENABLE_PCH)
  target_precompile_headers(
    libmodle_cpu
    PRIVATE
    <algorithm>
    <cstddef>
    <cstdint>
    <filesystem>
    <limits>
    <memory>
    <string>
    <string_view>
    <utility>
    <vector>
    <absl/types/span.h>
    <fmt/format.h>)
endif()

install(
  TARGETS libmodle_cpu
  ARCHIVE
  INCLUDES
  LIBRARY)
