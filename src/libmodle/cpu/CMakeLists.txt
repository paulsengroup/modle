# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Boost CONFIG REQUIRED COMPONENTS headers)
find_package(bshoshany-thread-pool CONFIG REQUIRED)
find_package(concurrentqueue CONFIG REQUIRED)
find_package(cpp-sort CONFIG REQUIRED)
find_package(phmap CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(xxHash CONFIG REQUIRED)

add_library(libmodle_cpu)
add_library(Modle::libmodle_cpu ALIAS libmodle_cpu)
set_target_properties(
  libmodle_cpu
  PROPERTIES
    OUTPUT_NAME
      modle_cpu
)

target_sources(
  libmodle_cpu
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/collision_encoding_impl.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/context_manager_impl.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/register_contacts.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_simulate.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/simulation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/simulation_correct_moves.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/simulation_detect_collisions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/state_logger.cpp
)

target_include_directories(libmodle_cpu PUBLIC include/)

target_link_libraries(
  libmodle_cpu
  PRIVATE
    modle_project_warnings
    modle_project_options
    Modle::interval_tree
    Modle::io_compressed
    Modle::stats
  PUBLIC
    Modle::cmatrix
    Modle::common
    Modle::config
    Modle::libmodle_internal
    Modle::io_bed
    Modle::io_bigwig
    Modle::io_contact_matrix
    concurrentqueue::concurrentqueue
    fmt::fmt-header-only
)

target_link_system_libraries(
  libmodle_cpu
  PRIVATE
  cpp-sort::cpp-sort
  PUBLIC
  Boost::headers
  bshoshany-thread-pool::bshoshany-thread-pool
  phmap
  spdlog::spdlog_header_only
  xxHash::xxhash
)

install(TARGETS libmodle_cpu ARCHIVE INCLUDES LIBRARY)
