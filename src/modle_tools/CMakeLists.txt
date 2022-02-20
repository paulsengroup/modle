# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED)
find_package(CLI11 CONFIG REQUIRED)
find_package(cpp-sort CONFIG REQUIRED)
find_package(readerwriterqueue CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(tomlplusplus CONFIG REQUIRED)

add_executable(modle_tools main.cpp)

target_sources(
  modle_tools
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cli.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/modle_tools/modle_tools_config.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/eval.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/find_barrier_clusters.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/noisify.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/transform.cpp)

target_link_libraries(
  modle_tools
  PRIVATE project_warnings
          project_options
          Modle::cmatrix
          Modle::common
          Modle::interval_tree
          Modle::io_bed
          Modle::io_bigwig
          Modle::io_compressed
          Modle::io_cooler
          Modle::stats)

target_link_system_libraries(
  modle_tools
  PRIVATE
  absl::algorithm_container
  absl::btree
  absl::failure_signal_handler
  absl::flat_hash_map
  absl::flat_hash_set
  absl::span
  absl::strings
  absl::str_format
  absl::symbolize
  absl::time
  absl::variant
  Boost::filesystem
  Boost::headers
  CLI11::CLI11
  cpp-sort::cpp-sort
  fmt::fmt
  readerwriterqueue::readerwriterqueue
  tomlplusplus::tomlplusplus
  spdlog::spdlog
  thread-pool)

target_include_directories(modle_tools PRIVATE include)

if(ENABLE_PCH)
  target_precompile_headers(
    modle_tools
    PRIVATE
    <algorithm>
    <cassert>
    <cstddef>
    <cstdint>
    <filesystem>
    <stdexcept>
    <string>
    <string_view>
    <utility>
    <vector>
    <absl/container/flat_hash_set.h>
    <absl/types/span.h>
    <fmt/format.h>)
endif()

include(GNUInstallDirs)
install(
  TARGETS modle_tools
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${PROJECT_SOURCE_DIR}/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/modle_tools/)
