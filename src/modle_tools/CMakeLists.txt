# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl REQUIRED)
find_package(Boost REQUIRED)
find_package(CLI11 REQUIRED)
find_package(readerwriterqueue REQUIRED)
find_package(spdlog REQUIRED)

add_executable(modle_tools main.cpp)

target_sources(
  modle_tools
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cli.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/include/modle_tools/config.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/eval.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/filter_barriers.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/find_barrier_clusters.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/noisify.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/stats.cpp)

target_link_libraries(
  modle_tools
  PRIVATE project_warnings
          project_options
          absl::algorithm_container
          absl::btree
          absl::flat_hash_map
          absl::flat_hash_set
          absl::failure_signal_handler
          absl::strings
          absl::symbolize
          absl::time
          Boost::headers
          CLI11::CLI11
          readerwriterqueue::readerwriterqueue
          spdlog::spdlog
          Modle::common
          Modle::correlation
          Modle::libmodle_io)

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
