# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS headers)
find_package(bshoshany-thread-pool CONFIG REQUIRED)
find_package(CLI11 CONFIG REQUIRED)
find_package(cpp-sort CONFIG REQUIRED)
find_package(hictk CONFIG REQUIRED)
find_package(phmap CONFIG REQUIRED)
find_package(readerwriterqueue CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(tomlplusplus CONFIG REQUIRED)

add_executable(modle_tools main.cpp)

target_sources(
  modle_tools
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/annotate_barriers.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/modle_tools/modle_tools_config.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/eval.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/transform.cpp
)

target_link_libraries(
  modle_tools
  PRIVATE
    modle_project_warnings
    modle_project_options
    fmt::fmt-header-only
    Modle::cmatrix
    Modle::common
    Modle::config
    Modle::interval_tree
    Modle::io_bed
    Modle::io_bigwig
    Modle::io_chrom_sizes
    Modle::io_compressed
    Modle::io_contact_matrix
    Modle::io_json
    Modle::io_toml
    Modle::stats
)

target_link_system_libraries(
  modle_tools
  PRIVATE
  absl::span
  absl::strings
  absl::str_format
  absl::time
  Boost::headers
  bshoshany-thread-pool::bshoshany-thread-pool
  CLI11::CLI11
  hictk::libhictk
  cpp-sort::cpp-sort
  phmap
  readerwriterqueue::readerwriterqueue
  tomlplusplus::tomlplusplus
  spdlog::spdlog_header_only
)

target_include_directories(modle_tools PRIVATE include)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties(
    modle_tools
    PROPERTIES
      LINK_FLAGS_RELEASE
        -s
  )
endif()

include(GNUInstallDirs)
install(
  TARGETS
    modle_tools
  COMPONENT Runtime
  LIBRARY
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/scratch/"
  ARCHIVE
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/scratch/"
  RUNTIME
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
  PUBLIC_HEADER
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/scratch/"
  PRIVATE_HEADER
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/scratch/"
)

install(FILES "${PROJECT_SOURCE_DIR}/LICENSE" DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/modle/" COMPONENT Runtime)
