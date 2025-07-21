# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)

add_library(modle_config)

add_library(Modle::config ALIAS modle_config)

include("${PROJECT_SOURCE_DIR}/cmake/Versioning.cmake")

set(MODLE_VERSION_SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/version")
# gersemi: off
ConfigureVersioning(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${MODLE_VERSION_SOURCE_DIR}"
)
# gersemi: on

add_dependencies(modle_config _modle_check_git)

target_sources(
  modle_config
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include/modle/config/git.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/modle/config/version.hpp"
  PRIVATE
    "${MODLE_VERSION_SOURCE_DIR}/git.cpp"
    "${MODLE_VERSION_SOURCE_DIR}/version.cpp"
)

target_link_libraries(
  modle_config
  PRIVATE
    modle_project_options
    modle_project_warnings
    Modle::common
    fmt::fmt-header-only
)
target_link_system_libraries(modle_config PUBLIC absl::strings)

target_include_directories(modle_config PUBLIC include/)
