# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(cpp-sort CONFIG REQUIRED)

add_library(interval_tree INTERFACE)
add_library(Modle::interval_tree ALIAS interval_tree)

target_sources(interval_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/interval_tree_impl.hpp)

target_include_directories(interval_tree INTERFACE include/)

target_link_libraries(interval_tree INTERFACE modle_project_options modle_project_warnings Modle::common)

target_link_system_libraries(
  interval_tree
  INTERFACE
  absl::span
  cpp-sort::cpp-sort)
