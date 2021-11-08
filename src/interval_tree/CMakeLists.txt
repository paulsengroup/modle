# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl QUIET REQUIRED)
find_package(cpp-sort QUIET REQUIRED)

add_library(interval_tree INTERFACE)
add_library(Modle::interval_tree ALIAS interval_tree)

target_sources(interval_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/interval_tree_impl.hpp)

target_include_directories(interval_tree INTERFACE include/)

target_link_libraries(
  interval_tree
  INTERFACE project_options
            project_warnings
            absl::span
            cpp-sort::cpp-sort
            Modle::common)
