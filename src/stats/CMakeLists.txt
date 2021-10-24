# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Boost QUIET REQUIRED)
find_package(cpp-sort QUIET REQUIRED)

add_library(modle_stats INTERFACE)
add_library(Modle::stats ALIAS modle_stats)

target_sources(
  modle_stats INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/correlation_impl.hpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/descriptive_impl.hpp ${CMAKE_CURRENT_SOURCE_DIR}/test_impl.hpp)

target_include_directories(modle_stats INTERFACE include/)
target_link_libraries(
  modle_stats
  INTERFACE project_warnings
            project_options
            Boost::headers
            cpp-sort::cpp-sort)