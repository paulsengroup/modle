# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED)
find_package(cpp-sort CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(xxHash CONFIG REQUIRED)

add_library(cmatrix INTERFACE)
add_library(Modle::cmatrix ALIAS cmatrix)

target_sources(
  cmatrix INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/contacts_impl.hpp ${CMAKE_CURRENT_SOURCE_DIR}/contacts_safe_impl.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/contacts_unsafe_impl.hpp)

target_include_directories(cmatrix INTERFACE include/)

target_link_libraries(
  cmatrix
  INTERFACE project_warnings
            project_options
            Modle::common
            Modle::io_compressed
            Modle::stats)

target_link_system_libraries(
  cmatrix
  INTERFACE
  absl::span
  absl::strings
  Boost::filesystem
  Boost::headers
  cpp-sort::cpp-sort
  fmt::fmt
  xxHash::xxHash
  thread-pool)
